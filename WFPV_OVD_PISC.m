clear; %close all;

% Dataset IDs
IDData = {'DATA_01_TYPE01','DATA_02_TYPE02','DATA_03_TYPE02','DATA_04_TYPE02',...
    'DATA_05_TYPE02','DATA_06_TYPE02','DATA_07_TYPE02','DATA_08_TYPE02','DATA_09_TYPE02',...
    'DATA_10_TYPE02','DATA_11_TYPE02','DATA_12_TYPE02','DATA_S04_T01',...
    'TEST_S01_T01', 'TEST_S02_T01', 'TEST_S02_T02', 'TEST_S03_T02', ...
    'TEST_S04_T02', 'TEST_S05_T02', 'TEST_S06_T01', 'TEST_S06_T02',...
    'TEST_S07_T02', 'TEST_S08_T01'};

% Overall parameters
srate_original = 125;    % original sampling rate (Hz)
FFTres = 1024;           % FFT resolution - user tunable parm1
WFlength = 15;           % Wiener filter length - user tunable parm2
allD = size(IDData,2);   % number of recordings

% Filter parameters（4阶Butterworth带通滤波）
CutoffFreqHzHP = 1;      % 对应60 BPM（心率下限）
CutoffFreqHzLP = 3;      % 对应180 BPM（心率上限）
[b,a] = butter(4, [0.4 4]/(srate_original/2),'bandpass');

% metrics
fullBPM=[];              % 在线维特比估计HR序列（汇总）
fullBPM0=[];             % 真实HR序列（汇总）
myErrorN = zeros(1,allD);% 平均绝对误差（AAE）
myErrorStd = zeros(1,allD);% 误差标准差
myRelError = zeros(1,allD);% 相对误差

% Framework rule（8s窗口，2s步长）
window_sec = 8;          % 窗口时长（秒）
step_sec = 2;            % 步长（秒）
window_samples = window_sec * srate_original;  % 窗口采样点数（原始125Hz）
step_samples = step_sec * srate_original;      % 步长采样点数

% 在线维特比配置
online_params.D = 3;     % 延迟窗口数（3个窗口=6秒，可调整）
online_params.srate_down = 25;  % 下采样后采样率（固定25Hz）
tic;
% -------------------------- 循环处理每个数据集 --------------------------
for idnb = 1 : allD
    temp_srate = 25; % 125 / 5
    temp_FreqRange_full = linspace(0, temp_srate, FFTres);
    [~, temp_lowR] = min(abs(temp_FreqRange_full - CutoffFreqHzHP));
    [~, temp_highR] = min(abs(temp_FreqRange_full - CutoffFreqHzLP));
    base_FreqRange = temp_FreqRange_full(temp_lowR:temp_highR);
    online_params.state_num = length(base_FreqRange);

    % 1.3 生成转移概率矩阵
    trans_full = trans_prob(setdiff(1:23, idnb), base_FreqRange);  % 留一法生成转移矩阵
    trans_full = ((moving(trans_full',4))');  % 4点滑动平均平滑转移矩阵
    % 验证转移矩阵维度（避免后续错误）
    if size(trans_full,1) ~= online_params.state_num || size(trans_full,2) ~= online_params.state_num
        error(['数据集' IDData{idnb} '：trans_full维度错误！请检查基准FreqRange生成逻辑']);
    end

    if idnb > 13
        gt_data = load(['all_data/BPM_data/True' IDData{idnb}(5:end)]);
        current_BPM0 = gt_data.BPM0;
    else
        gt_data = load(['all_data/BPM_data/' IDData{idnb} '_BPMtrace']);
        current_BPM0 = gt_data.BPM0;
    end


    % -------------------------- 重新加载完整数据（避免临时变量覆盖） --------------------------
    clear sig temp_data temp_PPG1 temp_PPG2 temp_PPG_ave temp_PPG_FFT temp_FreqRange_full;
    load(['all_data/PPG_data/' IDData{idnb}]);
        % 确定通道索引
    if idnb>13
        ch1 = 1; ch2 = 2; ch3 = 3; ch4 = 4; ch5 = 5;
    else
        ch1 = 2; ch2 = 3; ch3 = 4; ch4 = 5; ch5 = 6;
    end

    % 计算总窗口数
    windowNb = floor((length(sig) - window_samples)/step_samples) + 1;

    % -------------------------- 变量初始化 --------------------------
    BPM_est = zeros(1, windowNb);  % 原始HR估计（无维特比）
    BPM_est_online = zeros(1, windowNb - online_params.D + 1);  % 在线维特比输出（带延迟）
    clear W1_FFTi W11_FFTi W2_FFTi W21_FFTi W1_PPG_ave_FFT_Clean W2_PPG_ave_FFT_Clean ...
          W11_PPG_ave_FFT_Clean W21_PPG_ave_FFT_Clean PPG_ave_FFT_FIN FreqRangePPG;
    
    % 在线维特比缓存（仅缓存最近D个窗口的关键数据）
    online_cache.obs = [];          % 观测向量缓存：[state_num, D]（每列一个窗口）
    online_cache.freq_range = [];   % 频率范围缓存：[state_num, D]
    online_cache.vOld = [];         % 前向得分缓存：[state_num, 1]（上一窗口的最优得分）
    online_cache.pTR_history = [];  % 路径记录缓存：[state_num, D]（记录每个窗口的转移路径）
    log_trans_full = log(trans_full);
    % -------------------------- 逐窗口处理（核心：在线维特比集成） --------------------------
    for i = 1 : windowNb
        % 1. 截取当前窗口数据
        curSegment = (i-1)*step_samples + 1 : (i-1)*step_samples + window_samples;
        curData = sig(:, curSegment);

        % 2. 信号分离（PPG+加速度）
        PPG1 = curData(ch1,:); PPG2 = curData(ch2,:);
        ACC_X = curData(ch3,:); ACC_Y = curData(ch4,:); ACC_Z = curData(ch5,:);

        % 3. 带通滤波（4阶Butterworth 1~3Hz）
        PPG1 = filter(b,a,PPG1); PPG2 = filter(b,a,PPG2);
        ACC_X = filter(b,a,ACC_X); ACC_Y = filter(b,a,ACC_Y); ACC_Z = filter(b,a,ACC_Z);

        % 4. PPG信号归一化与平均
        PPG_ave = 0.5*(PPG1 - mean(PPG1))/std(PPG1) + ...
                  0.5*(PPG2 - mean(PPG2))/std(PPG2);

        % 5. 下采样到25Hz
        PPG_ave = downsample(PPG_ave, 5);
        ACC_X = downsample(ACC_X, 5); ACC_Y = downsample(ACC_Y, 5); ACC_Z = downsample(ACC_Z, 5);

        % 6. FFT与频率范围截取（1~3Hz）
        PPG_ave_FFT = fft(PPG_ave, FFTres);
        FreqRange_full = linspace(0, online_params.srate_down, size(PPG_ave_FFT,2));
        % 截取有效频率范围（1~3Hz）
        [~, lowR] = min(abs(FreqRange_full - CutoffFreqHzHP));
        [~, highR] = min(abs(FreqRange_full - CutoffFreqHzLP));
        FreqRange = FreqRange_full(lowR:highR);
        PPG_ave_FFT = PPG_ave_FFT(lowR:highR);
        % 加速度FFT（用于维纳滤波）
        ACC_X_FFT = fft(ACC_X, FFTres); ACC_X_FFT = ACC_X_FFT(lowR:highR);
        ACC_Y_FFT = fft(ACC_Y, FFTres); ACC_Y_FFT = ACC_Y_FFT(lowR:highR);
        ACC_Z_FFT = fft(ACC_Z, FFTres); ACC_Z_FFT = ACC_Z_FFT(lowR:highR);

        % 7. 相位声码器（细化频率估计）
        FreqRangePPG(i,:) = FreqRange;
        if i > 1  % 需前一窗口的相位信息
            for ii = 1 : size(FreqRangePPG,2)
                curPhase = angle(PPG_ave_FFT(ii));
                prevPhase = angle(PPG_ave_FFTpr(ii));
                vocoder = zeros(1,20);  % 20个候选频率
                for n = 1:20
                    vocoder(n) = ((curPhase - prevPhase) + 2*pi*(n-1))/(2*pi*step_sec);
                end
                % 选择最接近当前频率的候选值
                [~, deltaidx] = min(abs(vocoder - FreqRange(ii)));
                FreqRangePPG(i,ii) = vocoder(deltaidx);
            end
        end
        % 相位声码器结果平滑（3点滑动平均）
        FreqRangePPG(i,:) = moving(FreqRangePPG(i,:), 3);
        % 保存当前相位用于下一窗口
        PPG_ave_FFTpr = PPG_ave_FFT;

        % 8. 双维纳滤波
        WC1 = WFlength; WC2 = WFlength;  % 维纳滤波窗口长度（15）

        % 8.1 Wiener 1（abs归一化）
        W1_FFTi(i,:) = abs(PPG_ave_FFT)/max(abs(PPG_ave_FFT));
        if i == 1
            W1_PPG_ave_FFT_ALL = W1_FFTi(i,:);
        else
            W1_PPG_ave_FFT_ALL = mean(W1_FFTi(max(1,i-WC1):i,:), 1);
        end
        W1_PPG_ave_FFT_ALL_norm = W1_PPG_ave_FFT_ALL/max(W1_PPG_ave_FFT_ALL);
        W1_ACC_X_FFT_norm = abs(ACC_X_FFT)/max(abs(ACC_X_FFT));
        W1_ACC_Y_FFT_norm = abs(ACC_Y_FFT)/max(abs(ACC_Y_FFT));
        W1_ACC_Z_FFT_norm = abs(ACC_Z_FFT)/max(abs(ACC_Z_FFT));
        WF1 = 1 - ((W1_ACC_X_FFT_norm + W1_ACC_Y_FFT_norm + W1_ACC_Z_FFT_norm)/3) ./ (W1_PPG_ave_FFT_ALL_norm);
        WF1(WF1 < 0) = -1;  % 限制负权重
        W1_PPG_ave_FFT_Clean(i,:) = abs(PPG_ave_FFT) .* WF1;

        % 8.2 Wiener 2（abs归一化）
        
        W2_FFTi(i,:) = abs(PPG_ave_FFT)/max(abs(PPG_ave_FFT));
        if i == 1
            W2_PPG_ave_FFT_ALL = W2_FFTi(i,:);
        else
            W2_PPG_ave_FFT_ALL = mean(W2_FFTi(max(1,i-WC2):i,:), 1);
        end
        W2_PPG_ave_FFT_ALL_norm = W2_PPG_ave_FFT_ALL/max(W2_PPG_ave_FFT_ALL);
        W2_ACC_X_FFT_norm = abs(ACC_X_FFT)/max(abs(ACC_X_FFT));
        W2_ACC_Y_FFT_norm = abs(ACC_Y_FFT)/max(abs(ACC_Y_FFT));
        W2_ACC_Z_FFT_norm = abs(ACC_Z_FFT)/max(abs(ACC_Z_FFT));
        WF2 = W2_PPG_ave_FFT_ALL_norm ./ [(W2_ACC_X_FFT_norm + W2_ACC_Y_FFT_norm + W2_ACC_Z_FFT_norm)/3 + W2_PPG_ave_FFT_ALL_norm];
        W2_PPG_ave_FFT_Clean(i,:) = abs(PPG_ave_FFT) .* WF2;

        % 8.3 维纳滤波结果归一化与融合
        W1_PPG_ave_FFT_Clean(i,:) = W1_PPG_ave_FFT_Clean(i,:)/std(W1_PPG_ave_FFT_Clean(i,:));
        W2_PPG_ave_FFT_Clean(i,:) = W2_PPG_ave_FFT_Clean(i,:)/std(W2_PPG_ave_FFT_Clean(i,:));
        PPG_ave_FFT_FIN(i,:) = W1_PPG_ave_FFT_Clean(i,:) + W2_PPG_ave_FFT_Clean(i,:);

        % ===================== 调试绘图探针 START =====================
        % 仅针对 Dataset 10，且仅绘制出问题的帧范围 (例如第 20 到 35 帧)
        % 你可以根据图表上出错的时间点修改这里的范围
        % if idnb == 10 && i >= 29 && i <= 37
        %     figure(100 + i); clf; set(gcf, 'Position', [100, 100, 1000, 800]);
        % 
        %     % 准备 X 轴 (BPM)
        %     freq_bpm_axis = FreqRange * 60;
        % 
        %     % 获取当前帧的真实心率 (Ground Truth)
        %     if i <= length(current_BPM0)
        %         true_bpm = current_BPM0(i);
        %     else
        %         true_bpm = 0;
        %     end
        % 
        %     % --- 子图 1: 加速度频谱 (噪声源) ---
        %     subplot(3,1,1);
        %     % 归一化并叠加三轴加速度
        %     acc_spec = (abs(ACC_X_FFT) + abs(ACC_Y_FFT) + abs(ACC_Z_FFT))/3;
        %     plot(freq_bpm_axis, acc_spec / max(acc_spec), 'k', 'LineWidth', 1.5);
        %     xline(true_bpm, '--r', 'Ground Truth', 'LineWidth', 2);
        %     title(['Frame ', num2str(i), ' - Accelerometer Noise']);
        %     xlabel('BPM'); ylabel('Norm. Mag'); grid on;
        %     xlim([40 200]);
        % 
        %     % --- 子图 2: 原始 PPG 频谱 (信号+噪声) ---
        %     subplot(3,1,2);
        %     raw_ppg_spec = abs(PPG_ave_FFT);
        %     plot(freq_bpm_axis, raw_ppg_spec / max(raw_ppg_spec), 'b', 'LineWidth', 1.5);
        %     xline(true_bpm, '--r', 'Ground Truth', 'LineWidth', 2);
        %     title(['Frame ', num2str(i), ' - Raw PPG Spectrum']);
        %     xlabel('BPM'); ylabel('Norm. Mag'); grid on;
        %     xlim([40 200]);
        % 
        %     % --- 子图 3: Wiener 滤波后的最终 PPG 频谱 (Viterbi 的输入) ---
        %     subplot(3,1,3);
        %     % 这是给 Viterbi 的实际输入
        %     clean_spec = PPG_ave_FFT_FIN(i, :); 
        %     % 归一化方便观察
        %     clean_spec = clean_spec - min(clean_spec); 
        %     clean_spec = clean_spec / max(clean_spec);
        % 
        %     plot(freq_bpm_axis, clean_spec, 'm', 'LineWidth', 2, 'DisplayName', 'Cleaned PPG (Wiener)');
        %     hold on;
        %     % 标出最大值位置（即如果不加 Viterbi，单纯找峰值会选哪里）
        %     [~, max_idx] = max(clean_spec);
        %     est_peak = freq_bpm_axis(max_idx);
        %     % plot(est_peak, 1, 'vk', 'MarkerFaceColor', 'k');
        %     % text(est_peak, 1.1, sprintf('Max: %.1f', est_peak));
        %     [pks_cl, locs_cl] = findpeaks(clean_spec, freq_bpm_axis, 'SortStr', 'descend', 'NPeaks', 1);
        %     if ~isempty(pks_cl)
        %         % 圈出主峰 (这就是真实的 HR)
        %         plot(locs_cl(1), pks_cl(1), 'bo', 'MarkerSize', 10, 'LineWidth', 2,'HandleVisibility','on','DisplayName', 'Candidate Peak');
        %         % text(locs_cl(1), pks_cl(1)+0.1, 'True HR', 'Color', 'b', 'FontSize', 10, 'HorizontalAlignment', 'center');
        %     end
        %     xline(true_bpm, '--r', 'Ground Truth', 'LineWidth', 2, 'DisplayName', 'Ground Truth');
        %     title(['Frame ', num2str(i), ' - Cleaned PPG (Input to Viterbi)']);
        %     xlabel('BPM'); ylabel('Likelihood'); grid on;
        %     xlim([40 200]);
        %     legend('Location', 'best');
        %     drawnow;
        % end
        % ===================== 调试绘图探针 END =====================



        % 8.4 原始HR估计（无维特比，用于对比）
        [~, max_idx] = max(PPG_ave_FFT_FIN(i,:));
        BPM_est(i) = 60 * FreqRangePPG(i, max_idx);

        % -------------------------- 9. 在线维特比核心逻辑--------------------------
        % 9.1 缓存当前窗口的观测向量与频率范围
        current_obs = PPG_ave_FFT_FIN(i,:)';  % 观测向量：[state_num, 1]
        online_cache.obs = [online_cache.obs, current_obs];  % 追加到缓存
        online_cache.freq_range = [online_cache.freq_range, FreqRangePPG(i,:)'];  % 频率范围缓存

        % 9.2 初始化前向得分（第一个窗口）
        if isempty(online_cache.vOld)
            logTR = log_trans_full;  % 转移概率转对数（避免数值下溢）
            logE = current_obs;        % 观测概率
            % 初始得分：对角线转移概率 * 观测概率
            online_cache.vOld = diag(logTR) .* logE;
            online_cache.pTR_history = zeros(online_params.state_num, 1);  % 初始化路径记录
            continue;
        end
        
        % 9.3 前向更新（逐状态计算最优得分与转移路径）
        logTR = log_trans_full;
        logE = current_obs;
        v = repmat(-inf, online_params.state_num, 1);  % 当前窗口得分
        pTR = zeros(online_params.state_num, 1);       % 当前窗口路径记录

        for state = 1 : online_params.state_num
            bestVal = -inf;
            bestPTR = 0;
            % 遍历所有前一状态，找最优转移路径
            for inner = 1 : online_params.state_num
                val = online_cache.vOld(inner) + logTR(inner, state);
                if val > bestVal
                    bestVal = val;
                    bestPTR = inner;
                end
            end
            pTR(state) = bestPTR;          % 记录当前状态的最优前向状态
            v(state) = logE(state) + bestVal;  % 更新当前状态得分
        end
        % 9.4 缓存前向得分与路径
        % 归一化 vOld，防止数值溢出（这对 Viterbi 路径选择无影响，因为是加性常数）
        v = v - max(v);
        online_cache.vOld = v;
        online_cache.pTR_history = [online_cache.pTR_history, pTR];  % 追加路径记录

        % 9.5 延迟输出（缓存满D个窗口时，回溯并输出最早窗口结果）
        if size(online_cache.obs, 2) == online_params.D
            % 9.5.1 回溯最优路径（从当前窗口反向找D步）
            [~, finalState] = max(online_cache.vOld);  % 当前窗口最优状态
            backtrack_path = zeros(1, online_params.D);  % 回溯路径
            backtrack_path(online_params.D) = finalState;

            % 反向回溯D步
            temp_state = finalState;
            for count = online_params.D - 1 : -1 : 1
                temp_state = online_cache.pTR_history(temp_state, count + 1);  % 取历史路径
                backtrack_path(count) = temp_state;
            end

            % 9.5.2 输出最早窗口的HR（第1个缓存窗口）
            output_idx = i - online_params.D + 1;  % 在线输出的索引
            state_idx = backtrack_path(1);         % 最早窗口的最优状态
            BPM_est_online(output_idx) = 60 * online_cache.freq_range(state_idx, 1);
            % =================================================================
            % 9.5.2 输出最早窗口的HR（含 首帧锚定 + Top-3 修正）
            % =================================================================
            output_idx = i - online_params.D + 1;  % 在线输出的索引
            state_idx = backtrack_path(1);         % Viterbi 回溯出的“原始”最优状态

            % 获取当前输出帧的频率轴
            current_freq_axis = online_cache.freq_range(:, 1);

            % 获取当前输出帧的原始频谱 (Clean PPG)，对应缓存中的第 1 列
            current_obs_frame = online_cache.obs(:, 1);

            % --- [修改点] 核心逻辑 ---
            if output_idx == 1 || output_idx == 2 || output_idx == 3
                % 【策略：首帧锚定】
                % 第1帧直接信任最大峰值，防止后续帧的强噪声通过回溯机制“篡改历史”。
                % 这为后续的 Top-3 修正提供了一个正确的“上一帧基准”。
                [~, max_idx_1st] = max(current_obs_frame);
                final_bpm = 60 * current_freq_axis(max_idx_1st);

                % (可选) 打印日志
                % fprintf('Frame 1 Anchored: %.1f BPM\n', final_bpm);

            else
                   
                % 获取当前 Viterbi 估计值 (原始输出，可能出错)
                raw_viterbi_bpm = 60 * current_freq_axis(state_idx);
                % 获取前一帧 Viterbi 估计值
                prev_viterbi_bpm = BPM_est_online(output_idx - 1); 
                % 获取当前帧的观测谱（已维纳滤波和指数增强的幅度谱）
                freq_range_bpm = 60 * current_freq_axis; % 频率/BPM值向量
                
                % 1. 寻找局部极大值（peaks）
                [pks, locs] = findpeaks(current_obs_frame);
    
                if ~isempty(pks)
                    
                    % 2. 筛选和排序 Top 3 峰值
                    [sorted_pks, sorted_idx] = sort(pks, 'descend');
                    sorted_locs = locs(sorted_idx);
            
                    K = 3; % 选取前 K 个峰值
                    candidates_count = min(K, length(sorted_pks));
            
                    % 提取 Top K 候选 BPM 和强度
                    candidates_bpm = freq_range_bpm(sorted_locs(1:candidates_count));
                    all_strengths = sorted_pks(1:candidates_count);
            
                    % 3. 计算得分（Score Calculation）
                    % 3.1. 距离（Distance）
                    if abs(raw_viterbi_bpm - prev_viterbi_bpm) < 10
                        dists = abs(candidates_bpm - raw_viterbi_bpm);
                    else
                        dists = abs(candidates_bpm - prev_viterbi_bpm);
                    end
                    dists = abs(candidates_bpm - raw_viterbi_bpm);
                    % 3.2. 距离得分 (Dist_Score): 距离越近，得分越高。
                    % 假设距离 15 BPM 以外的得分接近于 0。
                    max_relevant_dist = 15; 
                    dist_scores = max(0, 1 - dists / max_relevant_dist); % 简单线性衰减到 0
            
                    % 3.3. 峰值得分 (Peak_Score): 峰值强度越高，得分越高。
                    peak_scores = all_strengths / max(all_strengths); % 归一化到 [0, 1]
            
                    % 3.4. 综合得分 (Combined_Score): 距离得分和峰值得分的加权平均
                    dist_weight = 0.6; % 权重可调
                    peak_weight = 0.4; % 权重可调
                    combined_scores = dist_weight * dist_scores + peak_weight * peak_scores;
            
                    % 4. 选出最佳修正候选和第二佳修正候选
                    
                    % [修改点 1: 对综合得分进行降序排序，提取索引]
                    [sorted_scores, sorted_indices] = sort(combined_scores, 'descend');
                    
                    % 确保至少有 1 个候选峰
                    if candidates_count >= 1 
                        
                        % 4.1 提取最佳修正候选 (Top 1)
                        best_score = sorted_scores(1);
                        best_idx = sorted_indices(1);
                        best_bpm = candidates_bpm(best_idx);
                        min_dist = dists(best_idx); % 最佳候选峰与Viterbi估计的距离
                        
                        % 4.2 提取第二佳修正候选 (Top 2)
                        second_best_score = 0; % 默认值
                        second_best_bpm = raw_viterbi_bpm; % 默认值
                        second_best_dist = 0; % 默认值
                        third_best_score = 0; % 默认值
                        third_best_bpm = raw_viterbi_bpm; % 默认值
                        third_best_dist = 0; % 默认值
                        if candidates_count >= 2 % 确保有第二佳候选
                            second_best_score = sorted_scores(2);
                            second_best_idx = sorted_indices(2);
                            second_best_bpm = candidates_bpm(second_best_idx);
                            second_best_dist = dists(second_best_idx);
                            third_best_score = sorted_scores(3);
                            third_best_idx = sorted_indices(3);
                            third_best_bpm = candidates_bpm(third_best_idx);
                            third_best_dist = dists(third_best_idx);
                        end
                        
                    else
                        % 如果没有找到任何峰值，设置安全默认值
                        best_score = 0;
                        best_bpm = raw_viterbi_bpm;
                        min_dist = 0;
                        second_best_score = 0;
                        second_best_bpm = raw_viterbi_bpm;
                        second_best_dist = 0;
                    end
                    % -------------------------- 7. 修正条件判断 (增强逻辑) -------------------------- 
                    
                    final_bpm = raw_viterbi_bpm; % 默认不修正
                    
                    if output_idx > 1 % 只有第二帧及以后才能判断跳变
                        
                        % Viterbi 与前一帧的差异（跳变程度）
                        viterbi_jump = abs(raw_viterbi_bpm - prev_viterbi_bpm);
                        best_jump = abs(best_bpm - prev_viterbi_bpm);
                        second_best_jump = abs(second_best_bpm - prev_viterbi_bpm);
                        third_best_jump = abs(third_best_bpm - prev_viterbi_bpm);
                        % 参数定义 (用于新增的置信度拉动逻辑)
                        confidence_threshold = 0.65; % 修正置信度门限
                        absolute_offset_threshold = 10; % 绝对偏移门限 (BPM)
                        smooth_follow_threshold = 5; % 平滑跟随门限 (BPM)
                        
                        is_high_confidence = (best_score > confidence_threshold);
                        is_misaligned = (abs(best_bpm - raw_viterbi_bpm) > absolute_offset_threshold);
                        is_smooth_following = (viterbi_jump < smooth_follow_threshold);
                        
                        % 触发条件 1：Viterbi 发生大的跳变 (>12 BPM)，且修正峰得分够高 (>0.4)
                        if viterbi_jump > 12 && best_jump < 10
                             final_bpm = best_bpm;
                             % fprintf('Frame %d: Viterbi %.1f -> Corrected %.1f (跳变修正)\n', output_idx, raw_viterbi_bpm, final_bpm);
                        elseif viterbi_jump > 12 && best_jump > 10 && second_best_jump < 10 
                             final_bpm = second_best_bpm;
                        elseif is_high_confidence && is_misaligned && is_smooth_following
                             final_bpm = best_bpm;
                        end
                    end
                    
                end
                if abs(final_bpm - prev_viterbi_bpm) > 50
                   % fprintf("%d data %d frame is wrong ,diff bpm is %2f\n",idnb,output_idx,abs(final_bpm - prev_viterbi_bpm));
                   final_bpm = mean(BPM_est_online(output_idx-3:output_idx-1));
                end
                % final_bpm = raw_viterbi_bpm;
            end
            
            % 赋值最终结果
            BPM_est_online(output_idx) = final_bpm;
            


            % 9.5.3 滑动缓存（移除最早的1个窗口，保留最近D-1个）
            online_cache.obs = online_cache.obs(:, 2:end);
            online_cache.freq_range = online_cache.freq_range(:, 2:end);
            online_cache.pTR_history = online_cache.pTR_history(:, 2:end);
        end
        % -------------------------- 在线维特比逻辑结束 --------------------------
    end

    % -------------------------- 10. 加载Ground Truth并计算误差（与离线版一致） --------------------------
    if idnb > 13
        load(['all_data/BPM_data/True' IDData{idnb}(5:end)]);  % 测试集GT加载
    else
        load(['all_data/BPM_data/' IDData{idnb} '_BPMtrace']);  % 训练集GT加载
    end
    
    % 10.1 滑动平均平滑在线结果（与离线版一致，4点平均）
    BPM_est_online_smoothed = moving(BPM_est_online, 4);
    valid_len = min(length(BPM0), length(BPM_est_online_smoothed));
    BPM0_valid = BPM0(1:valid_len);
    BPM_est_online_valid = BPM_est_online_smoothed(1:valid_len);
    BPM_est_valid = BPM_est(1:valid_len);
    BPM_est_valid_smoothed = moving(BPM_est_valid,4);
    % BPM_est_online_valid = BPM_est_valid_smoothed;
    % fprintf("data %d valid len = %d\n",idnb,valid_len);
    % 10.2 计算误差指标
    myErrorN(idnb) = mean(abs(BPM0_valid - BPM_est_online_valid));  % AAE
    myRelError(idnb) = mean(abs(BPM0_valid - BPM_est_online_valid) ./ BPM0_valid) * 100;  % 相对误差（%）
    myErrorStd(idnb) = std(abs(BPM0_valid - BPM_est_online_valid));  % 误差标准差

    % 10.3 汇总结果（用于后续绘图）
    fullBPM = [fullBPM BPM_est_online_valid'];
    fullBPM0 = [fullBPM0 BPM0_valid'];

    % 10.4 保存在线结果（新建Results文件夹）
    if ~exist('Results', 'dir')
        mkdir('Results');
    end
    save(['Results/Result_' IDData{idnb} '_DATA_WFPV_OVD_PISC'], ...
         'BPM_est_online', 'BPM0', 'myErrorN', 'myErrorStd', 'myRelError');
    % if idnb==10,           figure; plot(BPM0_valid,'ro');hold on; plot(BPM_est_online_valid,'blue'); title(['Recording ' num2str(idnb)]); xlabel('Time (in frames)'); ylabel('HR(BPM)'); legend({'Ground truth', 'Estimates'}); end;
    % 打印当前数据集结果
    fprintf('数据集%s：AAE=%.2f BPM，sdAE=%.2f BPM，AAEP=%.2f%%\n', ...
            IDData{idnb}, myErrorN(idnb), myErrorStd(idnb), myRelError(idnb));
end
exe_time = toc;
% -------------------------- 11. 整体结果统计与绘图 --------------------------
% 11.1 打印整体误差
fprintf('\n==================== 整体结果统计 ====================\n');
fprintf('训练集（1-12）：AAE=%.2f(%.2f) BPM，AAEP=%.2f%%\n', mean(myErrorN(1:12)), mean(myErrorStd(1:12)), mean(myRelError(1:12)));
fprintf('测试集（13-23）：AAE=%.2f(%.2f) BPM，AAEP=%.2f%%\n', mean(myErrorN(13:end)), mean(myErrorStd(13:end)), mean(myRelError(13:end)));
fprintf('全数据集：平均AAE=%.2f(%.2f) BPM，AAEP=%.2f%%\n', ...
        mean(myErrorN), mean(myErrorStd), mean(myRelError));
fprintf('运行时间为 %2f s\n', exe_time);
% 11.2 Bland-Altman图（对比GT与在线估计）
BlandAltman(fullBPM0', fullBPM', {'Ground Truth HR (BPM)', 'Estimated HR (BPM)'});
title('在线维特比解码HR估计：Bland-Altman图', 'FontSize', 12);

% % 11.3 相关性分析
% % figure;
% % scatter(fullBPM0, fullBPM, 30, 'b', 'filled');
% % hold on;
% % % 线性拟合
% % p = polyfit(fullBPM0, fullBPM, 1);
% % fit_line = polyval(p, fullBPM0);
% % plot(fullBPM0, fit_line, 'r-', 'LineWidth', 1.5);
% % % 对角线（理想情况）
% % xlim([min(fullBPM0), max(fullBPM0)]);
% % ylim([min(fullBPM), max(fullBPM)]);
% % plot(xlim(), ylim(), 'k--', 'LineWidth', 1);
% % xlabel('Ground Truth HR (BPM)', 'FontSize', 11);
% % ylabel('Online Viterbi Estimated HR (BPM)', 'FontSize', 11);
% % title(sprintf('HR估计相关性：R²=%.4f，相关系数=%.3f', ...
% %         corrcoef(fullBPM0, fullBPM)(1,2)^2, corrcoef(fullBPM0, fullBPM)(1,2)), 'FontSize', 12);
% % legend('估计值', sprintf('拟合线：y=%.3fx+%.3f', p(1), p(2)), '理想线', 'Location', 'best');
% % grid on;
% 
% % 11.4 GT心率分布直方图
% figure;
% histogram(fullBPM0, 20);
% xlabel('Ground Truth HR (BPM)', 'FontSize', 11);
% ylabel('频次', 'FontSize', 11);
% title('Ground Truth HR分布直方图', 'FontSize', 12);
% grid on;