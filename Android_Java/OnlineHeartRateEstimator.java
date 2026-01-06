package com.thinkrace.mt4sdk.function;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import android.util.Log;


public class OnlineHeartRateEstimator {

    // --- Updated Parameters for 50Hz Input ---
    private static final int SRATE_ORIGINAL = 50; // 输入采样率
    private static final int SRATE_DOWN = 25;     // 算法内部采样率
    private static final int DOWNSAMPLE_FACTOR = SRATE_ORIGINAL / SRATE_DOWN; // 50/25 = 2

    private static final int FFT_RES = 1024;
    private static final int DELAY_D = 3;

    // Filter Coefficients (User Provided for 50Hz)
    private double[] b_filt;
    private double[] a_filt;

    // State Cache
    private double[] vOld;
    private int[][] pTR_history;
    private double[][] obs_history;
    private double[][] freq_range_history;
    private int processed_frames = 0;

    // Post Process Cache
    private double last_output_bpm = 0;
    private List<Double> bpm_history = new ArrayList<>();

    private List<Double> smoothingWindow = new ArrayList<>();
    private static final int SMOOTH_WINDOW_SIZE = 5; // 5点中值滤波
    private double smoothedBPM = 0;
    public OnlineHeartRateEstimator() {
        initCoefficients();
        reset();
    }

    private void initCoefficients() {
        // 使用您提供的 50Hz 系数
//        this.b_filt = new double[]{0.0015, 0, -0.0062, 0, 0.0093, 0, -0.0062, 0, 0.0015};
//        this.a_filt = new double[]{1.0000, -6.7349, 19.9541, -33.9953, 36.4462, -25.1872, 10.9588, -2.7447, 0.3029};
        this.b_filt = new double[]{0.033559, 0.0, -0.067119, 0.0, 0.033559};
        this.a_filt = new double[]{1.000000, -3.272717, 4.090226, -2.327539, 0.511327};
    }

    public void reset() {
        vOld = null;
        pTR_history = new int[DELAY_D + 5][];
        obs_history = new double[DELAY_D + 5][];
        freq_range_history = new double[DELAY_D + 5][];
        processed_frames = 0;
        last_output_bpm = 0;
        bpm_history.clear();
    }

    /**
     * 核心执行方法
     * @param rawPPG 输入长度应为 8s * 50Hz = 400 点
     */
    public HeartRateResult execute(double[] rawPPG, double[] rawAccX, double[] rawAccY, double[] rawAccZ) {
        processed_frames++;
        // ================= 调试与佩戴检测逻辑 START =================
        // 1. 计算 PPG 信号的峰峰值 (Peak-to-Peak Amplitude)
        double ppgMax = -Double.MAX_VALUE;
        double ppgMin = Double.MAX_VALUE;
        for (double v : rawPPG) {
            if (v > ppgMax) ppgMax = v;
            if (v < ppgMin) ppgMin = v;
        }
        double amplitude = ppgMax - ppgMin;

        // 2. 打印调试日志 (请在 Logcat 中搜索 "AlgoDebug" 查看)
        Log.d("AlgoDebug", "Frame: " + processed_frames + " | Raw Amplitude: " + amplitude);

        // 3. 信号强度阈值判断 (解决未佩戴仍有数值的问题)
        // 这里的 100.0 是一个经验值。
        // 请查看 Logcat，未佩戴时 Amplitude 是多少（通常很小，比如 < 50），
        // 佩戴时是多少（通常 > 1000）。将此阈值设为两者中间。
        if (amplitude < 200.0) {
            Log.w("AlgoDebug", "Signal too weak (Not worn?), returning 0.");
            return new HeartRateResult(0, new double[rawPPG.length]); // 返回 0 和空波形
        }
            // 1. 预处理 (Filter -> Normalize -> Downsample to 25Hz)
        double[] ppg_proc = preprocess(rawPPG);
        double[] accX_proc = preprocess(rawAccX);
        double[] accY_proc = preprocess(rawAccY);
        double[] accZ_proc = preprocess(rawAccZ);

        // 2. FFT
        DSPUtils.Complex[] ppg_fft = DSPUtils.fft(ppg_proc, FFT_RES);
        DSPUtils.Complex[] accX_fft = DSPUtils.fft(accX_proc, FFT_RES);
        DSPUtils.Complex[] accY_fft = DSPUtils.fft(accY_proc, FFT_RES);
        DSPUtils.Complex[] accZ_fft = DSPUtils.fft(accZ_proc, FFT_RES);

        // 3. 计算用于维纳滤波的幅度谱
        // 频率轴计算 (0 - 25Hz)
        // 1024点FFT，每个bin宽 25/1024 Hz
        double[] ppg_mag = new double[FFT_RES];
        double[] accX_mag = new double[FFT_RES];
        double[] accY_mag = new double[FFT_RES];
        double[] accZ_mag = new double[FFT_RES];

        for(int i=0; i<FFT_RES; i++) {
            ppg_mag[i] = ppg_fft[i].abs();
            accX_mag[i] = accX_fft[i].abs();
            accY_mag[i] = accY_fft[i].abs();
            accZ_mag[i] = accZ_fft[i].abs();
        }

        // 4. 维纳滤波 (核心去噪)
        // 返回的是净化后的幅度谱 (Magnitude Spectrum)
        double[] cleaned_mag_full = performWienerFilterFull(ppg_mag, accX_mag, accY_mag, accZ_mag);

        // 5. 【新增】时域信号重构 (IFFT)
        // 使用 净化后的幅度 + 原始相位 重构复数数组
        DSPUtils.Complex[] cleaned_complex = new DSPUtils.Complex[FFT_RES];
        for(int i=0; i<FFT_RES; i++) {
            double phase = ppg_fft[i].phase();
            double mag = cleaned_mag_full[i];
            // Re = Mag * cos(phase), Im = Mag * sin(phase)
            cleaned_complex[i] = new DSPUtils.Complex(mag * Math.cos(phase), mag * Math.sin(phase));
        }

        // 执行 IFFT
        DSPUtils.Complex[] time_domain_complex = DSPUtils.ifft(cleaned_complex);

        // 取实部，截取有效长度 (因为FFT补零到了1024，实际数据只有 ppg_proc.length)
        double[] wiener_time_domain = new double[ppg_proc.length];
        for(int i=0; i<ppg_proc.length; i++) {
            wiener_time_domain[i] = time_domain_complex[i].re;
        }

        // 6. 截取有效频率范围 (1Hz - 3Hz) 用于 Viterbi 心率计算
        int idx_low = (int) Math.floor(1.0 / (25.0 / FFT_RES));
        int idx_high = (int) Math.ceil(3.0 / (25.0 / FFT_RES));
        int len = idx_high - idx_low + 1;

        double[] freq_axis_roi = new double[len];
        double[] obs_roi = new double[len]; // Region of Interest
        // 寻找 ROI 内的最大峰值位置
        int maxPeakIndex = -1;
        double maxPeakVal = -1;
        for (int i = 0; i < len; i++) {
            int original_idx = idx_low + i;
            freq_axis_roi[i] = original_idx * (25.0 / FFT_RES);
            obs_roi[i] = cleaned_mag_full[original_idx];
            if (obs_roi[i] > maxPeakVal) {
                maxPeakVal = obs_roi[i];
                maxPeakIndex = original_idx;
            }
        }
        if (maxPeakIndex == idx_low) {
            Log.w("AlgoDebug", "Peak at lower boundary (Noise at 58.6 BPM), returning 0.");
            return new HeartRateResult(0, wiener_time_domain); // 虽然心率无效，但波形可能还能看
        }

        for (int i = len - 1; i >= 0; i--) {
            int current_idx = idx_low + i;
            double current_freq = current_idx * (25.0 / FFT_RES);

            // 只检查 > 90 BPM (1.5 Hz) 的频率，怀疑它们是倍频
            if (current_freq > 1.5) {
                // 寻找对应的半频索引 (Half Frequency Index)
                int half_idx = (int)Math.round(current_idx / 2.0);

                // 如果半频在我们的搜索范围内
                if (half_idx >= idx_low) {
                    int roi_idx_half = half_idx - idx_low;
                    int roi_idx_curr = i;

                    // 如果半频处有显著能量 (例如是当前高频能量的 40% 以上)
                    // 说明当前的高频可能只是谐波
                    if (obs_roi[roi_idx_half] > 0.4 * obs_roi[roi_idx_curr]) {
                        // 惩罚高频能量 (打折)
                        obs_roi[roi_idx_curr] *= 0.5;
                    }
                }
            }
        }

        // 归一化观测概率 (针对 ROI)
        double maxVal = DSPUtils.max(obs_roi);
        if (maxVal > 1e-9) {
            for (int i = 0; i < len; i++) obs_roi[i] /= maxVal;
        }

        // Viterbi 估算
        double rawBpm = runViterbi(obs_roi, freq_axis_roi);

        // 执行平滑
        double finalBpm = performSmoothing(rawBpm);

        // 返回平滑后的结果
        return new HeartRateResult(finalBpm, wiener_time_domain);
    }

    private double[] preprocess(double[] raw) {
        // 1. Filter
        double[] filtered = DSPUtils.filter(b_filt, a_filt, raw);
        // 2. Normalize (z-score)
        double mean = DSPUtils.mean(filtered);
        double std = DSPUtils.std(filtered, mean);
        if (std == 0) std = 1;
        for (int i = 0; i < filtered.length; i++) {
            filtered[i] = (filtered[i] - mean) * 0.5 / std;
        }
        // 3. Downsample (50 -> 25, Factor = 2)
        return DSPUtils.downsample(filtered, DOWNSAMPLE_FACTOR);
    }

    // 维纳滤波完整版 (不截断频率，用于IFFT)
    private double[] performWienerFilterFull(double[] ppg, double[] accX, double[] accY, double[] accZ) {
        int N = ppg.length;
        double[] ppg_norm = DSPUtils.normalizeMax(ppg);
        double[] accX_norm = DSPUtils.normalizeMax(accX);
        double[] accY_norm = DSPUtils.normalizeMax(accY);
        double[] accZ_norm = DSPUtils.normalizeMax(accZ);

        // WF1 Stage
        double[] WF1 = new double[N];
        for (int i = 0; i < N; i++) {
            double acc_avg = (accX_norm[i] + accY_norm[i] + accZ_norm[i]) / 3.0;
            double val = 1 - (acc_avg / (ppg_norm[i] + 1e-10));
            WF1[i] = val < 0 ? 0 : val; // 为了IFFT稳定性，负值通常置0或者很小的值
        }

        double[] ppg_clean_1 = new double[N];
        for (int i = 0; i < N; i++) ppg_clean_1[i] = ppg[i] * WF1[i]; // abs(PPG) * WF1

        // WF2 Stage
        double[] WF2 = new double[N];
        double[] ppg_clean_1_norm = DSPUtils.normalizeMax(ppg_clean_1);

        for(int i=0; i<N; i++){
            double acc_avg = (accX_norm[i] + accY_norm[i] + accZ_norm[i]) / 3.0;
            double signal = ppg_clean_1_norm[i];
            double noise = acc_avg;
            WF2[i] = signal / (signal + noise + 1e-10);
        }

        double[] ppg_final = new double[N];
        for (int i = 0; i < N; i++) {
            // 这里将两个阶段的结果融合 (对应 Matlab: W1... + W2...)
            // 注意：Matlab代码中 PPG_ave_FFT_FIN = W1_Clean + W2_Clean
            double w1_out = ppg_clean_1[i];
            double w2_out = ppg[i] * WF2[i]; // 或者是 ppg_clean_1 * WF2, 视具体实现
            // 按照 Matlab 逻辑: PPG_ave_FFT_FIN(i,:) = W1... + W2...
            // 但 Matlab 中 W1和W2经过了 std 归一化。这里简化为直接叠加幅度
            ppg_final[i] = w1_out + w2_out;
        }
        return ppg_final;
    }

    private double runViterbi(double[] observation, double[] freq_axis) {
        // (同前次回答，无需修改)
        // 略去以节省篇幅，逻辑完全一致
        // ...

        // 为确保完整性，此处应保留 runViterbi 的完整代码
        // 如果您在IDE中拼接，请使用上一个回答中的 runViterbi 和 postProcess 方法
        // 仅仅需要注意：freq_axis 的范围是 1Hz-3Hz
        int stateNum = observation.length;
        int ptr = (processed_frames - 1) % DELAY_D;

        if (vOld == null) {
            vOld = new double[stateNum];
            for (int i = 0; i < stateNum; i++) vOld[i] = observation[i];
            pTR_history = new int[DELAY_D][stateNum];
            obs_history = new double[DELAY_D][stateNum];
            freq_range_history = new double[DELAY_D][stateNum];
            return -1;
        }

        double[] vNew = new double[stateNum];
        int[] ptrNew = new int[stateNum];

        for (int curr = 0; curr < stateNum; curr++) {
            double maxVal = -Double.MAX_VALUE;
            int bestPrev = 0;
            for (int prev = 0; prev < stateNum; prev++) {
                double dist = Math.abs(curr - prev);
//                double logTrans = -(dist * dist) / (2.0 * 3.0 * 3.0);
                double sigma = 1.0;
                double logTrans = -(dist * dist) / (2.0 * sigma * sigma);
                double score = vOld[prev] + logTrans + Math.log(observation[curr] + 1e-10);
                if (score > maxVal) {
                    maxVal = score;
                    bestPrev = prev;
                }
            }
            vNew[curr] = maxVal;
            ptrNew[curr] = bestPrev;
        }

        double maxV = DSPUtils.max(vNew);
        for(int i=0; i<stateNum; i++) vNew[i] -= maxV;
        vOld = vNew;

        for(int d=0; d<DELAY_D-1; d++) {
            pTR_history[d] = pTR_history[d+1];
            obs_history[d] = obs_history[d+1];
            freq_range_history[d] = freq_range_history[d+1];
        }
        pTR_history[DELAY_D-1] = ptrNew;
        obs_history[DELAY_D-1] = observation;
        freq_range_history[DELAY_D-1] = freq_axis;

        if (processed_frames < DELAY_D) return -1;

        int bestStateIdx = 0;
        double bestScore = -Double.MAX_VALUE;
        for(int i=0; i<stateNum; i++) {
            if (vOld[i] > bestScore) {
                bestScore = vOld[i];
                bestStateIdx = i;
            }
        }

        int tempState = bestStateIdx;
        for (int d = DELAY_D - 1; d >= 0; d--) {
            if (pTR_history[d] != null) tempState = pTR_history[d][tempState];
        }

        double raw_bpm = freq_range_history[0][tempState] * 60.0;
        return postProcess(raw_bpm, freq_range_history[0], obs_history[0]);
    }

    /**
     * 对最终输出进行平滑处理
     */
    private double performSmoothing(double rawBpm) {
        if (rawBpm <= 0) return 0;

        // 1. 中值滤波 (去除突发的 137 这种尖峰)
        smoothingWindow.add(rawBpm);
        if (smoothingWindow.size() > SMOOTH_WINDOW_SIZE) {
            smoothingWindow.remove(0);
        }

        List<Double> sorted = new ArrayList<>(smoothingWindow);
        Collections.sort(sorted);
        double medianBpm;
        if (sorted.size() % 2 == 0) {
            medianBpm = (sorted.get(sorted.size()/2 - 1) + sorted.get(sorted.size()/2)) / 2.0;
        } else {
            medianBpm = sorted.get(sorted.size()/2);
        }

        // 2. 指数平滑 (EMA) - 让数值变化看起来更连贯
        // alpha 越小，平滑度越高，响应越慢。0.3 是一个平衡值。
        if (smoothedBPM == 0) {
            smoothedBPM = medianBpm;
        } else {
            smoothedBPM = smoothedBPM * 0.7 + medianBpm * 0.3;
        }

        return smoothedBPM;
    }

    private double postProcess(double raw_viterbi_bpm, double[] freq_axis, double[] observation) {
        // 如果是第一帧，直接返回
        if (last_output_bpm == 0) {
            last_output_bpm = raw_viterbi_bpm;
            return raw_viterbi_bpm;
        }

        // 仅当 Viterbi 结果极其离谱（比如一帧跳了 20 BPM）才尝试修正
        // 之前的阈值是 10，对于 50Hz 采样有点太敏感了，建议改为 20
        boolean trigger_correction = Math.abs(raw_viterbi_bpm - last_output_bpm) > 20;

        double final_bpm = raw_viterbi_bpm;

        if (trigger_correction) {
            List<Integer> peaks = DSPUtils.findPeaks(observation);
            Collections.sort(peaks, (p1, p2) -> Double.compare(observation[p2], observation[p1]));
            int K = Math.min(3, peaks.size());
            double best_cand = raw_viterbi_bpm;
            double min_dist = Double.MAX_VALUE;
            for (int k = 0; k < K; k++) {
                int idx = peaks.get(k);
                double val = freq_axis[idx] * 60.0;
                double dist = Math.abs(val - last_output_bpm);
                if (dist < min_dist) { min_dist = dist; best_cand = val; }
            }
            if (Math.abs(raw_viterbi_bpm - last_output_bpm) > 10 && min_dist < 10) final_bpm = best_cand;
        }

        // 注意：这里的 bpm_history 平滑逻辑我已经移到了外部的 performSmoothing
        // 所以这里直接更新 last_output_bpm 并返回 final_bpm 即可
        last_output_bpm = final_bpm;
        return final_bpm;
    }
//    private double postProcess(double raw, double[] freq, double[] obs) {
//        // (同前次回答)
//        if (last_output_bpm == 0) { last_output_bpm = raw; return raw; }
//        boolean trigger = Math.abs(raw - last_output_bpm) > 10;
//        double final_bpm = raw;
//        if (trigger) {
//            List<Integer> peaks = DSPUtils.findPeaks(obs);
//            Collections.sort(peaks, (p1, p2) -> Double.compare(obs[p2], obs[p1]));
//            int K = Math.min(3, peaks.size());
//            double best_cand = raw;
//            double min_dist = Double.MAX_VALUE;
//            for (int k = 0; k < K; k++) {
//                int idx = peaks.get(k);
//                double val = freq[idx] * 60.0;
//                double dist = Math.abs(val - last_output_bpm);
//                if (dist < min_dist) { min_dist = dist; best_cand = val; }
//            }
//            if (Math.abs(raw - last_output_bpm) > 10 && min_dist < 10) final_bpm = best_cand;
//        }
//        bpm_history.add(final_bpm);
//        if (bpm_history.size() > 4) bpm_history.remove(0);
//        double s = 0; for(double b : bpm_history) s+=b; s/=bpm_history.size();
//        last_output_bpm = s;
//        return s;
//    }
}

