package com.thinkrace.mt4sdk.function;

public class HeartRateResult {
    /**
     * 估算的心率值 (BPM)
     */
    public double bpm;

    /**
     * 维纳滤波后的时域信号 (25Hz)
     * 这是通过 IFFT 还原回来的波形，去除了运动伪影噪声。
     * 长度 = 窗口时间(8s) * 25Hz = 约 200 点 (实际由FFT/IFFT截取决定)
     */
    public double[] wienerFilteredSignal;

    public HeartRateResult(double bpm, double[] wienerFilteredSignal) {
        this.bpm = bpm;
        this.wienerFilteredSignal = wienerFilteredSignal;
    }
}