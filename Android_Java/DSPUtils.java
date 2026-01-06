package com.thinkrace.mt4sdk.function;

import java.util.ArrayList;
import java.util.List;

public class DSPUtils {

    // ... (mean, std, filter, downsample, findPeaks, max, normalizeMax 保持不变) ...

    public static double[] filter(double[] b, double[] a, double[] x) {
        double[] y = new double[x.length];
        for (int n = 0; n < x.length; n++) {
            double val = 0;
            for (int k = 0; k < b.length; k++) if (n - k >= 0) val += b[k] * x[n - k];
            for (int k = 1; k < a.length; k++) if (n - k >= 0) val -= a[k] * y[n - k];
            y[n] = val;
        }
        return y;
    }

    public static double mean(double[] data) {
        double s=0; for(double d:data)s+=d; return s/data.length;
    }
    public static double std(double[] data, double m) {
        double s=0; for(double d:data)s+=(d-m)*(d-m); return Math.sqrt(s/(data.length-1));
    }
    public static double[] downsample(double[] data, int factor) {
        int len = (int)Math.ceil((double)data.length/factor);
        double[] out = new double[len];
        for(int i=0;i<len;i++) out[i] = data[i*factor];
        return out;
    }
    public static double max(double[] data) {
        double m = -Double.MAX_VALUE; for(double d:data) if(d>m)m=d; return m;
    }
    public static double[] normalizeMax(double[] data) {
        double m = max(data); double[] out = new double[data.length];
        if(m==0) return out;
        for(int i=0;i<data.length;i++) out[i] = Math.abs(data[i])/m;
        return out;
    }
    public static List<Integer> findPeaks(double[] data) {
        List<Integer> p = new ArrayList<>();
        for(int i=1;i<data.length-1;i++) if(data[i]>data[i-1] && data[i]>data[i+1]) p.add(i);
        return p;
    }

    public static class Complex {
        public final double re;
        public final double im;
        public Complex(double re, double im) { this.re = re; this.im = im; }
        public Complex add(Complex b) { return new Complex(this.re + b.re, this.im + b.im); }
        public Complex sub(Complex b) { return new Complex(this.re - b.re, this.im - b.im); }
        public Complex mult(Complex b) { return new Complex(this.re * b.re - this.im * b.im, this.re * b.im + this.im * b.re); }
        public double abs() { return Math.hypot(re, im); }
        public double phase() { return Math.atan2(im, re); }
        // 新增共轭方法
        public Complex conjugate() { return new Complex(re, -im); }
        public Complex scale(double scalar) { return new Complex(re * scalar, im * scalar); }
    }

    public static Complex[] fft(double[] realInput, int n) {
        Complex[] x = new Complex[n];
        for (int i = 0; i < n; i++) {
            if (i < realInput.length) x[i] = new Complex(realInput[i], 0);
            else x[i] = new Complex(0, 0);
        }
        return fftRecursive(x);
    }

    // 实现 IFFT
    // IFFT(x) = conjugate(FFT(conjugate(x))) / N
    public static Complex[] ifft(Complex[] x) {
        int n = x.length;
        Complex[] conjInput = new Complex[n];
        for (int i = 0; i < n; i++) {
            conjInput[i] = x[i].conjugate();
        }

        Complex[] fftResult = fftRecursive(conjInput);

        Complex[] result = new Complex[n];
        for (int i = 0; i < n; i++) {
            result[i] = fftResult[i].conjugate().scale(1.0 / n);
        }
        return result;
    }

    private static Complex[] fftRecursive(Complex[] x) {
        int n = x.length;
        if (n == 1) return new Complex[]{x[0]};
        // 简单起见不检查 n power of 2，调用者需保证
        Complex[] even = new Complex[n / 2];
        Complex[] odd = new Complex[n / 2];
        for (int i = 0; i < n / 2; i++) {
            even[i] = x[2 * i];
            odd[i] = x[2 * i + 1];
        }
        Complex[] q = fftRecursive(even);
        Complex[] r = fftRecursive(odd);
        Complex[] y = new Complex[n];
        for (int k = 0; k < n / 2; k++) {
            double kth = -2 * k * Math.PI / n;
            Complex wk = new Complex(Math.cos(kth), Math.sin(kth));
            Complex times = wk.mult(r[k]);
            y[k] = q[k].add(times);
            y[k + n / 2] = q[k].sub(times);
        }
        return y;
    }
}
