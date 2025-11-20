package com.alpsmts.core;

import java.util.Arrays;

public class Vec {
    public static double[] addInPlace(double[] a, double[] b) {
        for (int i = 0; i < a.length; i++) a[i] += b[i];
        return a;
    }
    public static double[] sub(double[] a, double[] b) {
        double[] r = new double[a.length];
        for (int i = 0; i < a.length; i++) r[i] = a[i] - b[i];
        return r;
    }
    public static void scaleInPlace(double[] a, double c) {
        for (int i = 0; i < a.length; i++) a[i] *= c;
    }
    public static double norm2(double[] a) {
        double s = 0.0;
        for (double v : a) s += v * v;
        return Math.sqrt(s);
    }
    public static double dot(double[] a, double[] b) {
        double s = 0.0;
        for (int i = 0; i < a.length; i++) s += a[i] * b[i];
        return s;
    }
    public static double[] zeros(int d) { return new double[d]; }
    public static double[] copy(double[] a) { return Arrays.copyOf(a, a.length); }
}
