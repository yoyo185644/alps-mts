package com.alpsmts.core;

public class PositionalEncoding {
    private final int Dp;
    private final double[] freqs; // Dp/2
    private final double scale;

    public PositionalEncoding(int Dp) {
        if (Dp % 2 != 0) throw new IllegalArgumentException("Dp must be even.");
        this.Dp = Dp;
        this.freqs = new double[Dp / 2];
        double f0 = 2 * Math.PI / 8.0;
        double growth = 2.0;
        for (int k = 0; k < freqs.length; k++) freqs[k] = f0 * Math.pow(growth, k);
        this.scale = Math.sqrt(2.0 / Dp);
    }

    public double[] encode(int t) {
        double[] out = new double[Dp];
        int idx = 0;
        for (int k = 0; k < freqs.length; k++) {
            double z = freqs[k] * t;
            out[idx++] = scale * Math.cos(z);
            out[idx++] = scale * Math.sin(z);
        }
        return out;
    }
}
