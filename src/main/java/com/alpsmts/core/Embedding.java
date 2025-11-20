package com.alpsmts.core;

public class Embedding {
    public final RFF rff;
    public final PositionalEncoding pe;
    public final int Dx, Dp, D;

    public Embedding(RFF rff, PositionalEncoding pe) {
        this.rff = rff;
        this.pe = pe;
        this.Dx = rff.transform(new double[rff.V]).length;
        this.Dp = pe.encode(0).length;
        this.D = Dx * Dp;
    }

    public double[] phiXT(double[] x, int t) {
        double[] u = rff.transform(x); // Dx
        double[] p = pe.encode(t);     // Dp
        double[] out = new double[D];
        int idx = 0;
        for (int i = 0; i < Dx; i++) {
            for (int j = 0; j < Dp; j++) out[idx++] = u[i] * p[j];
        }
        return out;
    }

    public double[] windowMu(double[][] P, int a, int b) {
        int L = b - a + 1;
        double inv = 1.0 / Math.sqrt(L);
        double[] mu = new double[D];
        for (int d = 0; d < D; d++) mu[d] = (P[b][d] - P[a - 1][d]) * inv;
        return mu;
    }
}
