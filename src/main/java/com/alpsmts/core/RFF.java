package com.alpsmts.core;

import java.util.Random;

public class RFF {
    public final int V;
    public final int Dx;
    public final double sigma;
    private final double[][] omega; // V x Dx
    private final double[] b;       // Dx
    private final double scale;

    public RFF(int V, int Dx, double sigma, long seed) {
        this.V = V;
        this.Dx = Dx;
        this.sigma = sigma;
        this.omega = new double[V][Dx];
        this.b = new double[Dx];
        this.scale = Math.sqrt(2.0 / Dx);
        Random rnd = new Random(seed);
        double std = 1.0 / sigma;
        for (int j = 0; j < Dx; j++) {
            // 给每一个纬度都生成符合高斯分布的随机矩阵 Dx*V
            for (int i = 0; i < V; i++) omega[i][j] = std * rnd.nextGaussian();
            // 给每一个纬度都一个bias 1*Dx
            b[j] = rnd.nextDouble() * 2.0 * Math.PI;
        }
    }

    public double[] transform(double[] x) {
        double[] out = new double[Dx];
        for (int j = 0; j < Dx; j++) {
            double z = b[j];
            for (int i = 0; i < V; i++) z += omega[i][j] * x[i];
            out[j] = scale * Math.cos(z);
        }
        return out;
    }
}
