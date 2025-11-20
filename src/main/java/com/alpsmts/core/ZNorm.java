package com.alpsmts.core;

public class ZNorm {
    // Multivariate z-norm Euclidean distance between Q and X[a..a+L-1]
    // Q: Lq x V, Window: L x V (must be equal length here; if not, caller should align/truncate/pad)
    public static double distance(double[][] Q, double[][] X, int a, int L) {
        int Lq = Q.length;
        if (L != Lq) {
            throw new IllegalArgumentException("ZNorm distance requires equal lengths. got L=" + L + ", Lq=" + Lq);
        }
        int V = Q[0].length;
        final double EPS = 1e-8;

        // Means and stds
        double[] meanQ = new double[V];
        double[] meanW = new double[V];
        for (int i = 0; i < Lq; i++) {
            for (int v = 0; v < V; v++) {
                meanQ[v] += Q[i][v];
                meanW[v] += X[a + i][v];
            }
        }
        for (int v = 0; v < V; v++) {
            meanQ[v] /= Lq;
            meanW[v] /= Lq;
        }
        double[] stdQ = new double[V];
        double[] stdW = new double[V];
        for (int i = 0; i < Lq; i++) {
            for (int v = 0; v < V; v++) {
                double dq = Q[i][v] - meanQ[v];
                double dw = X[a + i][v] - meanW[v];
                stdQ[v] += dq * dq;
                stdW[v] += dw * dw;
            }
        }
        for (int v = 0; v < V; v++) {
            stdQ[v] = Math.sqrt(stdQ[v] / Math.max(1, Lq - 1)) + EPS;
            stdW[v] = Math.sqrt(stdW[v] / Math.max(1, Lq - 1)) + EPS;
        }

        // Sum squared diffs after z-norm
        double sum = 0.0;
        for (int i = 0; i < Lq; i++) {
            for (int v = 0; v < V; v++) {
                double qn = (Q[i][v] - meanQ[v]) / stdQ[v];
                double wn = (X[a + i][v] - meanW[v]) / stdW[v];
                double diff = qn - wn;
                sum += diff * diff;
            }
        }
        return Math.sqrt(sum);
    }
}
