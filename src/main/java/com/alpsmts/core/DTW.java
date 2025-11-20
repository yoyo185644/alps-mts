package com.alpsmts.core;

public class DTW {
    // Multivariate DTW with Sakoe-Chiba band.
    // Cost between frames is Euclidean norm across variables.
    public static double bandedDTW(double[][] Q, double[][] X, int a, int L, int band) {
        int Lq = Q.length;
        int V = Q[0].length;
        int T = L; // window length

        double INF = 1e300;
        double[][] dp = new double[Lq + 1][T + 1];
        for (int i = 0; i <= Lq; i++) {
            for (int j = 0; j <= T; j++) dp[i][j] = INF;
        }
        dp[0][0] = 0.0;

        for (int i = 1; i <= Lq; i++) {
            int jStart = Math.max(1, i - band);
            int jEnd   = Math.min(T, i + band);
            for (int j = jStart; j <= jEnd; j++) {
                double cost = frameDist(Q[i - 1], X[a + j - 1], V);
                double best = dp[i - 1][j];
                if (dp[i][j - 1] < best) best = dp[i][j - 1];
                if (dp[i - 1][j - 1] < best) best = dp[i - 1][j - 1];
                dp[i][j] = cost + best;
            }
        }
        return dp[Lq][T];
    }

    private static double frameDist(double[] q, double[] x, int V) {
        double s = 0.0;
        for (int v = 0; v < V; v++) {
            double d = q[v] - x[v];
            s += d * d;
        }
        return Math.sqrt(s);
    }
}
