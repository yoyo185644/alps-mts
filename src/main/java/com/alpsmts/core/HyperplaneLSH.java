package com.alpsmts.core;

import java.util.*;

public class HyperplaneLSH {
    public final int dim;
    public final int bits;
    private final double[][] planes; // bits x dim

    public HyperplaneLSH(int dim, int bits, long seed) {
        this.dim = dim;
        this.bits = bits;
        this.planes = new double[bits][dim];
        Random rnd = new Random(seed);
        for (int b = 0; b < bits; b++) {
            double norm = 0.0;
            for (int d = 0; d < dim; d++) {
                double val = rnd.nextGaussian();
                planes[b][d] = val;
                norm += val * val;
            }
            norm = Math.sqrt(norm) + 1e-12;
            for (int d = 0; d < dim; d++) planes[b][d] /= norm;
        }
    }

    public long[] signature(double[] v) {
        int words = (bits + 63) / 64;
        long[] sig = new long[words];
        for (int b = 0; b < bits; b++) {
            double s = 0.0;
            double[] w = planes[b];
            for (int d = 0; d < dim; d++) s += w[d] * v[d];
            int word = b >>> 6;
            int off = b & 63;
            if (s >= 0) sig[word] |= (1L << off);
        }
        return sig;
    }

    public static int hamming(long[] a, long[] b) {
        int d = 0;
        for (int i = 0; i < a.length; i++) d += Long.bitCount(a[i] ^ b[i]);
        return d;
    }

    public List<long[]> multiProbe(long[] base, int bits, int R, int maxProbes) {
        List<long[]> probes = new ArrayList<>();
        probes.add(base);
        if (R <= 0) return probes;

        List<Integer> bitPositions = new ArrayList<>();
        for (int b = 0; b < bits; b++) bitPositions.add(b);
        Collections.shuffle(bitPositions, new Random(42));

        int added = 1;
        // 1-bit flips
        for (int i = 0; i < bitPositions.size() && added < maxProbes; i++) {
            long[] c1 = base.clone();
            flipBit(c1, bitPositions.get(i));
            probes.add(c1);
            added++;
        }
        if (R >= 2) {
            outer:
            for (int i = 0; i < bitPositions.size() && added < maxProbes; i++) {
                for (int j = i + 1; j < bitPositions.size() && added < maxProbes; j++) {
                    long[] c2 = base.clone();
                    flipBit(c2, bitPositions.get(i));
                    flipBit(c2, bitPositions.get(j));
                    probes.add(c2);
                    added++;
                    if (added >= maxProbes) break outer;
                }
            }
        }
        return probes;
    }

    private static void flipBit(long[] sig, int bitIndex) {
        int word = bitIndex >>> 6;
        int off = bitIndex & 63;
        sig[word] ^= (1L << off);
    }
}
