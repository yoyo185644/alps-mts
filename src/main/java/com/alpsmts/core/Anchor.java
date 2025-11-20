package com.alpsmts.core;

public class Anchor {
    public final int start;     // 1-based
    public final int length;    // L
    public final double[] mu;   // embedding

    public Anchor(int start, int length, double[] mu) {
        this.start = start;
        this.length = length;
        this.mu = mu;
    }
}
