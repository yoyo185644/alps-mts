package com.alpsmts.core;

public class BucketIndex {
    public final int L;
    public final HyperplaneLSH lsh;
    public final InvertedLists inv;

    public BucketIndex(int L, HyperplaneLSH lsh) {
        this.L = L;
        this.lsh = lsh;
        this.inv = new InvertedLists();
    }
}
