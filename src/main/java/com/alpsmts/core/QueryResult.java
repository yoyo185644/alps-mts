package com.alpsmts.core;

public class QueryResult {
    public final int start;
    public final int length;
    public final double embedDistance;
    public final double zNormDistance;
    public final Double dtwDistance; // nullable if DTW disabled

    public QueryResult(int start, int length, double embedDistance, double zNormDistance, Double dtwDistance) {
        this.start = start;
        this.length = length;
        this.embedDistance = embedDistance;
        this.zNormDistance = zNormDistance;
        this.dtwDistance = dtwDistance;
    }
}
