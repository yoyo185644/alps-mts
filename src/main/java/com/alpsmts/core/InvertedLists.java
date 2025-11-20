package com.alpsmts.core;

import java.util.*;

public class InvertedLists {
    private final Map<String, List<Anchor>> map = new HashMap<>();

    public static String keyOf(long[] sig) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < sig.length; i++) {
            if (i > 0) sb.append(':');
            sb.append(Long.toHexString(sig[i]));
        }
        return sb.toString();
    }

    public void add(long[] sig, Anchor a) {
        String k = keyOf(sig);
        map.computeIfAbsent(k, kk -> new ArrayList<>()).add(a);
    }

    public List<Anchor> get(long[] sig) {
        return map.getOrDefault(keyOf(sig), Collections.emptyList());
    }
}
