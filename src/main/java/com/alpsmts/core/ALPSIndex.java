package com.alpsmts.core;

import java.util.*;

public class ALPSIndex {
    // Hyperparams
    private final int Dx, Dp;
    private final double sigma;
    private final int Lmin, Lmax;
    private final double gamma;
    private final int tau;
    private final int lshBits, lshProbes;
    private final int candidatesK;
    private final int refineR;
    private final double refineLenFrac;

    // Re-ranking
    private final boolean useZNorm;
    private final boolean useDTW;
    private final int rerankM;
    private final double dtwBandRatio; // % of max(Lq, L)

    // Data dims
    private int N, V;

    // Components
    private RFF rff;
    private PositionalEncoding pe;
    private Embedding emb;

    // Prefix P[t] in R^D
    private double[][] P;

    private List<BucketIndex> buckets = new ArrayList<>();

    private double[][] X; // original series

    public ALPSIndex(int Dx, int Dp, double sigma,
                     int Lmin, int Lmax, double gamma, int tau,
                     int lshBits, int lshProbes, int candidatesK,
                     int refineR, double refineLenFrac,
                     boolean useZNorm, boolean useDTW, int rerankM, double dtwBandRatio) {
        this.Dx = Dx;
        this.Dp = Dp;
        this.sigma = sigma;
        this.Lmin = Lmin;
        this.Lmax = Lmax;
        this.gamma = gamma;
        this.tau = tau;
        this.lshBits = lshBits;
        this.lshProbes = lshProbes;
        this.candidatesK = candidatesK;
        this.refineR = refineR;
        this.refineLenFrac = refineLenFrac;
        this.useZNorm = useZNorm;
        this.useDTW = useDTW;
        this.rerankM = rerankM;
        this.dtwBandRatio = dtwBandRatio;
    }

    public void build(double[][] X) {
        this.X = X;
        this.N = X.length;
        this.V = X[0].length;

        // feature maps
        this.rff = new RFF(V, Dx, sigma, 12345L);
        this.pe  = new PositionalEncoding(Dp);
        this.emb = new Embedding(rff, pe);

        // prefix P
        this.P = new double[N + 1][emb.D];
        for (int t = 1; t <= N; t++) {
            double[] phi = emb.phiXT(X[t - 1], t);
            for (int d = 0; d < emb.D; d++) P[t][d] = P[t - 1][d] + phi[d];
        }

        // buckets
        List<Integer> Lset = geometricBuckets(Lmin, Lmax, gamma);
        // 固定随机种子，让每次构建可复现
        Random rnd = new Random(7);
        for (int L : Lset) {
            HyperplaneLSH lsh = new HyperplaneLSH(emb.D, lshBits, rnd.nextLong());
            BucketIndex bi = new BucketIndex(L, lsh);
            int step = Math.max(1, L / tau);
            for (int a = 1; a + L - 1 <= N; a += step) {
                double[] mu = emb.windowMu(P, a, a + L - 1);
                long[] sig = lsh.signature(mu);
                bi.inv.add(sig, new Anchor(a, L, mu));
            }
            buckets.add(bi);
        }
    }

    public QueryResult query(double[][] Q) {
        int Lq = Q.length;
        double[] muQ = muOfQuery(Q);

        Set<Integer> chosen = nearestBuckets(Lq);

        // Coarse candidates (keep K with smallest lb)
        PriorityQueue<Cand> topK = new PriorityQueue<>((c1, c2) -> Double.compare(c2.lb, c1.lb)); // max-heap by lb
        for (BucketIndex bi : buckets) {
            if (!chosen.contains(bi.L)) continue;
            long[] sigQ = bi.lsh.signature(muQ);
            List<long[]> probes = bi.lsh.multiProbe(sigQ, lshBits, 2, lshProbes);
            Set<String> seen = new HashSet<>();
            for (long[] probe : probes) {
                for (Anchor a : bi.inv.get(probe)) {
                    String key = a.start + "#" + a.length;
                    if (!seen.add(key)) continue;
                    double lb = Vec.norm2(Vec.sub(muQ, a.mu)); // coarse dist
                    topK.add(new Cand(a.start, a.length, lb));
                    if (topK.size() > candidatesK) topK.poll(); // drop worst
                }
            }
        }


        // Fallback: if no candidates from LSH, do a light brute-force scan on nearest buckets
        if (topK.isEmpty()) {
            for (BucketIndex bi : buckets) {
                if (!chosen.contains(bi.L)) continue;
                int step = Math.max(1, bi.L / Math.max(1, tau));
                for (int a = 1; a + bi.L - 1 <= N; a += step) {
                    double[] mu = emb.windowMu(P, a, a + bi.L - 1);
                    double d = Vec.norm2(Vec.sub(muQ, mu));
                    topK.add(new Cand(a, bi.L, d));
                    if (topK.size() > candidatesK) topK.poll(); // keep best K by smallest lb
                }
            }
        }

        // Local refinement: gather top M (by embed distance) during scanning
        PriorityQueue<Ref> bestRef = new PriorityQueue<>((r1, r2) -> Double.compare(r2.embedD, r1.embedD)); // max-heap
        while (!topK.isEmpty()) {
            Cand c = topK.poll();
            int L0 = c.L;
            int a0 = c.a;
            int dL = Math.max(1, (int)Math.round(L0 * refineLenFrac));
            for (int da = -refineR; da <= refineR; da++) {
                for (int dl = -dL; dl <= dL; dl++) {
                    int a = a0 + da;
                    int L = L0 + dl;
                    if (a < 1 || a + L - 1 > N || L <= 0) continue;
                    double[] mu = emb.windowMu(P, a, a + L - 1);
                    double d = Vec.norm2(Vec.sub(muQ, mu));
                    bestRef.add(new Ref(a, L, d));
                    if (bestRef.size() > Math.max(rerankM, 50)) bestRef.poll();
                }
            }
        }

        if (bestRef.isEmpty()) return null;

        // Re-ranking stage
        // First, pick the best by z-norm if enabled; otherwise by embed distance
        List<Ref> pool = new ArrayList<>();
        while (!bestRef.isEmpty()) pool.add(bestRef.poll());
        Collections.sort(pool, Comparator.comparingDouble(r -> r.embedD)); // ascending

        Ref best = pool.get(0);
        double bestEmbed = best.embedD;
        double bestZ = Double.POSITIVE_INFINITY;
        Double bestDTW = null;

        Ref bestByZ = best;
        boolean anyZ = false;
        if (useZNorm) {
            for (int i = 0; i < Math.min(pool.size(), rerankM); i++) {
                Ref r = pool.get(i);
                if (r.L != Lq) continue; // require equal length for z-norm
                double z = ZNorm.distance(Q, X, r.a - 1, r.L);
                if (z < bestZ) {
                    bestZ = z;
                    bestByZ = r;
                    anyZ = true;
                }
            }
            if (!anyZ) {
                bestZ = Double.NaN; // no equal-length candidate; skip z-norm
            }
        } else {
            bestZ = Double.NaN;
        }

        Ref bestByDTW = bestByZ;
        if (useDTW) {
            int Lcand = bestByZ.L;
            int band = (int)Math.round(Math.max(Lq, Lcand) * dtwBandRatio);
            double d = DTW.bandedDTW(Q, X, bestByZ.a - 1, Lcand, Math.max(1, band));
            bestDTW = d;

            // Try a few more top candidates for DTW (keep smallest)
            for (int i = 1; i < Math.min(pool.size(), 10); i++) {
                Ref r = pool.get(i);
                int bandI = (int)Math.round(Math.max(Lq, r.L) * dtwBandRatio);
                double di = DTW.bandedDTW(Q, X, r.a - 1, r.L, Math.max(1, bandI));
                if (di < bestDTW) {
                    bestDTW = di;
                    bestByDTW = r;
                }
            }
            // Report z-norm distance only if lengths equal AND we actually computed it
            Double zOut = (!Double.isNaN(bestZ) && bestByDTW.L == Lq && anyZ) ?
                          (bestByDTW == bestByZ ? bestZ : ZNorm.distance(Q, X, bestByDTW.a - 1, bestByDTW.L))
                          : Double.NaN;
            return new QueryResult(bestByDTW.a, bestByDTW.L, bestByDTW.embedD, zOut, bestDTW);
        } else {
            return new QueryResult(bestByZ.a, bestByZ.L, bestByZ.embedD, bestZ, null);
        }
    }

    /* -------- helpers -------- */

    private double[] muOfQuery(double[][] Q) {
        int Lq = Q.length;
        double[] sum = new double[emb.D];
        for (int i = 0; i < Lq; i++) {
            double[] phi = emb.phiXT(Q[i], i + 1);
            for (int d = 0; d < emb.D; d++) sum[d] += phi[d];
        }
        Vec.scaleInPlace(sum, 1.0 / Math.sqrt(Lq));
        return sum;
    }

    private List<Integer> geometricBuckets(int Lmin, int Lmax, double gamma) {
        List<Integer> Ls = new ArrayList<>();
        int L = Lmin;
        while (L <= Lmax) {
            Ls.add(L);
            L = (int)Math.max(L + 1, Math.round(L * gamma));
        }
        if (Ls.get(Ls.size() - 1) != Lmax) Ls.add(Lmax);
        return Ls;
    }

    private Set<Integer> nearestBuckets(int Lq) {
        int best = -1, bestDiff = Integer.MAX_VALUE;
        int second = -1, secondDiff = Integer.MAX_VALUE;
        for (BucketIndex bi : buckets) {
            int diff = Math.abs(bi.L - Lq);
            if (diff < bestDiff) {
                second = best; secondDiff = bestDiff;
                best = bi.L; bestDiff = diff;
            } else if (diff < secondDiff) {
                second = bi.L; secondDiff = diff;
            }
        }
        Set<Integer> s = new HashSet<>();
        if (best != -1) s.add(best);
        if (second != -1) s.add(second);
        return s;
    }

    static class Cand {
        final int a, L;
        final double lb;
        Cand(int a, int L, double lb) { this.a = a; this.L = L; this.lb = lb; }
    }

    static class Ref {
        final int a, L;
        final double embedD;
        Ref(int a, int L, double embedD) { this.a = a; this.L = L; this.embedD = embedD; }
    }
}
