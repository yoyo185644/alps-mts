package com.alpsmts.demo;

import com.alpsmts.core.*;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;
import java.util.concurrent.ThreadLocalRandom;

public class Main {
    public static void main(String[] args) {
//        int N = 5000;
//        int V = 5;
//        double[][] X = synthSeries(N, V);
//        System.out.println("Synthetic series generated: N=" + N + ", V=" + V);

        String filePath = "src/main/java/com/alpsmts/data/UEA/"; // 可以改成 "D:/test/data.csv" 或 "/home/user/data.csv"
        String fileName = "Crop_TRAIN.csv";

        double[][] X = readCSV(filePath, fileName);
        // Hyper-params
        int Dx = 64;
        int Dp = 8;
        double sigma = 2.0;
        int Lmin = 32;
        int Lmax = 512;
        double gamma = 2.0;
        int tau = 8;
        int lshBits = 64;
        int lshProbes = 8;
        int candidatesK = 64;
        int refineR = 8;
        double refineLenFrac = 0.1;

        boolean useZNorm = true;
        boolean useDTW = true;
        int rerankM = 50;
        double dtwBandRatio = 0.1;

        ALPSIndex index = new ALPSIndex(Dx, Dp, sigma,
                Lmin, Lmax, gamma, tau,
                lshBits, lshProbes, candidatesK,
                refineR, refineLenFrac,
                useZNorm, useDTW, rerankM, dtwBandRatio);

        System.out.println("Building index...");
        index.build(X);
        System.out.println("Index built.\n");

        // Query
        int qLen = 150;
        int qStart = 1000;
        double[][] Q = sliceWindow(X, qStart, qLen);

        System.out.println("Running query with L_q=" + qLen + "...");
        QueryResult res = index.query(Q);
        if (res != null) {
            System.out.printf(Locale.US,
                    "Best match: start=%d, length=%d, embedDist=%.6f, zNorm=%.6f, dtw=%s\n",
                    res.start, res.length, res.embedDistance, res.zNormDistance,
                    (res.dtwDistance == null ? "N/A" : String.format(Locale.US, "%.6f", res.dtwDistance)));
            int errStart = Math.abs(res.start - qStart);
            int errLen = Math.abs(res.length - qLen);
            System.out.printf("Sanity vs true: Δstart=%d, Δlen=%d\n", errStart, errLen);
        } else {
            System.out.println("No result.");
        }
    }

    static double[][] readCSV(String filePath, String fileName) {
        // 1. 指定文件路径（绝对路径 或 相对路径）


        List<double[]> dataList = new ArrayList<>();

        try (BufferedReader br = new BufferedReader(new FileReader(filePath + fileName))) {
            String line;
            while ((line = br.readLine()) != null) {
                // 用逗号分隔
                String[] values = line.split(",");

                // 转成 double[]
                double[] row = new double[values.length];
                for (int i = 0; i < values.length; i++) {
                    try {
                        row[i] = Double.parseDouble(values[i].trim());
                    } catch (NumberFormatException e) {
                        row[i] = Double.NaN; // 如果不是数字（如表头），用 NaN 占位
                    }
                }
                dataList.add(row);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        // 转成 double[][]
        double[][] dataArray = dataList.toArray(new double[0][]);

        // 打印结果（前几行）
        System.out.println("读取到的数据:");
        for (int i = 0; i < Math.min(5, dataArray.length); i++) {
            System.out.println(Arrays.toString(dataArray[i]));
        }
        // 4. 打印结果
        System.out.println("总行数: " + dataArray.length);
        System.out.println("最大列数: " + dataArray[0].length);
        return dataArray;
    }


    static double[][] synthSeries(int N, int V) {
        ThreadLocalRandom rnd = ThreadLocalRandom.current();
        double[][] X = new double[N][V];
        for (int t = 0; t < N; t++) {
            for (int v = 0; v < V; v++) {
                double trend = 0.01 * t;
                double season = Math.sin(2 * Math.PI * (t / 100.0 + v * 0.17));
                double noise = rnd.nextGaussian() * 0.2;
                X[t][v] = 0.5 * trend + 1.0 * season + noise;
            }
        }
        for (int t = 1000; t < 1150; t++) {
            for (int v = 0; v < V; v++) {
                X[t][v] += Math.sin(2 * Math.PI * (t / 12.0)) * (1.5 + 0.1 * v);
            }
        }
        return X;
    }

    static double[][] sliceWindow(double[][] X, int start, int L) {
        double[][] Q = new double[L][X[0].length];
        for (int i = 0; i < L; i++) {
            System.arraycopy(X[start + i], 0, Q[i], 0, X[0].length);
        }
        return Q;
    }
}
