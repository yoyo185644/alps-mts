# ALPS-MTS (Java, Maven)

Any-Length Prefix Sketch for Multivariate Time Series:
- Multivariate, any-length subsequence search
- Compact index (geometric length buckets + LSH)
- Fast query with local refinement
- Re-ranking with z-normalized Euclidean and banded DTW

## Build

```bash
cd alps-mts
mvn -q -DskipTests package
# fat jar:
java -jar target/alps-mts-1.0-SNAPSHOT-jar-with-dependencies.jar
```

## Code layout

- `com.alpsmts.core`: core algorithms
- `com.alpsmts.demo.Main`: runnable demo (synthetic data)

You can adapt `Main` to load your own datasets.
