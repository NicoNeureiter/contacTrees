package contactrees.test.util;

import java.util.HashMap;

public class Counter<T> extends HashMap<T, Integer> {

    public void increment(T key) {
        this.put(key, this.getOrDefault(key, 0) + 1);
    }

    public int getCount(T key) {
        return this.getOrDefault(key, 0);
    }

    public void clear(T key) {
        if (this.containsKey(key)) {
            this.remove(key);
        }
    }

    public void resetAll() {
        this.clear();
    }

    @Override
    public String toString() {
        return "Counter{" + super.toString() + '}';
    }

    public HashMap<T, Double> getFrequencies() {
        HashMap<T, Double> frequencies = new HashMap<>();

        // Get the total count
        int totalCount = this.values().stream().mapToInt(Integer::intValue).sum();

        // Avoid division by zero
        if (totalCount == 0) {
            return frequencies;
        }

        // Put the normalized counts in the frequencies map
        for (T key : this.keySet()) {
            frequencies.put(key, this.get(key) / (double) totalCount);
        }

        return frequencies;
    }

}
