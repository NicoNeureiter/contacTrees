package contactrees.test.util;

import java.util.HashMap;
import java.util.concurrent.Callable;

public class DefaultHashMap<K, V> extends HashMap<K, V> {

    Callable<V> defaultFactory;

    public DefaultHashMap(Callable<V> defaultFactory) {
        this.defaultFactory = defaultFactory;
    }

    public DefaultHashMap(V defaultValue) {
        this.defaultFactory = (() -> defaultValue);
    }

    @Override
    public V get(Object key) {
        V returnValue = super.get(key);
        if (returnValue == null) {
            try {
                returnValue = defaultFactory.call();
            } catch (Exception e) {
                throw new RuntimeException(e);
            }
            this.put((K) key, returnValue);
        }
        return returnValue;
    }
}
