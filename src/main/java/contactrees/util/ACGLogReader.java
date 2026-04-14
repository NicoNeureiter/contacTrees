package contactrees.util;


import contactrees.ACGWithBlocks;

/**
 * @author Nico Neureiter
 */
public interface ACGLogReader extends Iterable<ACGWithBlocks> {

    int getACGCount();
    int getCorrectedACGCount();
}
