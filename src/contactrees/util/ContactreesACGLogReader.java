package contactrees.util;

import contactrees.ACGWithBlocks;
import contactrees.Block;
import contactrees.BlockSet;
import contactrees.ConversionGraph;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.antlr.v4.runtime.misc.Pair;

/**
 * Class representing ACG log files.  Includes methods for
 * querying the number of ACGs defined, included and excluded
 * by the given burn-in percentage, as well as implementing an
 * iterator over all ACGs included after burn-in.  The iterator
 * automatically displays a progress bar on stdout.
 *
 * @author Nico Neureiter
 */
public class ContactreesACGLogReader implements ACGLogReader {
    File logFile;
    BufferedReader reader;

    List<String> preamble, postamble;
    String nextLine;

    int nACGs, burnin;
    
    ArrayList<Block> blocks;

    /**
     * Construct and initialize the reader.  The Preamble is
     * read and the list of loci constructed immediately.
     *
     * @param logFile ACG log file.
     * @throws IOException
     */
    public ContactreesACGLogReader(File logFile, double burninPercentage) throws IOException {
        this.logFile = logFile;

        reader = new BufferedReader(new FileReader(logFile));

        preamble = new ArrayList<>();
        skipPreamble();

        nACGs = 0;
        while (true) {
            if (getNextTreeString() == null)
                break;

            nACGs += 1;
        }
        burnin = (int) Math.round(nACGs*burninPercentage/100);

        postamble = new ArrayList<>();
        readPostamble();

        blocks = new ArrayList<>();
        extractBlocks();
    }


    /**
     * Internal method for skimming the preamble at the start
     * of the log, before we get to the tree section.
     *
     * @throws IOException
     */
    private void skipPreamble() throws IOException {
        boolean recordPreamble = preamble.isEmpty();

        while(true) {
            nextLine = reader.readLine();

            if (nextLine == null)
                throw new IOException("Reached end of file while searching for first tree.");

            nextLine = nextLine.trim();

            if (nextLine.toLowerCase().startsWith("tree"))
                break;

            if (recordPreamble)
                preamble.add(nextLine);
        }
    }

    /**
     * Internal method for extracting postamble following trees.
     *
     * @throws IOException
     */
    private void readPostamble() throws IOException {
        while (true) {
            if (nextLine == null)
                break;

            postamble.add(nextLine);

            nextLine = reader.readLine();
        }
    }

    /**
     * @return Everything read from the log file up until the first tree line.
     */
    public String getPreamble() {
        StringBuilder sb = new StringBuilder();
        for (String line : preamble)
            sb.append(line).append("\n");

        return sb.toString();
    }

    /**
     * @return Everything read from the log file following the last tree line.
     */
    public String getPostamble() {
        StringBuilder sb = new StringBuilder();
        for (String line : postamble)
            sb.append(line).append("\n");

        return sb.toString();
    }

    /**
     * Retrieve list of blocks from preamble or postamble.
     */
    private void extractBlocks() {        
        List<String> prepost = new ArrayList<>();
        prepost.addAll(preamble);
        prepost.addAll(postamble);

        final String PREFIX = "blockSet "; 
        for (String line : prepost) {
            line = line.trim();
            
            if (line.startsWith(PREFIX) && line.endsWith(";")) {
                line = line.substring(PREFIX.length(), line.length()-1);
                
                for (String blockID : line.split(" ")) {
                    blocks.add(new Block(blockID));
                }
            }
        }
    }

    /**
     * Rewind to the beginning of the file.
     *
     * @throws IOException
     */
    private void reset() throws IOException {
        reader.close();
        reader = new BufferedReader(new FileReader(logFile));
        skipPreamble();
    }

    /**
     * @return the next available tree string or null if none exists
     * @throws IOException
     */
    private String getNextTreeString() throws IOException {
        StringBuilder sb = new StringBuilder();

        while (true) {
            if (nextLine == null || nextLine.trim().toLowerCase().equals("end;"))
                return null;

            sb.append(nextLine.trim());
            if (nextLine.trim().endsWith(";"))
                break;

            nextLine = reader.readLine();
        }
        nextLine = reader.readLine();

        String treeString = sb.toString();

        return treeString.substring(treeString.indexOf("("));
    }

    /**
     * Skip burn-in portion of log.
     *
     * @throws IOException
     */
    private void skipBurnin() throws IOException {
        for (int i=0; i<burnin; i++)
            getNextTreeString();
    }

    /**
     * @return total number of ACGs defined by file.
     */
    public int getACGCount() {
        return nACGs;
    }

    /**
     * @return number of ACGs excluded as burn-in
     */
    public int getBurnin() {
        return burnin;
    }

    /**
     * @return number of ACGs excluding burn-in
     */
    public int getCorrectedACGCount() {
        return nACGs - burnin;
    }

    /**
     * Retrieve an iterator for iterating over the ACGs and corresponding 
     * blockSets represented by this log file. Important points
     *
     * 1. The iterator only iterates over as many (non-burnin) ACGs as exist
     * in the file when the ACGLogFileReader is constructed.  This is to avoid
     * problems associated with summarising ongoing analyses.
     *
     * 2. The iterator reuses a single ConversionGraph object during
     * the iteration.  This means that if you want to collect these
     * graphs as the iteration progresses you'll need to use
     * ConversionGraph::copy.
     *
     * @return ConversionGraph iterator
     */
    @Override
    public Iterator<ACGWithBlocks> iterator() {
        try {
            reset();
            skipBurnin();
        } catch (IOException e) {
            throw new IllegalStateException(e.getMessage());
        }

        ACGWithBlocks acgWithBlocks;
//        ACGWithBlocks acgWithBlocks = new ACGWithBlocks();
        try {
            acgWithBlocks = ACGWithBlocks.newFromNewick(blocks);
        } catch (Exception e) {
            throw new IllegalStateException(e.getMessage());
        }

        return new Iterator<ACGWithBlocks>() {

            boolean lineConsumed = true;
            String nextLine = null;

            int current = 0;

            private String getNextLineNoConsume() {
                if (lineConsumed) {
                    try {
                        nextLine = getNextTreeString();
                        lineConsumed = false;
                    } catch (IOException e) {
                        throw new IllegalStateException(e.getMessage());
                    }
                }

                return nextLine;
            }

            private void printProgressBar() {

                if (current==0) {
                    System.out.println("0%             25%            50%            75%           100%");
                    System.out.println("|--------------|--------------|--------------|--------------|");
                }

                if (current < getCorrectedACGCount()-1) {
                    if (current % (int) Math.ceil(getCorrectedACGCount() / 61.0) == 0) {
                        System.out.print("\r");
                        for (int i = 0; i < Math.round(61.0 * current / getCorrectedACGCount()); i++)
                            System.out.print("*");
                        System.out.flush();
                    }
                } else {
                    System.out.print("\r");
                    for (int i=0; i<61; i++)
                        System.out.print("*");
                    System.out.println();
                }
            }

            @Override
            public boolean hasNext() {
                return current<getCorrectedACGCount() && getNextLineNoConsume() != null;
            }

            @Override
            public ACGWithBlocks next() {
                String result = getNextLineNoConsume();
                lineConsumed = true;
                
                for (Block block : acgWithBlocks.blockSet)
                    block.removeAllMoves();
                acgWithBlocks.removeAllConversions();
                
                acgWithBlocks.fromExtendedNewick(result);

                printProgressBar();
                current += 1;

                return acgWithBlocks;
            }
        };
    }
}

