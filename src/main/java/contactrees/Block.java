/**
 *
 */
package contactrees;

import java.io.PrintStream;
import java.util.Arrays;
import java.util.LinkedHashSet;
import java.util.Set;

import org.w3c.dom.Node;

import beast.base.core.Description;
import beast.base.inference.Operator;
import beast.base.inference.StateNode;


/**
 * @author Nico Neureiter
 */
@Description("A block which is affected by a marginal tree and defined through a set of conversions wihtin a network.")
public class Block extends StateNode {

    // The indices (in the ConversionList) of the conversions affecting this block, in insertion order.
    protected LinkedHashSet<Integer> convIDs = new LinkedHashSet<>(), convIDsStored = new LinkedHashSet<>();

    // Version of the current state (see ConversionGraph#getVersion()).
    protected long version = ConversionGraph.nextVersion(), storedVersion = version;

    @Override
    public void initAndValidate() {
        convIDs = new LinkedHashSet<>();
        convIDsStored = new LinkedHashSet<>();
        version = ConversionGraph.nextVersion();
    }

    /**
     * @return the version of the current state of this block.
     */
    public long getVersion() {
        return version;
    }

    /**
     * Add a conversion over which this block moves in its history.
     * @param The new conversion.
     */
    public void addMove(Conversion conv) {
        startEditing(null);
        assert conv != null;
        boolean added = convIDs.add(conv.getID());
        assert added;
    }


    /**
     * Remove a conversion move from this blocks history.
     * @param The conversion to be removed.
     */
    public void removeMove(Conversion conv) {
        startEditing(null);
        convIDs.remove(conv.getID());
    }

    /**
     * Remove all conversion moves from this blocks history.
     */
    public void removeAllMoves() {
        startEditing(null);
        convIDs.clear();
    }

    /**
     * Evaluate whether this block was affected by the given conversion.
     * @param The conversion in question.
     * @return "true" iff the block moved over the conversion.
     */
    public boolean isAffected(Conversion conv) {
        return convIDs.contains(conv.getID());
    }

    public int countMoves() {
        return convIDs.size();
    }

    /**
     * Number of conversions over which this block moves.
     * @return Number of conversions.
     */
    public int size() {
        return convIDs.size();
    }

    public Set<Integer> getConversionIDs(){
        return convIDs;
    }


    /*
     * Loggable implementations
     */

    /**
     * init() and close() can stay empty.
     */
    @Override
    public void init(PrintStream out) {}
    @Override
    public void close(PrintStream out) {}

    /*
     * StateNode implementations
     */
    @Override
    public Block copy() {
        Block copy = new Block();

        copy.setID(getID());
        copy.convIDs = new LinkedHashSet<>(convIDs);
        copy.convIDsStored = new LinkedHashSet<>(convIDsStored);

        return copy;
    }

    @Override
    public void assignTo(StateNode other) {
        final Block block = (Block) other;
        block.setID(getID());
        block.convIDs = new LinkedHashSet<>(convIDs);
        block.convIDsStored = new LinkedHashSet<>(convIDsStored);
        block.version = ConversionGraph.nextVersion();
    }

    @Override
    public void assignFrom(StateNode other) {
        final Block block = (Block) other;
        setID(block.getID());
        convIDs = new LinkedHashSet<>(block.convIDs);
        convIDsStored = new LinkedHashSet<>(block.convIDsStored);
        version = ConversionGraph.nextVersion();
    }

    @Override
    protected void store() {
        storedVersion = version;
        convIDsStored.clear();
        convIDsStored.addAll(convIDs);
    }

    @Override
    public void restore() {
        LinkedHashSet<Integer> tmp = convIDs;
        convIDs = convIDsStored;
        convIDsStored = tmp;
        version = storedVersion;
    }

    @Override
    public void startEditing(Operator operator) {
        super.startEditing(operator);
        version = ConversionGraph.nextVersion();
    }

    @Override
    public void assignFromFragile(StateNode other) {
        assignFrom(other);
    }

    // public double getArrayValue(int dim) {
    //     return convIDs.size();
    // }

    @Override
    public void setEverythingDirty(boolean isDirty) {
        setSomethingIsDirty(isDirty);
    }

    @Override
    public void fromXML(Node node) {
        // Initialize the ID sets
        convIDs = new LinkedHashSet<>();
        convIDsStored = new LinkedHashSet<>();

        // Get the string representation of the `convIDs` list
        String str = node.getTextContent();

        // Remove the square brackets at the start and end.
        assert str.startsWith("[");
        assert str.endsWith("]");
        str = str.substring(1, str.length() - 1);

        // Parse the string list and put the IDs into the `convIDs` set.
        if (str.length() > 0) {
            for (String convID : str.split("\\s*,\\s*")) {
                boolean added = convIDs.add(Integer.parseInt(convID));
                assert added;
            }
        }
        version = ConversionGraph.nextVersion();
    }

    public Block() {
        super();
    }

    @Override
    public String toString() {
        return convIDs.toString();
    }

    /**
     * Create a new block with specified BEASTObject-ID.
     * Used for parsing log files.
     * @param ID
     */
    public Block(String ID) {
        this();
        setID(ID);
    }

    @Override
    public int getChecksum() {
        Integer[] hashValues = new Integer[2];
        hashValues[0] = countMoves();
        hashValues[1] = 0;
        for (int cID : getConversionIDs())
            hashValues[1] = 31 * hashValues[1] + cID;

        return Arrays.deepHashCode(hashValues);
    }

	@Override
	public void log(long sample, PrintStream out) {
		out.print(this.toString() + "\t");
	}

}
