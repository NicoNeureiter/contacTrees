/**
 *
 */
package contactrees;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;

import org.w3c.dom.Node;

import beast.base.core.Description;
import beast.base.inference.Operator;
import beast.base.inference.StateNode;
import contactrees.util.Util;


/**
 * @author Nico Neureiter
 */
@Description("A block which is affected by a marginal tree and defined through a set of conversions wihtin a network.")
public class Block extends StateNode {

    // The indices (in the ConversionList) of the conversions affecting this block.
    List<Integer> convIDs, convIDsStored;

    // The same IDs as a set, for O(1) membership queries.
    protected HashSet<Integer> convIDSet = new HashSet<>();

    // Version of the current state (see ConversionGraph#getVersion()).
    protected long version = ConversionGraph.nextVersion(), storedVersion = version;

    @Override
    public void initAndValidate() {
        convIDs = new ArrayList<Integer>();
        convIDsStored = new ArrayList<>();
        convIDSet.clear();
        version = ConversionGraph.nextVersion();
    }

    private void rebuildIDSet() {
        convIDSet.clear();
        convIDSet.addAll(convIDs);
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
        assert convIDs != null;
        assert !convIDSet.contains(conv.getID());

        convIDs.add(conv.getID());
        convIDSet.add(conv.getID());
    }


    /**
     * Remove a conversion move from this blocks history.
     * @param The conversion to be removed.
     */
    public void removeMove(Conversion conv) {
        startEditing(null);
        Integer cID = conv.getID();
        convIDs.remove(cID);
        convIDSet.remove(cID);
        assert !convIDs.contains(cID);
    }

    /**
     * Remove all conversion moves from this blocks history.
     */
    public void removeAllMoves() {
        startEditing(null);
        convIDs.clear();
        convIDSet.clear();
    }

    /**
     * Evaluate whether this block was affected by the given conversion.
     * @param The conversion in question.
     * @return "true" iff the block moved over the conversion.
     */
    public boolean isAffected(Conversion conv) {
        return convIDSet.contains(conv.getID());
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

    public List<Integer> getConversionIDs(){
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
        copy.convIDs = Util.deepCopyIntegers(convIDs);
        copy.convIDsStored = Util.deepCopyIntegers(convIDsStored);
        copy.rebuildIDSet();

        return copy;
    }

    @Override
    public void assignTo(StateNode other) {
        final Block block = (Block) other;
        block.setID(getID());
        block.convIDs = Util.deepCopyIntegers(convIDs);
        block.convIDsStored = Util.deepCopyIntegers(convIDsStored);
        block.rebuildIDSet();
        block.version = ConversionGraph.nextVersion();
    }

    @Override
    public void assignFrom(StateNode other) {
        final Block block = (Block) other;
        setID(block.getID());
        convIDs = Util.deepCopyIntegers(block.convIDs);
        convIDsStored = Util.deepCopyIntegers(block.convIDsStored);
        rebuildIDSet();
        version = ConversionGraph.nextVersion();
    }

    @Override
    protected void store() {
        storedVersion = version;
        convIDsStored.clear();
        for (int i=0; i<convIDs.size(); i++) {
            convIDsStored.add(convIDs.get(i));
        }
    }

    @Override
    public void restore() {
        List<Integer> tmp = convIDs;
        convIDs = convIDsStored;
        convIDsStored = tmp;
        rebuildIDSet();
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
        // Initialize arrays
        convIDs = new ArrayList<Integer>();
        convIDsStored = new ArrayList<>();

        // Get the string representation of the `convIDs` list
        String str = node.getTextContent();

        // Remove the square brackets at the start and end.
        assert str.startsWith("[");
        assert str.endsWith("]");
        str = str.substring(1, str.length() - 1);

        // Parse the string list and put convIDs into the `convIDs` array.
        if (str.length() > 0) {
            for (String convID : str.split("\\s*,\\s*")) {
                convIDs.add(Integer.parseInt(convID));
            }
        }
        rebuildIDSet();
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
