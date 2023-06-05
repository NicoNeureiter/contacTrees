/**
 *
 */
package contactrees;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
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
    protected boolean aheadOfMTree;

    @Override
    public void initAndValidate() {
        convIDs = new ArrayList<Integer>();
        convIDsStored = new ArrayList<>();
    }

    /**
     * Add a conversion over which this block moves in its history.
     * @param The new conversion.
     */
    public void addMove(Conversion conv) {
        startEditing(null);
        assert conv != null;
        assert convIDs != null;
        assert !convIDs.contains(conv.getID());

        convIDs.add(conv.getID());
    }


    /**
     * Remove a conversion move from this blocks history.
     * @param The conversion to be removed.
     */
    public void removeMove(Conversion conv) {
        startEditing(null);
        Integer cID = conv.getID();
        convIDs.remove(cID);
        assert !convIDs.contains(cID);
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

    @Override
    public int getDimension() {
        return 1;
    }

    /*
     * StateNode implementations
     */
    @Override
    public Block copy() {
        Block copy = new Block();

        copy.setID(getID());
        copy.convIDs = Util.deepCopyIntegers(convIDs);
        copy.convIDsStored = Util.deepCopyIntegers(convIDsStored);

        return copy;
    }

    @Override
    public void assignTo(StateNode other) {
        final Block block = (Block) other;
        block.setID(getID());
        block.convIDs = Util.deepCopyIntegers(convIDs);
        block.convIDsStored = Util.deepCopyIntegers(convIDsStored);
    }

    @Override
    public void assignFrom(StateNode other) {
        aheadOfMTree = true;
        final Block block = (Block) other;
        setID(block.getID());
        convIDs = Util.deepCopyIntegers(block.convIDs);
        convIDsStored = Util.deepCopyIntegers(block.convIDsStored);
    }

    @Override
    protected void store() {
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

        // Marginal trees are not stored/restored
        //  => restore leads to out-dated MTree
        aheadOfMTree = true;
    }

    @Override
    public void startEditing(Operator operator) {
        aheadOfMTree = true;
        super.startEditing(operator);
    }

    public void updatedMarginalTree() {
        aheadOfMTree = false;
    }

    public boolean isAheadOfMarginalTree() {
        return aheadOfMTree;
    }

    @Override
    public void assignFromFragile(StateNode other) {
        assignFrom(other);
    }

    @Override
    public double getArrayValue(int dim) {
        return convIDs.size();
    }

    @Override
    public void setEverythingDirty(boolean isDirty) {
        setSomethingIsDirty(isDirty);
    }

    @Override
    public void setSomethingIsDirty(final boolean isDirty) {
        aheadOfMTree = true;
        super.setSomethingIsDirty(isDirty);
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
    }

    @Override
    public int scale(double scale) {
        throw new UnsupportedOperationException("StateNode method scale() not applicable to Block.");
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
