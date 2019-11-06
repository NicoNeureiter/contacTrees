/**
 * 
 */
package contactrees;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import org.w3c.dom.Node;

import beast.core.Description;
import beast.core.Operator;
import beast.core.StateNode;
import contactrees.util.Util;


/**
 * @author Nico Neureiter <nico.neureiter@gmail.com>
 */
@Description("A block which is affected by a marginal tree and defined through a set of conversions wihtin a network.")
public class Block extends StateNode { 

    // The indices (in the ConversionList) of the conversions affecting this block. 
    List<Integer> convIDs, convIDsStored;

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
        convIDs.add(conv.getID());
    }
    
    
    /**
     * Remove a conversion move from this blocks history.
     * @param The conversion to be removed.
     */
    public void removeMove(Conversion conv) {
        startEditing(null);
        convIDs.remove((Integer) conv.getID());
    }

    /**
     * Evaluate whether this block was affected by the given conversion.
     * @param The conversion in question.
     * @return "true" iff the block moved over the conversion.
     */
    public boolean isAffected(Conversion conv) {
        // TODO resolve overshadowed conversions somewhere
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
        final Block block = (Block) other;
        setID(block.getID());
        convIDs = Util.deepCopyIntegers(block.convIDs);
        convIDsStored = Util.deepCopyIntegers(block.convIDsStored);
    }

    @Override
    protected void store() {
        convIDsStored = Util.deepCopyIntegers(convIDs);
    }

    @Override
    public void restore() {
//        System.out.println("Block restore");
        List<Integer> tmp = convIDs;
        convIDs = convIDsStored;
        convIDsStored = tmp;
    }
    
    @Override
    public void startEditing(Operator operator) {
        if (state != null)
            super.startEditing(operator);
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
    public void fromXML(Node node) {
        // TODO Auto-generated method stub
    }

    @Override
    public int scale(double scale) {
        throw new UnsupportedOperationException("StateNode method scale() not applicable to Block."); 
    }

}
