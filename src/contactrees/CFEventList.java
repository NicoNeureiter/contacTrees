package contactrees;

import beast.evolution.tree.Node;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import javax.management.InvalidAttributeValueException;

/**
 * Maintains an ordered list of events which make up the clonal frame.
 *
 * @author Nico Neureiter
 */
public class CFEventList {

    /**
     * Class of events types on clonal frame.
     */
    public enum EventType {COALESCENCE, SAMPLE }

    /**
     * Class of events on clonal frame.
     */
    public static class Event {
        EventType type;
        double t;
        int lineages;
        Node node;

        /**
         * Construct event object corresponding to chosen CF node.
         *
         * @param node chosen CF node
         */
        public Event(Node node) {
            this.node = node;

            if (node.isLeaf())
                type = EventType.SAMPLE;
            else
                type = EventType.COALESCENCE;
            
            t = node.getHeight();
        }

        /**
         * Construct an event object with only the time defined.
         * Useful for conducting binary searches over the event list.
         *
         * @param t age of the event object
         */
        public Event(double t) {
            this.t = t;
        }

        public double getHeight() {
            return t;
        }
        
        public EventType getType() {
            return type;
        }

        public Node getNode() {
            return node;
        }
        
        /**
         * @return number of lineages _above_ this event.
         */
        public int getLineageCount() {
            return lineages;
        }
        
        @Override
        public String toString() {
            return "t: " + t + ", k: " + lineages + ", type: " + type;
        }
    } 

    /**
     * Ancestral conversion graph this list belongs to.
     */
    private final ConversionGraph acg;

    /**
     * List of events on clonal frame.
     */
    private final List<Event> events;
    private boolean dirty;

    public CFEventList(ConversionGraph acg) {
        this.acg = acg;
        
        events = new ArrayList<>();
        dirty = true;
    }

    /**
     * Obtain ordered list of events that make up the clonal frame.  Used
     * for ARG probability density calculations and for various state proposal
     * operators.
     * 
     * @return List of events.
     */
    public List<Event> getCFEvents() {
        updateEvents();
        return events;
    }

    /**
     * Obtain the number of events in the list.
     *
     * @return Number of events
     */
    public int countEvents() {
        return events.size();
    }

    /**
     * Mark the event list as dirty.
     */
    public void makeDirty() {
        dirty = true;
    }

    /**
     * Assemble sorted list of events on clonal frame and a map from nodes
     * to these events.
     */
    public void updateEvents() {
        // TODO Fix dirty logic
//        if (!dirty)
//            return;
        
        events.clear();
        
        // Create event list
        for (Node node : acg.getNodesAsArray()) {
            Event event = new Event(node);
            events.add(event);
        }
        
        // Sort events in increasing order of their heights
        Collections.sort(events, (Event o1, Event o2) -> {
            if (o1.t<o2.t)
                return -1;
            
            if (o2.t<o1.t)
                return 1;

            if (o1.getNode().getParent() == o2.getNode())
                return -1;

            if (o1.getNode() == o2.getNode().getParent())
                return 1;

            return 0;
        });
        
        // Compute lineage counts:
        int k=0;
        for (Event event : events) {
            if (event.type == EventType.SAMPLE)
                k += 1;
            else
                k -= 1;
            
            event.lineages = k;
        }


        dirty = false;
    }
    
    /**
     * Get the last event below the given height.
     * @param height
     * @return
     * @throws InvalidAttributeValueException 
     */
    public Event getEventAtHeight(double height) {
        updateEvents();
        int startIdx = 0;
        while ((startIdx < events.size()-1) && (events.get(startIdx+1).getHeight()<height))
            startIdx += 1;
        
        return events.get(startIdx);
    }
    
    public double getIntervalVolume(int i) {
        int k = events.get(i).lineages;
        double dt = events.get(i+1).t - events.get(i).t;
        return dt * k * (k-1);
    }

    public double[] getIntervalVolumes() {
        int nEvents = events.size();
        double[] volumes = new double[nEvents-1];
        for (int i=0; i<nEvents-1; i++) {
            volumes[i] = getIntervalVolume(i);
        }
        
        return volumes;
    }
    
}
