/**
 *
 */
package contactrees.util;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;

import beast.evolution.tree.Node;
import beast.util.Randomizer;

/**
 *
 *
 * @author Nico Neureiter
 */
public class Util {

	/**
	 *
	 */
	public Util() {}

	static public ArrayList<Integer> deepCopyIntegers(Iterable<Integer> original){
	    ArrayList<Integer> copy = new ArrayList<>();
	    for (Integer x : original) copy.add(x);
	    return copy;
	}

	static public <T> T sampleFrom(Collection<T> population) {
		int z = Randomizer.nextInt(population.size());
		for (T candidate : population) {
			if (z == 0)
				return candidate;
			z--;
		}
		assert population.size() == 0;
		throw new RuntimeException("Can not sample from empty set.");
	}

	static public void sortByHeight(List<Node> nodes, boolean reverse) {
		if (reverse) {
	        Collections.sort(nodes, (Node n1, Node n2) -> {
	            if (n1.getHeight()<n2.getHeight())
	                return -1;
	            if (n1.getHeight()>n2.getHeight())
	                return 1;
	            return 0;
	        });
		} else {
	        Collections.sort(nodes, (Node n1, Node n2) -> {
	            if (n1.getHeight()<n2.getHeight())
	                return 1;
	            if (n1.getHeight()>n2.getHeight())
	                return -1;
	            return 0;
	        });
		}

	}

	static public double logFactorial(int n) {
		double res = 0;
		for (int i=1; i<=n; i++) {
			res += Math.log(i);
		}
		return res;
	}

	static public double sum(double[] values) {
		double total = 0;
		for (double x : values) total += x;
		return total;
	}

	static public double sum(Iterable<Double> values) {
		double total = 0;
		for (double x : values) total += x;
		return total;
	}

	static public int sampleCategorical(double[] weights) {
		double weightSum = sum(weights);
		double u = weightSum * Randomizer.nextDouble();

		int i = 0;
		for (double w : weights) {
			if (u < w)
				return i;

			u -= w;
			i++;
		}

		throw new RuntimeException("Can not sample from an empty collection.");
	}

	static public int sampleCategorical(Iterable<Double> weights) {
		double weightSum = sum(weights);
		double u = weightSum * Randomizer.nextDouble();

		int i = 0;
		for (double w : weights) {
			if (u < w)
				return i;

			u -= w;
			i++;
		}

		throw new RuntimeException("Can not sample from an empty collection.");
	}

    static public double[] sampleSubset(double[] samples, int subsetSize) {
        assert subsetSize <= samples.length;

        double[] rndSubset = new double[subsetSize];
        int i = 0;

        for (int sampleIdx : Randomizer.shuffled(samples.length)) {
            rndSubset[i] = samples[sampleIdx];
            i++;
            if (i >= subsetSize)
                break;
        }

        return rndSubset;
    }

    static public <T> ArrayList<T> sampleSubset(List<T> samples, int subsetSize) {
        assert subsetSize <= samples.size();

        ArrayList<T> rndSubset = new ArrayList<>(subsetSize);
        int i = 0;

        for (int sampleIdx : Randomizer.shuffled(samples.size())) {
            rndSubset.set(i, samples.get(sampleIdx));
            i++;
            if (i >= subsetSize)
                break;
        }

        return rndSubset;
    }

    static public double logAddExp(double logP1, double logP2) {
        double logPMin, logPMax;
        if (logP1 < logP2) {
            logPMin = logP1;
            logPMax = logP2;
        } else {
            logPMin = logP2;
            logPMax = logP1;
        }

        return logPMax + Math.log1p(Math.exp(logPMin - logPMax));
    }

    static public double logSumExp(List<Double> logProbs) {
        double logPMax = Double.NEGATIVE_INFINITY;

        for (double logP : logProbs)
            if (logP > logPMax)
                logPMax = logP;

        double sumExp = 0;
        for (double logP : logProbs)
            sumExp += Math.exp(logP - logPMax);

        return logPMax + Math.log(sumExp);
    }

    static public double[] list2array(List<Double> list) {
        double[] array = new double[list.size()];
        for (int i=0; i<list.size(); i++) {
            array[i] = list.get(i);
        }
        return array;
    }

    static public Node getOtherChild(Node parent, Node child) {
        for (Node c : parent.getChildren())
            if (c != child)
                return c;

        return null;
    }

    static public Node getSibling(Node node) {
        return getOtherChild(node.getParent(), node);
    }


    /**
     * Alternative to ´getClosestRelatives´ where all relatives are collected.
     * @param node
     * @param height
     * @return
     */
    static public HashSet<Node> getClosestRelatives(Node node, double height) {
        return getClosestRelatives(node, height, Integer.MAX_VALUE);
    }

    /**
     * Gather the ´nClosest´ nodes at the same height with share the most recent MRCA with ´node´.
     * @param node: the reference node.
     * @param height: the height at which to search (can be anywhere between node.height and node.parent.height)
     * @param nClosest: the number of closest relatives to collect.
     * @return
     */
    static public HashSet<Node> getClosestRelatives(Node node, double height, int nClosest) {
        HashSet<Node> relatives = new HashSet<>();

        if (node.isRoot())
            throw new IllegalArgumentException();

        Node ancestor = node;
        while ((relatives.size() < nClosest) && (!ancestor.isRoot())) {
            // If we haven't gathered enough relatives, go further up in the tree
            relatives.addAll(
                    getDescendantsAtHeight(getSibling(ancestor), height)
            );

            ancestor = ancestor.getParent();
        }

        return relatives;
    }

    static public HashSet<Node> getDescendantsAtHeight(Node node, double height) {
        HashSet<Node> descendants = new HashSet<>();

        // Make sure we are not too far down already
        if (!node.isRoot()) {
            assert node.getParent().getHeight() > height;
        }

        if (node.getHeight() <= height) {
            // Target height reached -> root is the only descendant to return
            descendants.add(node);
            return descendants;
        } else {
            // Still above target height -> call recursively on children
            for (Node child : node.getChildren()) {
                descendants.addAll(getDescendantsAtHeight(child, height));
            }
        }

        return descendants;
    }

    static public HashSet<Node> getLeaves(Node node) {
        HashSet<Node> leaves = new HashSet<>();

        if (node.isLeaf()) {
            leaves.add(node);
        } else {
            for (Node child : node.getChildren()) {
                leaves.addAll(getLeaves(child));
            }
        }

        return leaves;
    }

    static public Node mrca(Node a, Node b) {
        while (a != b) {
            if (a.getHeight() < b.getHeight())
                a = a.getParent();
            else
                b = b.getParent();
        }

        return a;
    }

    static public double phylogeneticDistance(Node a, Node b) {
        double nodeDistance = 0;
        double patristicDistance = 0;

        while (a != b) {
            // Move the younger node to its parent
            if (a.getHeight() < b.getHeight())
                a = a.getParent();
            else
                b = b.getParent();

            nodeDistance += 1;
        }

        return nodeDistance;
    }

    static public double phylogeneticDistancePatristic(Node a, Node b, double aHeight, double bHeight) {
        double patristicDistance = 0;

        while (a != b) {
            // Move the younger node to its parent
            if (aHeight < bHeight) {
                a = a.getParent();
                patristicDistance += a.getHeight() - aHeight;
                aHeight = a.getHeight();
            } else {
                b = b.getParent();
                patristicDistance += b.getHeight() - bHeight;
                bHeight = b.getHeight();
            }

        }

        return patristicDistance;
    }
}
