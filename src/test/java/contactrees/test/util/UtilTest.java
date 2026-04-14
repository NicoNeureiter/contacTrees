package contactrees.test.util;

import java.util.HashSet;

import org.junit.Test;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import contactrees.test.ContactreesTest;
import contactrees.util.Util;


public class UtilTest extends ContactreesTest {

    @Test
    public void test() {
        Tree tree = acg2;
        Node me = node2_1;
        double height = 0.5;

        HashSet<Node> expected = new HashSet<>();
        HashSet<Node> actual;

        actual = Util.getClosestRelatives(me, height, 1);
        expected.add(node2_2);
        assert actual.equals(expected);

        actual = Util.getClosestRelatives(me, height, 2);
        expected.add(node2_3);
        assert actual.equals(expected);

        actual = Util.getClosestRelatives(me, height, 3);
        expected.add(node2_4);
        expected.add(node2_5);
        assert actual.equals(expected);

        actual = Util.getClosestRelatives(me, height, 4);
        assert actual.equals(expected);

        actual = Util.getClosestRelatives(me, height);
        assert actual.equals(expected);
    }

}
