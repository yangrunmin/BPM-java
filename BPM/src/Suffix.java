

/**
 * Defines suffix.
 * <p/>
 * When a new tree is added to the table, we step through all the currently defined suffixes
 * from the active point to the end point.  This structure defines a <b>Suffix</b> by its final character.
 * In the canonical representation, we define that last character by starting at a node in the tree,
 * and following a string of characters, represented by <b>first_char_index</b> and <b>last_char_index</b>.
 * The two indices point into the input string.  Note that if a suffix ends at a node,
 * there are no additional characters needed to characterize its last character position.
 * When this is the case, we say the node is <b>explicit</b>,
 * and set <b>first_char_index > last_char_index<b> to flag that.
 */
public class Suffix {
    private Node originNode;
    private int beginIndex;
    private int endIndex;

    public Suffix(Node originNode, int beginIndex, int endIndex) {
        this.originNode = originNode;
        this.beginIndex = beginIndex;
        this.endIndex = endIndex;
    }

    /**
     * test if a node is explicit
     * */
    public boolean isExplicit() {
        return beginIndex > endIndex;
    }

    /**
     * test if a node is implicit
     * */
    public boolean isImplicit() {
        return endIndex >= beginIndex;
    }

    /**
     * match the suffix along the tree from Suffix.originNode
     * */
    public void canonize() {
        if (!isExplicit()) {
            Edge edge = originNode.findEdge(originNode.charAt(beginIndex));
            
            int edgeSpan = edge.getSpan();
            while (edgeSpan <= getSpan()) {
                beginIndex += edgeSpan + 1;
                originNode = edge.getEndNode();
                if (beginIndex <= endIndex) {
                    edge = edge.getEndNode().findEdge(originNode.charAt(beginIndex));
                    edgeSpan = edge.getSpan();
                }
            }
        }
    }

    /**
     * get the span
     * */
    public int getSpan() {
        return endIndex - beginIndex;
    }

    /**
     * get the origin node
     * */
    public Node getOriginNode() {
        return originNode;
    }

    /**
     * get the begin index of the suffix
     * */
    public int getBeginIndex() {
        return beginIndex;
    }
    
    /**
     * increase begin index by 1
     * */
    public void incBeginIndex() {
        beginIndex++;
    }

    /**
     * change the origin node of the suffix 
     * */
    public void changeOriginNode() {
        originNode = originNode.getSuffixNode();
    }

    /**
     * get the end index of the suffix
     * */
    public int getEndIndex() {
        return endIndex;
    }

    /**
     * increase the end index by 1
     * */
    public void incEndIndex() {
        endIndex++;
    }
}
