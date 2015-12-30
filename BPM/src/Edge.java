

/**
 * Represents a string of characters.
 * <p/>
 * The suffix tree is made up of edges connecting nodes. Each edge represents a string of characters starting
 * at first_char_index and ending at last_char_index. 
 */
public class Edge {
    protected int beginIndex;
    protected int endIndex;
    protected Node startNode = null;
    protected Node endNode = null;
 
    public Edge(int beginIndex, int endIndex, Node startNode) {
        this.beginIndex = beginIndex;
        this.endIndex = endIndex;
        this.startNode = startNode;
        this.endNode = new Node(startNode, null);
    }
    
    public Edge(int beginIndex, int endIndex, Node startNode, Node endNode){
    	 this.beginIndex = beginIndex;
         this.endIndex = endIndex;
         this.startNode = startNode;
         this.endNode = endNode;
    }

    /**
     * When a suffix ends on an implicit node, adding a new character means I have to split an existing edge.
     * This function is called to split an edge at the point defined by the Suffix argument.
     * The existing edge loses its parent, as well as some of its leading characters.
     * The newly created edge descends from the original parent, and now has the existing edge as a child.
     * <p/>
     * Since the existing edge is getting a new parent and starting character,
     * its hash table entry will no longer be valid.  That's why it gets removed at the start of the function.
     * After the parent and start char have been recalculated, it is re-inserted.
     * <p/>
     * The number of characters stolen from the original node and given to the new node is equal to the number
     * of characters in the suffix argument, which is last - first + 1;
     *
     * @param s
     * @return
     */
    public Node splitEdge(Suffix suffix){
    	remove();
    	Node breakNode = new Node(suffix.getOriginNode(), null);
    	Edge newEdge = new Edge(beginIndex, beginIndex + suffix.getSpan(), suffix.getOriginNode(), breakNode);
        newEdge.insert();
        breakNode.setSuffixNode(suffix.getOriginNode());
        beginIndex += suffix.getSpan() + 1;
        startNode = breakNode;
        insert();
        return breakNode;
    }

    /**
     * insert the edge to the list associated with the startNode
     * Each (internal) node maintains a set of edges starting from it 
     * */
    public void insert() {
        startNode.addEdge(beginIndex, this);
    }

    /**
     * remove the edge from the set of edges associated with the startNode
     * */
    public void remove() {
        startNode.removeEdge(beginIndex);
    }

    /**
     * get the span of an edge 
     * @return the span of an edge
     * */
    public int getSpan() {
        return endIndex - beginIndex;
    }

    /**
     * get the length of an edge
     * @return the number of characters on an edge
     * */
    public int getLength(){
    	return endIndex - beginIndex + 1;
    }
    
    /**
     * get the begin index of an edge
     * @return the begin index of the edge
     * */
    public int getBeginIndex() {
        return beginIndex;
    }

    /**
     * get the end index of an edge
     * @return the end index of the edge
     * */
    public int getEndIndex() {
        return endIndex;
    }

    /**
     * set the end index of an edge
     * @param endIndex: the value to be set
     * */
    public void setEndIndex(int endIndex) {
        this.endIndex = endIndex;
    }

    /**
     * get the startNode (incoming node) of an edge
     * @return the start node of the edge
     * */
    public Node getStartNode() {
        return startNode;
    }

    /**
     * set the startNode of an edge
     * @param startNode: the value to be set
     * */
    public void setStartNode(Node startNode) {
        this.startNode = startNode;
    }

    /**
     * get the endNode of an edge
     * @return the end node of the edge
     * */
    public Node getEndNode() {
        return endNode;
    }
    
    
    /**
     * newly added, test if an edge has an edge node (internal node if yes)
     * @return true if yes
     * */
    public boolean hasEndNode(){
    	return endNode != null;
    }  
   
    /**
     * newly added, get an element on the edge
     * @param j: index, [beginIndex, endIndex]
     * @return the element
     * */
    public char getItemAt(int j){
    	String text = startNode.getSuffixTree().getText();
    	return text.charAt(j);
    }
    
    /**
     * get an element on the edge
     * @param j: an offset ranging from 0 to length-1
     * @return the element
     * */
    public char getItemAt2(int j){
    	String text = startNode.getSuffixTree().getText();
    	return text.charAt(beginIndex + j);
    }
    
    /*@Override
    public String toString() {
        return endNode.toString();
    }*/
    
    public String toString(){
    	return startNode.getSuffixTree().getText().substring(beginIndex, endIndex+1);
    }
}