

import java.util.Iterator;
import java.util.Map;
import java.util.Collection;

/**
 * Each node contains several attributes: 
 * @attribute suffixNode: a suffix node from it
 * @edges: a set of edges starting from it
 * @name: label
 * @leafIndex: -1 for non-leaf; y if it is a leaf for suffix [y...n]
 * */
public class Node {
    private SuffixTree suffixTree;
    private Node suffixNode = null;
    private Edge edges[] = null;
    private boolean visited = false;
    private int minPeptideMass = 0;
    private int maxPeptideMass = 0;

    /**
     * constructor
     * @param node
     * 			a node
     * @param suffixNode
     * 			the suffix node out from this node
     * */
    public Node(Node node, Node suffixNode) {
        this(node.suffixTree, suffixNode);
    }

    /**
     * constructor
     * @param suffixTree
     * 			an instance of suffixTree
     * @param suffixNode
     * 			a node
     * */
    public Node(SuffixTree suffixTree, Node suffixNode) {
        this.suffixTree = suffixTree;
        this.suffixNode = suffixNode;
        edges = new Edge[21];
        for (int i=0; i<21; i++)
        	edges[i] = null;
    }

    /**
     * newly added, constructor
     * @param node
     * 			a node
     * */
    public Node(Node node){
    	this.suffixTree = node.suffixTree;
    	this.suffixNode = null;
    }
    
    /**
     * get the character at the indicated position
     * @param index
     * 			indicated position
     * */
    public char charAt(int index) {
        return suffixTree.getText().charAt(index);
    }

    /**
     * add an edge into the set of edges(HashMap) associated with the node
     * use the first character on the edge as the key, allowing fast searching
     * @param charIndex: the index of a (typically, first) character
     * @param edge: the edge to be added 
     * @note: all the edges out from a node start with a different character 
     * */
    public void addEdge(int charIndex, Edge edge) {
    	edges[getCharacterIndex(charAt(charIndex))] = edge;
    }

    /**
     * remove an edge from the HapMap associated with the node
     * @param charIndex: the index of the character, used to specify an edge
     * @note: all the edges out from a node start with a different character 
     * */
    public void removeEdge(int charIndex) {
    	edges[getCharacterIndex(charAt(charIndex))] = null;
    }
    
    public void removeEdge(char ch){
    	edges[getCharacterIndex(ch)] = null;
    }
    
    public void removeAll(){
    	for (int i=0; i<21; i++)
    		edges[i] = null;
    	edges = null;
    }

    /**
     * find an edge out from a node
     * @param ch
     * 			the first character on the edge
     * @return the edge
     * */
    public Edge findEdge(char ch) {
    	return edges[getCharacterIndex(ch)];
    }
    
    /**
     * find an edge
     * @param index
     * 			the given edge index
     * @return the required edge
     * */
    public Edge getEdge(int index){
    	if (index<0 || index>20){
    		System.out.println("invalid index");
    		System.exit(0);
    	}
    	
    	return edges[index];
    }

    /**
     * get the index corresponding to one of the 22 characters
     * @param ch
     * 			the given character
     * @return the index (from 0 to 21) associated with the character
     * */
    private int getCharacterIndex(char ch){
    	if (ch == '#')
    		return 19;
    	else if (ch == '$')
    		return 20;
    	else
    		return MassTable.aminoAcidIndex[ch-'A'];
    }
    
    
    /**
     * get the suffix node out from the node
     * @return the suffix node out from this node
     * */
    public Node getSuffixNode() {
        return suffixNode;
    }

    /**
     * set the suffix node out from a node
     * @param suffixNode: the suffix node to be set
     * */
    public void setSuffixNode(Node suffixNode) {
        this.suffixNode = suffixNode;
    }

    /**
     * newly added, test if a node has a suffix link from it
     * @return true for yes
     * */
    public boolean hasSuffixNode(){
    	return suffixNode!=null;
    }
    
    
    /**
     * newly added, get the suffix tree
     * @return suffixTree
     * */
    public SuffixTree getSuffixTree(){
    	return suffixTree;
    }
    
    
    public void setVisitStatus(boolean boolvalue){
    	visited = boolvalue;
    }
    
    public boolean getVisitStatus(){
    	return visited;
    }
    
    public void setMinPeptideMass(int mass){
    	minPeptideMass = mass;
    }
    
    public int getMinPeptideMass(){
    	return minPeptideMass;
    }
    
    public void setMaxPeptideMass(int mass){
    	maxPeptideMass = mass;
    }
    
    public int getMaxPeptideMass(){
    	return maxPeptideMass;
    }
    
    public boolean inErrorTolerance(int refMass){
    	int error = DoSearch2.getScaledMassTolerance();
    	if ((refMass<= maxPeptideMass) && (refMass>=minPeptideMass))
    		return true;
    	return false;
	}
    
  
    /*@Override
    public String toString() {
        return ((Integer) name).toString();
    }*/
}

