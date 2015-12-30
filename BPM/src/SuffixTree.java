/**
 * Refactored java-code originally based on Mark Nelson's C++ implementation of Ukkonen's algorithm.
 * http://illya-keeplearning.blogspot.com/search/label/suffix%20tree
 */


import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Queue;
import java.util.ArrayList;
import java.util.Vector;

public class SuffixTree {
    private String text;
    private Node root;
    private int seqIndex = 0;
    private int leafCreatedThisStep = 0;
    private ProteinDatabase database = null;
    
    public SuffixTree(String text) {
        this.text = text;
        root = new Node(this, null);
        Suffix active = new Suffix(root, 0, -1);
        for (int i = 0; i < text.length(); i++)
            addPrefix(active, i);
  
    }
  
    public SuffixTree(String text, ProteinDatabase database) {
        this.text = text;
        root = new Node(this, null);
        this.database = database;
        Suffix active = new Suffix(root, 0, -1);
        for (int i = 0; i < text.length(); i++)
        	addPrefix(active, i);
       
    }

    private void addPrefix(Suffix active, int endIndex) {
        Node lastParentNode = null;
        Node parentNode;

        while (true) {
            Edge edge;
            parentNode = active.getOriginNode();

            // Step 1 is to try and find a matching edge for the given node.
            // If a matching edge exists, we are done adding edges, so we break out of this big loop.
            if (active.isExplicit()) {
                edge = active.getOriginNode().findEdge(text.charAt(endIndex));
                if (edge != null)
                    break;
            } else {
                //implicit node, a little more complicated
                edge = active.getOriginNode().findEdge(text.charAt(active.getBeginIndex()));
                int span = active.getSpan();
                if (text.charAt(edge.getBeginIndex() + span + 1) == text.charAt(endIndex))
                    break;
                parentNode = edge.splitEdge(active);
            }

            // We didn't find a matching edge, so we create a new one, add it to the tree at the parent node position,
            // and insert it into the hash table.  When we create a new node, it also means we need to create
            // a suffix link to the new node from the last node we visited.
            //Edge newEdge = new Edge(endIndex, text.length() - 1, parentNode);
            LeafEdge newEdge = new LeafEdge(endIndex, text.length()-1, parentNode, null);
            newEdge.setLeafIndex(leafCreatedThisStep);//newEdge.createSuffixPosition(getSeqIndex()+1, getLeafCreatedThisStep()+1);
            newEdge.setSeqNum(seqIndex);
            leafCreatedThisStep++;
            determineSuffixPos();
            newEdge.insert();
            updateSuffixNode(lastParentNode, parentNode);
            lastParentNode = parentNode;

            // This final step is where we move to the next smaller suffix
            if (active.getOriginNode() == root)
                active.incBeginIndex();
            else
                active.changeOriginNode();
            active.canonize();
        }
        updateSuffixNode(lastParentNode, parentNode);
        active.incEndIndex(); 
        active.canonize();
    }

    private void updateSuffixNode(Node node, Node suffixNode) {
        if ((node != null) && (node != root)) {
            node.setSuffixNode(suffixNode);
        }
    }
    
    private void determineSuffixPos(){
    	if (leafCreatedThisStep == (database.getSeqLength(seqIndex)+1)){
    		seqIndex++;
    		leafCreatedThisStep = 0;
    	}
    }
    
    
    
    /**
     * newly added, get the root of a suffix tree
     * @return the root of the suffix tree
     * */
    public Node getRoot(){
    	return root;
    }

    /**
     * get the text associated with the suffix tree
     * */
    public String getText() {
        return text;
    }

    /**
     * gets the current sequence index
     * @return current sequence index
     * */
    public int getSeqIndex(){
    	return seqIndex;
    }
    
    
    /**
     * gets the number of leaf created in this step (for a string)
     * @return leaf created during this step
     * */
    public int getLeafCreatedThisStep(){
    	return leafCreatedThisStep;
    }
    
    /**
     * increases the number of leaf created in the current step
     * */
    public void increaseLeafCreated(){
    	leafCreatedThisStep++;
    }
   
    
    /**
	 * search a string in the suffix tree
	 * @param target
	 * 			the string to be searched
	 * @return none
	 */
	public void search(String target){
		ArrayList<SuffixPosition> startPosList = new ArrayList<SuffixPosition>();
		LinkedList<Node> stack = new LinkedList<Node>();
		Edge matchEdge = findMatchEdge(target);
		if (matchEdge == null){
			System.out.println("no mismatch!");
			return;
		}
		Node node = matchEdge.getEndNode();
		if (node == null){
			System.out.println("Details of matching:");
			startPosList.add(((LeafEdge)matchEdge).getSuffixPosition());
			System.out.println(startPosList);
			return;
		}
		
		stack.add(node);
		while (stack.size() > 0){
			LinkedList<Node> childNodes = new LinkedList<Node>();
			for (Node node2 : stack){
				for (int i=0; i<21; i++){
					Edge edge = node2.getEdge(i);
					if (edge == null)
						continue;
					if (edge.hasEndNode())
						childNodes.push(edge.getEndNode());
					else
						startPosList.add(((LeafEdge)edge).getSuffixPosition());
				}
			}
			
			stack = childNodes;
		}
		
		System.out.println("Details of matching:");
		//sort the result in increasing order of its seqNum
		Collections.sort(startPosList, new SuffixPositionComparator());
		System.out.println(startPosList);
	}
	
	/**
	 * find the matched edge for a string to be searched
	 * @param target
	 * 			the string to be searched
	 * @return null is "target" is not found in the suffix tree
	 * 		   the edge
	 */
	public Edge findMatchEdge(String target){
		if (target.length() == 0)
			return null;
		
		int i = 0;
		Node node = root;
		Edge edge = null;
		while (i < target.length()){
			edge = node.findEdge(target.charAt(i));
			if (edge == null)
				return null; 
			i++;
			for (int j=edge.getBeginIndex()+1; j<=edge.getEndIndex(); j++){
				if (i == target.length())
					return edge;
				String item = ((Character)edge.getItemAt(j)).toString();
				char chr = target.charAt(i);
				String tmp = String.valueOf(chr);
				if (!tmp.equals(item))
					return null;
				i++;
			}
			
			node = edge.getEndNode();
		}
		return edge;	
	}
	
}
