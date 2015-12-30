import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Vector;

public class PatternSearch {
	private SuffixTree tree;
	private MassTable massTable;
	private Vector<MatchPos> list1;
	private Vector<MatchPos> list2;
	private boolean first = true;
	public static double errorTolerance = 0.05;
	private int scaledErrorTolerance = 5;
	private ArrayList<Node> nodePool;
	private PeptideFileHandler pf;

	/**
	 * constructor, create two lists
	 * @param tree suffix tree
	 * @param mt mass table
	 * */
	public PatternSearch(SuffixTree tree, MassTable mt, PeptideFileHandler pf){
		this.tree = tree;
		massTable = mt;
		this.pf = pf;
		list1 = new Vector<MatchPos>();
		list2 = new Vector<MatchPos>();
		scaledErrorTolerance = (int)Math.round(errorTolerance*massTable.getScaleFactor());
		//System.out.println("search error:"+errorTolerance);
		initNodePool();
	}
	
	/**
	 * search a pattern against a suffix tree
	 * @param pattern 
	 * 			the given pattern consisting of several numbers
	 * @return all occurrences of the pattern in the database
	 * */
	public ArrayList<LeafEdge> searchPattern(Pattern pattern, int scanNum){
		
		first = true;
		if (pattern.length() <= 0)
			return null;
		for (int i=0; i<pattern.length(); i++){
			searchPattern(pattern.get(i));
			//no match found
			if (list2.size() == 0)
				return null;
			//swap list1 and list2
			Vector<MatchPos> tmpList = list1;
			list1 = list2;
			list2 = tmpList;
			list2.clear();
			first = false;
		}

		ArrayList<LeafEdge> leafSet = findLeafSet(scanNum);
		//Collections.sort(leafSet, new SuffixPositionComparator());
		
		list1.clear();
		list2.clear();
		
		
		return leafSet;
	}
	
	/**
	 * search the suffix tree with a given value
	 * @param mass 
	 * 			a given integer value
	 * @notes there may exist several strings the sum of whom is equal to "mass"
	 * */
	public void searchPattern(int mass){
		for (int i=-scaledErrorTolerance; i<=scaledErrorTolerance; i++){
			if ((mass+i)<=0 || (mass+i)>massTable.getMaxSize())
				continue;
			Iterator<String> ite = massTable.getMassString(mass+i);
			while (ite.hasNext()){
				searchPattern((String)ite.next());
			}	
		}
	}

	
	/**
	 * search a particular string against the suffix tree
	 * @param str 
	 * 			the given string to be searched
	 * */
	public void searchPattern(String str){
		if (first)
			searchPattern(str, null);
		else {
			for (int i=0; i<list1.size(); i++)
				searchPattern(str, list1.get(i));
		}
	}
	
	/**
	 * search a particular string from a given position of the suffix tree
	 * @param str 
	 * 			the string to be searched
	 * @param matchPos 
	 * 			a position in the suffix tree (start matching the string from this given position)
	 * @return the matched position, null if mismatch
	 * */
	public MatchPos searchPattern(String str, MatchPos matchPos){
		int i = 0, j = 0;
		Node node;
		Edge edge;
		MatchPos newMatchPos = null;
		
		//start matching from the root
		if (matchPos == null){
			node = tree.getRoot();
			edge = null;
		}
		else {//start matching from the middle of the tree
			edge = matchPos.getEdge();
			j = matchPos.getLength();
			int remaining = edge.getLength() - j;
			if (str.length() <= remaining) {
				for (int k=0; k<str.length(); k++){
					boolean result = compareChar(str.charAt(k), edge.getItemAt2(j+k));
					if (!result) return null;
				}
				newMatchPos = new MatchPos(edge, j+str.length());
				list2.add(newMatchPos);
				return newMatchPos;
			}
			else {
				for (int k=0; k<remaining; k++){
					boolean result = compareChar(str.charAt(k), edge.getItemAt2(j+k));
					if (!result) return null;
				}
				i = remaining;
				node = edge.getEndNode();
			}
		}
		
		while (i < str.length()){
			edge = node.findEdge(str.charAt(i));
			if (edge == null)
				return null;
			i++;
			for (int k=edge.getBeginIndex()+1; k<=edge.getEndIndex(); k++){
				if (i == str.length()){
					newMatchPos = new MatchPos(edge, k-edge.getBeginIndex());
					list2.add(newMatchPos);
					return newMatchPos;
				}
				boolean result = compareChar(str.charAt(i), edge.getItemAt(k));
				if (!result) return null;	
				i++;
			}
			node = edge.getEndNode();
		}
		
		newMatchPos = new MatchPos(edge, edge.getLength());
		list2.add(newMatchPos);
		return newMatchPos;
	}
	
	/**
	 * find all the leaves in the subtrees of the given matched positions
	 * @return ArrayList<SuffixPosition> 
	 * @notes SuffixPosition indicates on which sequence and offset in the sequence
	 * */
	public ArrayList<LeafEdge> findLeafSet(int scanNum){
		ArrayList<LeafEdge> leafSet = new ArrayList<LeafEdge>();
		for (int i=0; i<list1.size(); i++)
			findLeafSet(leafSet, list1.get(i), scanNum);
		
		return leafSet;
	}
	
	/**
	 * find all the leaves in the subtree of the given match position
	 * @param leafSet used to store the leaf found
	 * @param the given position in the suffix tree
	 * */
	public void findLeafSet(ArrayList<LeafEdge> leafSet, MatchPos matchPos, int scanNum){
		ArrayList<Node> stack = new ArrayList<Node>();
		Edge matchEdge = matchPos.getEdge();
		if (matchEdge == null)//mismatch
			return;
		
		Node node = matchEdge.getEndNode();
		if (node == null){
			if (inErrorTolerance((LeafEdge)matchEdge, scanNum))
				leafSet.add((LeafEdge)matchEdge);
			return;
		}
		
		if (node.getVisitStatus()) return;
		nodePool.add(node);
		node.setVisitStatus(true);
		if (inErrorTolerance(node, scanNum) == false) return;
		//nodePool.add(node);
		//node.setVisitStatus(true);
		stack.add(node);
		while (stack.size() > 0){
			ArrayList<Node> childNodes = new ArrayList<Node>();
			for (int k=0; k<stack.size(); k++){
				Node node2 = stack.get(k);
				for (int i=0; i<21; i++){
					Edge edge = node2.getEdge(i);
					if (edge == null)
						continue;
					Node endNode = edge.getEndNode();
					if (endNode != null){ 
						if (endNode.getVisitStatus() == false){
							if (inErrorTolerance(endNode, scanNum))
								childNodes.add(endNode);
							nodePool.add(endNode);
							endNode.setVisitStatus(true);
						}
					}
					else{
						if (inErrorTolerance((LeafEdge)edge, scanNum))
							leafSet.add(((LeafEdge)edge));
					}
				}
			}
			stack.clear();
			stack = null;
			stack = childNodes;
		}
	}
	
	/**
	 * sets the error tolerance when searching the suffix tree
	 * @param error
	 * 			the error
	 * */
	public static void setSearchErrorTolerance(double error){
		errorTolerance = error;
	}
	
	/**
	 * compare two characters, I and L are considered to be matched
	 * @param ch1 first char
	 * @param ch2 second char
	 * @return true if the two input characters match
	 * */
	private boolean compareChar(char ch1, char ch2){
		boolean result = false;
		if (ch1 == ch2)
			result = true;
		//else if ((ch1=='I' && ch2=='L') || (ch1=='L' && ch2=='I'))
		//	result = true;//'I' && 'L' have same mass
		
		return result;
	}
	
	public void initNodePool(){
		nodePool = new ArrayList<Node>();
	}
	
	public void clearNodePool(){
		for (int i=0; i<nodePool.size(); i++){
			nodePool.get(i).setVisitStatus(false);
		}
		//nodePool.removeAll(nodePool);
		nodePool.clear();
	}
	
	public boolean inErrorTolerance(LeafEdge edge, int scanNum){
		return edge.inErrorTolerance(pf.getMassValue(scanNum));
	}
	
	public boolean inErrorTolerance(Node node, int scanNum){
		return node.inErrorTolerance(pf.getMassValue(scanNum));
	}
}

/**
 * Indicates a position in the suffix tree
 * */
class MatchPos{
	private Edge edge;// an edge in the suffix tree
	private int length;//the number of characters that have been matched
	
	
	/**
	 *  create an empty object
	 */
	public MatchPos(){
		edge = null;
		length = 0;
	}
	
	/**
	 * create a new MatchPos object with given information
	 * */
	public MatchPos(Edge edge, int len){
		this.edge = edge;
		length = len;
	}
	
	/**
	 * get the matching edge
	 * @return the matching edge
	 * */
	public Edge getEdge(){
		return edge;
	}
	
	/**
	 * get the offset where the matching position lies in the edge
	 * @return the number of characters matching the edge
	 * */
	public int getLength(){
		return length;
	}
	
	
	public String toString(){
		return edge.toString()+" "+length;
	}
}