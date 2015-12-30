
public class LeafEdge extends Edge{
	
	 private int leafIndex = -1;
	 private int seqNum = -1;
	 private int peptideMass = 0;
	 private String pepString;
	 
	
	public LeafEdge(int beginIndex, int endIndex, Node startNode) {
       super(beginIndex, endIndex, startNode);
    }
    
    public LeafEdge(int beginIndex, int endIndex, Node startNode, Node endNode){
    	super(beginIndex, endIndex, startNode, endNode);
    }
    
    
    /**
     * newly added, get leaf index of an edge
     * @return leafIndex, -1 for a non-leaf edge
     * */
    public int getLeafIndex(){
    	return leafIndex;
    }

    /**
     * newly added, set the leaf index of an edge
     * @param index: the value to be set
     * */
    public void setLeafIndex(int index){
    	leafIndex = index;
    }
    
    /**
     * get the seqNum associated with an edge
     * @return the sequence number associated with the edge
     * */
    public int getSeqNum(){
    	return seqNum;
    }
    
    /**
     * set seqNum to be given value
     * @param num: value to be set
     * */
    public void setSeqNum(int num){
    	seqNum = num;
    }
    
    /**
	 * gets the suffixPosition associated with an edge
	 * */
	public SuffixPosition getSuffixPosition(){
		if (leafIndex == -1)
			return null;
		return new SuffixPosition(seqNum, leafIndex);
	}
	
	public void setPeptideMass(int mass){
		peptideMass = mass;
	}
	
	public int getPeptideMass(){
		return peptideMass;
	}
	
	public void setPeptideSequence(String str){
		pepString = str;
	}
	
	public String getPepString(){
		return pepString;
	}
	
	public boolean inErrorTolerance(int refMass){
		if (Math.abs(peptideMass-refMass) <= DoSearch2.getScaledMassTolerance())
			return true;
		return false;
	}

}
