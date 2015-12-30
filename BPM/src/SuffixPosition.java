import java.io.IOException;
import java.util.Vector;

/**
 * For each leaf in the generalized suffix tree, represents its position in the original string
 * 
 * seqNum: sequence (string) number where the suffix corresponding to a leaf originates from, sequence number begins with 1
 * posInSeq: the start position of the suffix in the sequence, ranging from 1,...,seqLength 
 * 
 */
public class SuffixPosition {
	private int seqNum;//seq index starts from 0
	private int posInSeq;
	private int peptideStartPos = 0;//index starts from 0
	private int peptideEndPos = 0;//index starts from 0, peptide ranges from [peptideStartPos, peptideEndPos)
	
	/**
	 * Creates a new suffixPosition, called each time a leaf is created
	 * @param seq
	 * 			a given the sequence number
	 *  @param pos
	 *  		the given position
	 * */
	public SuffixPosition(int seq, int pos){
		seqNum = seq;
		posInSeq = pos;
	}
	
	/**
	 * constructor
	 * */
	public SuffixPosition(int seq, int pos, int pepStart){
		seqNum = seq;
		posInSeq = pos;
		peptideStartPos = pepStart;
	}
	
	/**
	 * get "seqNum"
	 * @return seqNum
	 * */
	public int getSeqNum(){
		return seqNum;
	}
	
	/**
	 * get "posInSeq"
	 * @return posInSeq
	 * */
	public int getPosInSeq(){
		return posInSeq;
	}
	
	/**
	 * set the "seqNum" and "posInSeq" to given values
	 * @param seq
	 * 			given value which "seqNum" should be set to be
	 * @param pos
	 *  		given value which "posInSeq" should be set to be
	 */
	public void setSuffixPos(int seq, int pos){
		seqNum = seq;
		posInSeq = pos;
	}
	
	/**
	 * set "seqNum" to given value
	 * @param seq
	 * 			given value which "seqNum" should be set to be
	 */
	public void setSeqNum(int seq){
		seqNum = seq;
	}
	
	/**
	 * set "posInSeq" to given value
	 * @param pos
	 * 			given value which "posInSeq" should be set to be 
	 */
	public void setPosInSeq(int pos){
		posInSeq = pos;
	}
	
	/**
	 * set peptideStartPos to given value
	 * @param pepStart
	 * 			a given integer value
	 * */
	public void setPeptideStartPos(int pepStart){
		peptideStartPos = pepStart;
	}
	
	/**
	 * get peptideStartPos
	 * @return peptide start position
	 * */
	public int getPeptideStartPos(){
		return peptideStartPos;
	}
	
	/**
	 * get peptideEndPos
	 * @return get peptide end position
	 * */
	public int getPeptideEndPos(){
		return peptideEndPos;
	}
	
	/**
	 * find to which peptide should this suffixPos belong
	 * @param peptideTable
	 * 			peptide table
	 * @param proteinDatabase
	 * 			protein database
	 * @return peptide start position
	 * */
	/*public int compPeptidePos(Vector<Vector<Integer>> peptideTable, ProteinDatabase proteinDatabase){
		peptideStartPos = 0;
		Vector<Integer> pep = peptideTable.get(seqNum);
		for (int i=0; i<pep.size()-1; i++){
			if (posInSeq>= pep.get(i) && posInSeq<pep.get(i+1)){
				peptideStartPos = pep.get(i);
				peptideEndPos = pep.get(i+1);
				return peptideStartPos;
			}
		}
		
		peptideStartPos = pep.get(pep.size()-1);//the last peptide
		peptideEndPos = proteinDatabase.getSeqLength(seqNum);
		return peptideStartPos;
	}*/
	
	public void compPeptidePos(Vector<Vector<Integer>> peptideTable, ProteinDatabase proteinDatabase){
		peptideStartPos = 0;
		Vector<Integer> pep = peptideTable.get(seqNum);
		if (posInSeq >= pep.get(pep.size()-1)){
			peptideStartPos = pep.get(pep.size()-1);//the last peptide
			peptideEndPos = proteinDatabase.getSeqLength(seqNum);
		}
		else {
			int middle = this.binarySearch(posInSeq, pep);
			int tmp = pep.get(middle);
			if (tmp <= posInSeq){
				peptideStartPos = tmp;
				peptideEndPos = pep.get(middle+1);
			}
			else{
				peptideStartPos = pep.get(middle-1);
				peptideEndPos = tmp;
			}
		}
	}
	
	
	public int binarySearch(int num, Vector<Integer> pep){
		int low = 0; 
		int high = pep.size()-1;
		int middle = 0;
		while (low <= high){
			middle = (low+high)/2;
			if (pep.get(middle) == num)
				return middle;
			else if (pep.get(middle) < num)
				low = middle + 1;
			else
				high = middle - 1;
		}
		
		return middle;
	}
	
	public void compPeptidePos2(Vector<Vector<Integer>> peptideTable, ProteinDatabase proteinDatabase){
		peptideStartPos = 0;
		Vector<Integer> pep = peptideTable.get(seqNum);
		if (posInSeq >= pep.get(pep.size()-1)){
			peptideStartPos = pep.get(pep.size()-1);//the last peptide
			peptideEndPos = proteinDatabase.getSeqLength(seqNum);
		}
		else{
			int tmp = proteinDatabase.getPeptidePos(seqNum, posInSeq);
			peptideStartPos = pep.get(tmp);
			peptideEndPos = pep.get(tmp+1);
		}
		
	}
	
	
	/**
	 * get the peptide sequence
	 * @param pd
	 * 			the protein database
	 * @return the peptide sequence
	 * */
	public String getPeptideSequence(ProteinDatabase pd){
		return pd.getSubSeq(seqNum, peptideStartPos, peptideEndPos);
	}
	

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + peptideStartPos;
		result = prime * result + seqNum;
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		SuffixPosition other = (SuffixPosition) obj;
		if (peptideStartPos != other.peptideStartPos)
			return false;
		if (seqNum != other.seqNum)
			return false;
		return true;
	}

	public String toString(){
		return "<"+seqNum+","+posInSeq+"> ";
	}
	
	public static void main(String[] args) throws IOException{
		DatabaseFileHandler df = new DatabaseFileHandler();
		ProteinDatabase pd = df.loadDatabase("test.fasta");
		Vector<Vector<Integer>> peptideTable = pd.extractPeptideFromSeq();
		System.out.println(peptideTable.get(0));
		SuffixPosition sp = new SuffixPosition(0, 245);
		sp.compPeptidePos2(peptideTable, pd);
		System.out.println("start:"+sp.getPeptideStartPos()+"  end:"+sp.getPeptideEndPos());
		System.out.println("pepSeq:"+sp.getPeptideSequence(pd));
		SuffixPosition sp2 = new SuffixPosition(4,245);
		sp2.compPeptidePos(peptideTable, pd);
		if (sp.equals(sp2)) System.out.println("equal");
		
	}
}
