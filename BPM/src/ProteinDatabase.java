import java.io.IOException;
import java.io.PrintWriter;
import java.util.Random;
import java.util.Vector;


/**
 * Defines a protein database
 * 
 * size: number of strings contained in the database
 * sequence: the combined resulting string, a # is inserted between two strings, a '$' is appended to the end of the sequence 
 * seqLength: array used to store the length of each string
 * accumulatedLength: array used to store the accumulated length of strings
 * 
 * */
public class ProteinDatabase {
	
	private int size;
	private String sequence;//combined sequence
	private Vector<String> proteinID;
	private Vector<String> seqSet;//store each individual sequence
	private int[][] peptidePos;
	
	/**
	 * constructor
	 * */
	public ProteinDatabase(){
		size = 0;
		sequence = "";
		proteinID = new Vector<String>();
		seqSet = new Vector<String>();
	}
	
	/**
	 * gets the number of strings in the database
	 * @return number of strings 
	 * */
	public int size(){
		return seqSet.size();
	}
	
	/**
	 * get the length of a particular string
	 * @param seqIndex: the index of the query string
	 * @return the length of the query string
	 * */
	public int getSeqLength(int seqIndex){
		return seqSet.get(seqIndex).length();
	}
	
	/**
	 * gets an individual protein sequence
	 * @param index 
	 * 			the protein index, starts from 0
	 * @return the protein sequence
	 * */
	public String getIndividualSeq(int index){
		return seqSet.get(index);
	}
	
	/**
	 * add an individual sequence to the sequence set
	 * @param seq
	 * 			the sequence to be added
	 * */
	public void addIndividualSeq(String seq){
		seqSet.add(seq);
	}
	
	/**
	 * gets the combined sequence
	 * @return the combined sequence
	 * */
	public String getSequence(){
		return sequence;
	}
	
	/**
	 * set the "sequence"
	 * @param seq: the value to be set
	 * */
	public void setSequence(String seq){
		sequence = seq;
	}
	
	
	/**
	 * gets the protein name
	 * @param index
	 * 			the given protein index
	 * @return the protein name 
	 * */
	public String getProteinID(int index){
		return proteinID.get(index);
	}
	
	/**
	 * stores a protein name
	 * @param proteinName
	 * 			the given protein name, may contain several parts delimited by spaces
	 * @notes only keep the first part in proteinName
	 * */
	public void addProteinID(String proteinName){
		String delims = "[, ]+";
		String[] res = proteinName.split(delims);
		proteinID.add(res[0].substring(1));//only keep the first part (without ">" symbol) in proteinName
		
		for (int i=0; i<res.length; i++)
			res[i] = null;
		res = null;
		delims = null;
	}
	
	/**
	 * gets the subsequence indicated by [startPos, endPos) from a sequence 
	 * */
	public String getSubSeq(int seqNum, int startPos, int endPos){
		//String tmp = seqSet.get(seqNum);
		return seqSet.get(seqNum).substring(startPos, endPos);
	}
	
	
	/**
	 * Test function,given a query string, find its occurrences in the protein database
	 * @param queryString
	 * 			the query string
	 * @note output each occurrence in the manner <seqNum, offset, proteinID>
	 * 		 BOTH seqNum and offset index starts from 0, instead of 1
	 * */
	public void findAllOccurrence(String queryString){
		System.out.println("All occurrences of "+queryString);
		if (queryString.length() == 0)
			return;
		for (int i=0; i<this.size(); i++){
			String currentSeq = this.getIndividualSeq(i);
			for (int j=0; j<currentSeq.length(); j++){
				String tmp = currentSeq.substring(j);
				int pos = tmp.indexOf(queryString);
				if (pos == -1)
					break;
				System.out.print("<"+i+","+(j+pos)+","+proteinID.get(i)+">, ");
				j += pos;
			}
		}
		System.out.println();
	}
	
	/**
	 * gets a random subsequence of length <subSeqLength> from the protein database
	 * @param subSeqLength
	 * 			the length of the subsequence
	 * @param random
	 * 			an instance of java.util.Random, used to create random numbers
	 * @return a random subsequence
	 * */
	public String getRandomSeq(int subSeqLength, Random random){
		//choose a random sequence with length > subSeqLength
		int seqNum = random.nextInt(this.size());
		int seqLength = getSeqLength(seqNum);
		while (seqLength <= subSeqLength){
			seqNum = random.nextInt(this.size());
			seqLength = getSeqLength(seqNum);
		}
		
		//choose a random position within the selected sequence
		int startPos = random.nextInt(seqLength - subSeqLength);
		String subseq = getSubSeq(seqNum, startPos, startPos+subSeqLength);
		//output the position of the subsequence <seqNum, startPos>, both indexes starts from 0
		System.out.println("<"+seqNum+","+(startPos)+">  "+ subseq);
		return subseq;
	}
	
	
	/**
	 * extract peptide info from sequences based on Trypsin digest rule
	 * @return the peptide info contained in each sequence in the database
	 * */
	public Vector<Vector<Integer>> extractPeptideFromSeq() throws IOException{
		Vector<Vector<Integer>> peptideTable = new Vector<Vector<Integer>>(this.size());
		for (int i=0; i<this.size(); i++){
			Vector<Integer> peptideInfo = extractPeptideFromOneSeq(this.getIndividualSeq(i));
			peptideTable.add(peptideInfo);
		}
			
		return peptideTable;
	}
	
	/**
	 * extract peptides from a sequence based on Trypsin digest rule (after R or K, but not before P)
	 * @param seq
	 * 			a given protein sequence
	 * @return the list of peptides in the sequence, indicated by their start position in the sequence
	 * */
	public Vector<Integer> extractPeptideFromOneSeq(String seq){
		Vector<Integer> peptideInfo = new Vector<Integer>();
		peptideInfo.add(0);
		for (int i=0; i<seq.length()-1; i++){
			char current = seq.charAt(i);
			char next = seq.charAt(i+1);
			//if ((current=='R'||current=='K') && next!='P')
				//peptideInfo.add(i+1);
			if (current=='R'||current=='K')
				peptideInfo.add(i+1);
		}
		
		return peptideInfo;
	}
	
	/**
	 * extracts peptides from the database. 
	 * @return a 2-dimention array. int[i][j]: the peptide index for the j-th position in seqence-i. The exact 
	 * 			peptide start position can be found by querying peptideTable with this int[i][j] value.
	 * @notes:  extractPeptideFromSeq() vs extractPeptideFromSeq2(): the former only stores the start position 
	 * of each peptide; the latter stores the start position of a peptide for each given position in the sequence,
	 * computed only once, used to speed up the program.
	 * */
	public Vector<Vector<Integer>> extractPeptideFromSeq2() throws IOException{
		initPeptidePos();
		Vector<Vector<Integer>> peptideTable = new Vector<Vector<Integer>>(this.size());
		for (int i=0; i<this.size(); i++){
			Vector<Integer> peptideInfo = extractPeptideFromOneSeq2(i);
			peptideTable.add(peptideInfo);
		}
			
		return peptideTable;
	}
	
	/**
	 * extracts peptide information from one sequence
	 * @param seqIndex
	 *			sequence index
	 * @return the peptide information of the sequence 
	 * */
	public Vector<Integer> extractPeptideFromOneSeq2(int seqIndex){
		String seq = getIndividualSeq(seqIndex);
		Vector<Integer> peptideInfo = new Vector<Integer>();
		peptideInfo.add(0);
		for (int i=0; i<seq.length()-1; i++){
			peptidePos[seqIndex][i] = peptideInfo.size()-1;
			char current = seq.charAt(i);
			//char next = seq.charAt(i+1);
			//if ((current=='R'||current=='K') && next!='P')
				//peptideInfo.add(i+1);
			if (current=='R'||current=='K')
				peptideInfo.add(i+1);
		}
		
		peptidePos[seqIndex][seq.length()-1] = peptideInfo.size()-1;
		
		return peptideInfo;
	}
	
	/**
	 * gets the index of the corresponding peptide
	 * @param seqIndex
	 * 			sequence index
	 * @param offset
	 * 			a position in the sequence
	 * @return the index of the corresponding peptide
	 * */
	public int getPeptidePos(int seqIndex, int offset){
		return peptidePos[seqIndex][offset];
	}
	
	/**
	 * initialization
	 * */
	private void initPeptidePos(){
		int numSeq = size();
		peptidePos = new int[numSeq][];
		for (int i=0; i<numSeq; i++)
			peptidePos[i] = new int[getSeqLength(i)];
	}
	
	
	/**
	 * Test function, writes peptideTable to file
	 * @param peptideTable
	 * 			obtained peptideTable by calling extractPeptideFromSeq()
	 * @param outputFile
	 * 			output file name
	 * */
	public void writePeptideTable2File(Vector<Vector<Integer>> peptideTable, String outputFile) throws IOException{
		PrintWriter writer = new PrintWriter(outputFile);
		for (int i=0; i<peptideTable.size(); i++){
			writer.print(i+": ");
			writer.print(peptideTable.get(i));
			writer.println();
		}
		writer.close();
	}
	
	public void freeMem(Vector<Vector<Integer>> peptideTable){
		for (int i=0; i<peptideTable.size(); i++){
			peptideTable.get(i).clear();
		}
		peptideTable.clear();
		peptideTable = null;
		for (int i=0; i<size(); i++)
			peptidePos[i] = null;
		peptidePos = null;
	}
	
	
	public static void main(String[] args) throws IOException{
		DatabaseFileHandler df = new DatabaseFileHandler();
		ProteinDatabase pd = df.loadDatabase("test-index2.fasta");
		Vector<Vector<Integer>> peptideTable = pd.extractPeptideFromSeq2();
		System.out.println(peptideTable.get(0));
		//for (int i=245; i<245; i++)
			System.out.print(pd.getPeptidePos(0, 245)+" ");
	}
}

