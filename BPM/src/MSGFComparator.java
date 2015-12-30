
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.HashSet;
import java.util.Vector;


public class MSGFComparator {
	
	private MyResultHandler my;
	
	/**
	 * constructor
	 * */
	public MSGFComparator(){
		my = new MyResultHandler();
	}
	
	/**
	 * loads msgf results
	 * @param msgfFile
	 * 			the file output by msgf+
	 * @param numScan
	 * 			total number of scan
	 * @param useFilter
	 * 			<true> if remove results with K/R in the middle; false otherwise
	 * */
	public Vector<Vector<String>> loadMsgfResult(String msgfFile, int numScan, boolean useFilter) throws IOException{
		FileInputStream msgf = new FileInputStream(msgfFile);
		BufferedReader msgfReader = new BufferedReader(new InputStreamReader(msgf));
		
		Vector<Vector<String>> msgfResult = new Vector<Vector<String>>();
		for (int i=0; i<=numScan; i++){
			Vector<String> re = new Vector<String>();
			msgfResult.add(re);
		}
		
		String delims = "[\t]+";
		String line = msgfReader.readLine();//head line
		line = msgfReader.readLine();
		while (line != null){
			String[] re = line.split(delims);
			int scanNum = Integer.parseInt(re[1]);
			String peptide = handleUnknownCharacter(re[9]);
			if (useFilter == false)
				msgfResult.get(scanNum).add(peptide);
			else{
				String tmp = peptide.substring(0, peptide.length()-1);
				if (tmp.indexOf('K')==-1 && tmp.indexOf('R')==-1)
					msgfResult.get(scanNum).add(peptide);
			}
	
			line = msgfReader.readLine();
		}
		
		msgfReader.close();
		msgf.close();
		return msgfResult;
	}
	
	
	/**
	 * removes some of the annotation string from the peptide
	 * @param peptide
	 * 			the given peptide sequence
	 * @return the result peptide 
	 * */
	private String handleUnknownCharacter(String peptide){
		peptide = peptide.replace("+57.021", "");
		return peptide;
	}
	
	
	/**
	 * determines whether a peptide contains K/R in the middle (including the first position).
	 * @param peptide
	 * 			a given peptide string
	 * @return true if it contains K/R in the peptide except the last position
	 * Notes: only the last character can be K/R
	 * */
	private boolean containKRInMiddle(String peptide){
		
		String tmp = peptide.substring(0, peptide.length()-1);
		int i = tmp.indexOf('K');
		int j = tmp.indexOf('R');
		if (i==tmp.length()-1 || j==tmp.length()-1)
			return true;
		if (i!=-1 && tmp.charAt(i+1)!='P')
			return true;
		if (j!=-1 && tmp.charAt(j+1)!='P')
			return true;
		
		return false;
	}

	/**
	 * gets some of the statistics of the MSGF result, such as the number of scans with search results
	 * @param msgfFile
	 * 			the result file output by MSGF+
	 * @param numScan
	 * 			total number of scans
	 * */
	public void getMSGFStatistics(String msgfFile, int numScan) throws IOException{
		Vector<Vector<String>> resultWithoutFilter = loadMsgfResult(msgfFile, numScan, false);
		Vector<Vector<String>> resultWithFilter = loadMsgfResult(msgfFile, numScan, true);
		
		int numScanWithResults = getScansWithResult(resultWithoutFilter, numScan);
		int numScanWithResultsFilter = getScansWithResult(resultWithFilter, numScan);

		
		System.out.println("**********MSGF result statistics**********");
		System.out.println("#scans with results without filtering:"+ numScanWithResults);
		System.out.println("#scans with results with filtering:"+ numScanWithResultsFilter);
		System.out.println("**********MSGF result statistics end**********");
		
	}
	
	/**
	 * gets the number of scans with at least one search result
	 * @param result
	 * 			the obtained search result
	 * @param numScan
	 * 			total number of scans
	 * @return the number of scans with at least one search result
	 * */
	private int getScansWithResult(Vector<Vector<String>> result, int numScan){
		int count = 0;
		for (int i=0; i<=numScan; i++)
			if (result.get(i).size() > 0)
				count++;
		return count;
	}
	

	/**
	 * compares our results with that of MSGF
	 * @param msgfFile
	 * 			the file output by MSGF+
	 * @param myFile
	 * 			the file output by my program
	 * @param matchFile
	 * 			the file name for storing the scan numbers where our result and that of MSGF+ are the same
	 * @param mismatchFile
	 * 			the file name for storing the scan numbers where our result and that of MSGF+ are different
	 * @param numScan
	 * 			the total number of scan
	 * */
	public void compareResults(String msgfFile, String myFile, String matchFile, String mismatchFile, int numScan, String peptideFile) throws IOException{	
		Vector<Vector<String>> msgfResult = loadMsgfResult(msgfFile, numScan, true);
		Vector<Vector<String>> myResult  = my.loadMyResult(myFile, numScan);
		PrintWriter matchWriter = new PrintWriter(matchFile);
		PrintWriter mismatchWriter = new PrintWriter(mismatchFile);
		
		boolean match = false;
		int num_match = 0;
		int overlap = 0;
		for (int i=0; i<=numScan; i++){
			if (msgfResult.get(i).size() == 0)
				continue;
			if (myResult.get(i).size() == 0)
				continue;
			overlap++;
			match = false;
			for (int j=0; j<myResult.get(i).size(); j++){
				String tmp1 = myResult.get(i).get(j);
				for (int k=0; k<msgfResult.get(i).size(); k++){
					String tmp2 = msgfResult.get(i).get(k);
					if (tmp1.equals(tmp2))
						match = true;
				}
			}
			if (match){
				num_match++;
				matchWriter.println(i);
			}
			else
				mismatchWriter.println(i);
		}
		
		double matchRate = num_match*1.0/getScansWithResult(msgfResult, numScan);
		System.out.println("Comparison result between ours and MSGF+:");
		System.out.println("#scans with search results in MSGF+:"+getScansWithResult(msgfResult, numScan));
		System.out.println("#scans with search results in ours+:"+getScansWithResult(myResult, numScan));
		System.out.println("Total overlap:"+overlap+"  #matches:"+num_match);
		System.out.println("match rate:"+matchRate);
		
		//getComparisionDetails(peptideFile, msgfResult, myResult);
		
		matchWriter.close();
		mismatchWriter.close();
	}
	
	
	public void getComparisionDetails(String peptideFile, Vector<Vector<String>> msgfResult, Vector<Vector<String>> myResult) throws IOException{
		MassTable mt = new MassTable();
		PeptideFileHandler pf = new PeptideFileHandler(mt);
		HashSet[] tagSet = pf.extractPatterns(peptideFile);
		int numScan = pf.getNumScan();
		PrintWriter writer_noTag = new PrintWriter("noTag");
		PrintWriter writer_msgfOnly = new PrintWriter("msgfOnly");
		
		
		for (int i=0; i<=numScan; i++){
			if (msgfResult.get(i).size()>0 && myResult.get(i).size()==0){
				if (tagSet[i+1].size() == 0)//in tagSet, index starts from 1
					writer_noTag.println(i);
				else//has tag but no search result
					writer_msgfOnly.println(i);
			}
		}
		
		writer_noTag.close();
		writer_msgfOnly.close();
	}
	
	
	
	/**
	 * args[0]: MSGF+ output result
	 * args[1]: Our output result
	 */
	
	public static void main(String[] args) throws IOException{
		MSGFComparator mh = new MSGFComparator();
		mh.compareResults(args[0], args[1], "matchFile", "mismatchFile", 38000, "deng_fei_data");
	}
	
	
	
	/*	
	public static void main(String[] args) throws IOException{
		
		MSGFComparator mh = new MSGFComparator();
		//mh.removeComplexPeptide("HCD-1-20.tsv", "HCD-1-20-simple");
		//String[] msgfResult = mh.loadMsgfResult("HCD_simple", 38000);
		//Vector<Vector<String>> myResult = mh.loadMyResult("output-1-15", 38000);
		mh.compareResults("MSGF-1-29-Q.tsv", "output-2-6", "matchFile", "mismatchFile", 38000, "deng_fei_data");
		//mh.getMSGFStatistics("HCD-1-29-PQ.tsv", 38000);
		//System.out.println("hello world!");
		System.out.println("\n\n");
		
	}*/

}
