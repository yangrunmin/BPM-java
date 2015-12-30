
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.HashSet;
import java.util.Vector;


/**
 * This project tries to compare X!Tandem output with ours
 * Note that the X!Tandem output has been converted into a text file by using "Peptide Hit Results Processor"
 * */
public class TandemComparator {
	
	private MyResultHandler my;
	
	public TandemComparator(){
		my = new MyResultHandler();
	}
	
	/**
	 * loads tandem output results and keeps only those without K/R in the middle
	 * @param tandemOutputFile
	 * 			a given tandem output file obtained by using "Peptide Hit Results Processor"
	 * @param numScan
	 * 			total number of scan
	 * @return the results for each scan
	 * @notes: the scan number in the file output by Tandem stars from 1
	 * */
	public Vector<Vector<String>> loadTandemResult(String tandemOutputFile, int numScan) throws IOException{
		FileInputStream fin = new FileInputStream(tandemOutputFile);
		BufferedReader in = new BufferedReader(new InputStreamReader(fin));
		
		Vector<Vector<String>> tandemResult = new Vector<Vector<String>>();
		for (int i=0; i<=numScan; i++){
			Vector<String> re = new Vector<String>();
			tandemResult.add(re);
		}
		
		String delims = "[\t ]+";
		String line = in.readLine();//header line
		line = in.readLine();
		while (line != null){
			String[] re = line.split(delims);
			int scanNum = Integer.parseInt(re[1]) - 1;//make it start from 0
			if (containSpecialCharacter(re[8])==false)//remove peptides with special characters (#, $ or @)
			{
				String pepString = this.extractPeptide(re[8]);
				if (this.containKRInMiddle(pepString)==false)
					tandemResult.get(scanNum).add(pepString);
			}
			line = in.readLine();
		}
		
		int count = 0;
		for (int i=0; i<=numScan; i++)
			if (tandemResult.get(i).size() != 0)
				count++;
		System.out.println("#scans in tandem results:"+count);
		
		in.close();
		fin.close();
		return tandemResult;
	}
	
	private boolean containSpecialCharacter(String peptide){
		//if (peptide.indexOf('#')==-1 || peptide.indexOf('@')==-1 || peptide.indexOf('$')==-1)
			//return false;
		if (peptide.contains("#"))
			return true;
		if (peptide.contains("$"))
			return true;
		if (peptide.contains("@"))
			return true;
		return false;
	}
	
	/**
	 * extracts useful peptide string from the given string
	 * @param peptide
	 * 			a peptide string in the tandem output, e.g. K.HGLEVIYMIEPIDEYC*VQQLK.E
	 * @return the part between the two '.' with special characters excluded
	 * */
	private String extractPeptide(String peptide){
		String str = peptide;
		//delete annotation characters from the peptide
		str = str.replace("*", "");
		str = str.replace("#", "");
		str = str.replace("@", "");
		str = str.replace("$", "");
		
		//delete the first and last two characters
		str = str.substring(2, str.length()-2);
		
		return str;
	}
	
	/**
	 * determines whether a peptide contains K/R in the middle (including the first position).
	 * @param peptide
	 * 			a given peptide string
	 * @return true if it contains K/R in the peptide except the last position
	 * Notes: only the last character can be K/R
	 * */
	/*private boolean containKRInMiddle(String peptide){
		boolean result = true;
		String tmp = peptide.substring(0, peptide.length()-1);
		if (tmp.indexOf('K')==-1 && tmp.indexOf('R')==-1)
			result = false;
		
		return result;
	}*/

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
	 * compares tandem results with our results
	 * @param tandemOutputFile
	 * 			a given tandem output file obtained by using "Peptide Hit Results Processor"
	 * @param myFile
	 * 			result file returned by my program
	 * @param matchFile
	 * 			a file name for storing matching scan numbers
	 * @param mismatchFile
	 * 			a file name for storing mismatching scan numbers
	 * */
	public void compareWithMyResult(String tandemOutputFile, String myFile, String matchFile, String mismatchFile, int numScan) throws IOException{
		
		Vector<Vector<String>> tandemResult = loadTandemResult(tandemOutputFile, numScan);
		Vector<Vector<String>> myResult  = my.loadMyResult(myFile, numScan);
		PrintWriter matchWriter = new PrintWriter(matchFile);
		PrintWriter mismatchWriter = new PrintWriter(mismatchFile);
		
		boolean match = false;
		int num_match = 0;
		int total = 0;
		for (int i=0; i<=numScan; i++){
			if (tandemResult.get(i).size() == 0)
				continue;
			if (myResult.get(i).size() == 0)
				continue;
			total++;
			match = false;
			for (int j=0; j<myResult.get(i).size(); j++){
				String tmp1 = myResult.get(i).get(j);
				for (int k=0; k<tandemResult.get(i).size(); k++){
					String tmp2 = tandemResult.get(i).get(k);
					if (tmp1.equals(tmp2)){
						//num_match++;
						match = true;
						//break;
					}
				}
			}
			if (match){
				num_match++;
				matchWriter.println(i);
			}
			else
				mismatchWriter.println(i);
		}
		
		double matchRate = num_match*1.0/total;
		System.out.println("total:"+total+"  match:"+num_match);
		System.out.println("match rate:"+ matchRate);
		
		matchWriter.close();
		mismatchWriter.close();
	}
	
	
	public void countScanWithoutPattern(String tandemOutputFile, String myFile, int numScan) throws IOException{
		MassTable mt = new MassTable();
		PeptideFileHandler pf = new PeptideFileHandler(mt);
		HashSet[] tagSet = pf.extractPatterns("deng_fei_data");
		Vector<Vector<String>> tandemResult = loadTandemResult(tandemOutputFile, numScan);
		Vector<Vector<String>> myResult = my.loadMyResult(myFile, numScan);
		
		int total = 0;
		for (int i=0; i<tagSet.length; i++)
			if (tandemResult.get(i).size() > 0)
				total++;
		System.out.println("total:"+total);
		
		int count = 0;
		for (int i=0; i<tagSet.length; i++)
			if (tandemResult.get(i).size()>0 && tagSet[i+1].size()==0)
				count++;
		System.out.println("#scans appeared in msgf search result but have no patterns generated in our program:"+count);
	}
	
	
	private void writeToFile(String tandemOutputFile, String outFile, int numScan) throws IOException{
		PrintWriter writer = new PrintWriter(outFile);
		Vector<Vector<String>> tandemResult = loadTandemResult(tandemOutputFile, numScan);
		for (int i=0; i<=numScan; i++)
			if (tandemResult.get(i).size() != 0)
				writer.println(i+"\t"+tandemResult.get(i));
		writer.close();
		
	}
	
	public static void main(String[] args) throws IOException{
		TandemComparator tc = new TandemComparator();
		tc.compareWithMyResult("HCDOutput.2014_01_26_16_25_02.t.txt", "output-1-27", "matchFile", "mismatchFile", 38000);
		//tc.countScanWithoutPattern("HCDOutput.2014_01_26_16_25_02.t.txt", "output-1-27", 38000);
	}
}

