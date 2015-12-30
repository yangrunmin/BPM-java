import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.Vector;


public class MyResultHandler {

	/**
	 * loads my results
	 * @param myFile
	 * 			my result file
	 * @param numScan
	 * 			total number of scans
	 * @return Vector<Vector<String>>, a scan may contains more than one result
	 * @notes: In the file ouput by my program, the scan number starts from 0
	 * */
	public Vector<Vector<String>> loadMyResult(String myFile, int numScan) throws IOException{
		FileInputStream my = new FileInputStream(myFile);
		BufferedReader myReader = new BufferedReader(new InputStreamReader(my));
		
		Vector<Vector<String>> myResult = new Vector<Vector<String>>();
		for (int i=0; i<=numScan; i++){
			Vector<String> re = new Vector<String>();
			myResult.add(re);
		}
		
		String delims = "[\t ]+";
		String line = myReader.readLine();
		while (line != null){
			String[] re = line.split(delims);
			int scanNum = Integer.parseInt(re[0]);
			myResult.get(scanNum).add(re[1]);
			line = myReader.readLine();
		}
		
		getResultStatistics(myResult);
		
		myReader.close();
		my.close();
		
		return myResult;
	}
	
	/**
	 * gets some of the statistics of the output of my program, e.g. total scans with results
	 * @param myResults
	 * 			my results
	 * */
	private void getResultStatistics(Vector<Vector<String>> myResult){
		//System.out.println("**********MyResult statistics**********\n");
		int size = myResult.size();
		int numScanWithResults = 0;
		for (int i=0; i<size; i++)
			if (myResult.get(i).size() > 0)
				numScanWithResults++;
		//System.out.println("Total number of scans with at least one search result:"+numScanWithResults);
		
		//System.out.println("**********MyResult statistics ends**********\n");
	}
	
	/**
	 * Test function, 
	 * @throws IOException 
	 * */
	public void compareTwoResults(String result1, String result2, int numScan) throws IOException{
		boolean[] m1 = loadMatchResult(result1, numScan);
		boolean[] m2 = loadMatchResult(result2, numScan);
		PrintWriter f1 = new PrintWriter("gap-only");
		PrintWriter f2 = new PrintWriter("ordinary-only");
		
		int count_BOTH = 0;
		for (int i=0; i<=numScan; i++){
			if (m1[i] && m2[i])
				count_BOTH ++;
			else if (m1[i] && !m2[i])
				f1.println(i);
			else if(!m1[i] && m2[i])
				f2.println(i);
		}
		System.out.println("in both:"+count_BOTH);
		
		f1.close();
		f2.close();
		
	}
	
	/**
	 * Test function, given a match file, load the match result (match or not for each scan)
	 * @param result
	 * 			the match result file, where each number a line
	 * @param numScan
	 * 			total number of scans
	 * @return boolean[numScan+1], boolean[i]=true if it appears in the file
	 * */
	private boolean[] loadMatchResult(String result, int numScan) throws IOException{
		FileInputStream my = new FileInputStream(result);
		BufferedReader myReader = new BufferedReader(new InputStreamReader(my));
		boolean[] matchResult = new boolean[numScan+1];
		for (int i=0; i<=numScan; i++)
			matchResult[i] = false;
		String line = myReader.readLine();
		while (line != null){
			int scanNum = Integer.parseInt(line);
			matchResult[scanNum] = true;
			line = myReader.readLine();
		}
		myReader.close();
		
		return matchResult;
	}
	
	
	public static void main (String[] args) throws IOException{
		MyResultHandler my = new MyResultHandler();
		//Vector<Vector<String>> myResult = my.loadMyResult("output-1-27", 38000);
		//System.out.println(myResult.get(37798));
		my.compareTwoResults("match-gap-8-2.0", "match-ordinary-8-2.0", 38000);
		
		
		System.out.println("hello world!");
	}
}
