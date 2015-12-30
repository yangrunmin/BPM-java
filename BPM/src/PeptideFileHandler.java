import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Vector;


public class PeptideFileHandler {
	private int numScan;
	private int[] massValue;
	private MassTable mt;
	//public static int fixedPatternLength = 8;
	public static boolean gapMode = true;	
	/**
	 * constructor
	 * */
	public PeptideFileHandler(MassTable massTable){
		this.mt = massTable;
	}
	
	/**
	 * returns the number of scans in the dataset
	 * @return the number of scans
	 * */
	public int getNumScan(){
		return numScan;
	}
	
	
	/**
	 * get the mass value (from the peptide dataset) associated with each scan
	 * @param scanNum
	 * 			the number of a scan
	 * */
	public int getMassValue(int scanNum){
		return massValue[scanNum];
	}
	
	/**
	 * Sets the parameter controlling the length of pattern length
	 * @param length
	 * 			a given integer number, default 3 
	 * */
	//public void setFixPatternLength(int length){
		//this.fixedPatternLength = length;
	//}
	
	
	
	/**
	 * gets the number of scans in the dataset(obtained from de novo sequencing)
	 * @param peptideFile
	 * 			file name for the dataset
	 * @return the max scanNum in the dataset
	 * */
	public void computeNumScan(String peptideFile) throws IOException{
		int max_ScanNum = 0;
		FileInputStream fis = new FileInputStream(peptideFile);
		BufferedReader in = new BufferedReader(new InputStreamReader(fis));
		
		String delims = "[, \t]+";
		String line = in.readLine();//the first line lists the attributes of each column
		line = in.readLine();
		while (line != null){
			String[] re = line.split(delims);
			int scanNum = Integer.parseInt(re[0]);
			if (max_ScanNum < scanNum)
				max_ScanNum = scanNum;			
			line = in.readLine();
		}
		
		in.close();
		fis.close();
		numScan = max_ScanNum;
	}
	
	/**
	 * preprocess a peptide to resolve unexpected amino acids
	 * @param pepStr
	 * 			the given peptide string
	 * @return the processed peptide string
	 * */
	private String preprocessPeptide(String pepStr){
		pepStr = pepStr.replace("(+57.02)", "");
		pepStr = pepStr.replace('B', 'A');
		pepStr = pepStr.replace('J', 'A');
		pepStr = pepStr.replace('O', 'A');
		pepStr = pepStr.replace('U', 'A');
		pepStr = pepStr.replace('X', 'A');
		pepStr = pepStr.replace('Z', 'A');
		
		return pepStr;
	}
	
	/**
	 * constructs patterns for each peptide in the dataset, one entry for each scan, 
	 * multiple patterns are counted only once for a scan
	 * @param peptideFile
	 * 			the file name for the dataset
	 * @param numScan
	 * 			total number of scans in the dataset
	 * @return the patterns extracted from the peptides in each scan
	 * */
	public HashSet[] extractPatterns(String peptideFile) throws IOException{
		FileInputStream fis = new FileInputStream(peptideFile);
		BufferedReader in = new BufferedReader(new InputStreamReader(fis));
		
		computeNumScan(peptideFile);
		massValue = new int[numScan+1];
		HashSet[] tagSet = new HashSet[numScan+1];
		for (int i=0; i<=numScan; i++)
			tagSet[i] = new HashSet<Pattern>();
		
		String delims = "[, \t]+";
		String line = in.readLine();//the first line contains attributes for each column
		line = in.readLine();
		while (line != null){
			line = this.preprocessPeptide(line);
			String[] re = line.split(delims);
			int scanNum = Integer.parseInt(re[0]);
			double tmp = Double.parseDouble(re[7]);
			massValue[scanNum] = (int)Math.round(tmp*mt.getScaleFactor());
			int[] confidence = new int[re[1].length()];
			for (int i=0; i<re[1].length(); i++)
				confidence[i] = Integer.parseInt(re[i+9]);
			
			Peptide pep = new Peptide(scanNum, re[1], confidence, mt);
			pep.constructPattern(tagSet, gapMode);
		
			line = in.readLine();
			
			
		}
		
		if (gapMode) 
			printTagSet(tagSet, "gapTag");
		else
			printTagSet(tagSet, "ordinaryTag");
		
		in.close();
		fis.close();
	
		return tagSet;
	}
	
	/**
	 * split the obtained patterns into subpatterns with fixed length indicated by <fixedPatternLength>
	 * @param tagSet
	 * 			the obtained list of patterns
	 * */
	/*public HashSet[] splitPatterns2FixedLength(HashSet[] tagSet){
		HashSet[] newTagSet = new HashSet[numScan+1];
		for (int i=0; i<=numScan; i++)
			newTagSet[i] = new HashSet<Pattern>();
		
		for (int i=0; i<=numScan; i++){
			Iterator ite = tagSet[i].iterator();
			while (ite.hasNext()){
				Pattern pattern = (Pattern)ite.next();
				Vector<Pattern> subPatternList = pattern.splitPattern(fixedPatternLength);
				//tagSet[i].remove(pattern);
				//pattern = null;
				for (int j=0; j<subPatternList.size(); j++)
					newTagSet[i].add(subPatternList.get(j));
			}
		}
		
		return newTagSet;
	}*/
	
	/**
	 * Test function. Get the average length of tags
	 * */
	public void getTagSetDetails(HashSet[] tagSet){
		System.out.println("======gapped tag details==========");
		int patternLength = Peptide.PATTERN_LENGTH;
		int[] len = new int[patternLength+1];
		for (int i=0; i<=patternLength; i++)
			len[i] = 0;
		
		int numTag = 0;
		for (int i=0; i<=numScan; i++){
			Iterator ite = tagSet[i].iterator();
			while (ite.hasNext()){
				Pattern pattern = (Pattern) ite.next();
				numTag++;
				len[pattern.length()]++;
			}
		}
		int sum = 0;
		for (int i=1; i<=patternLength; i++){
			sum += i*len[i];
			System.out.println("#tags with length-"+i+": "+len[i]);
		}
		
		double ave = sum*1.0/numTag;
		System.out.println("#tags:"+numTag);
		System.out.println("avarage length of tags:"+ave);
		
		
		System.out.println("============end here=================");
	}
	
	/**
	 * generate a sub dataset by selecting the top peptides (reported by PEAKS)
	 * @param peptideFile
	 * 			the original data file
	 * @param outFile
	 * 			output data file
	 * @param numData
	 * 			number of peptides selected for each spectrum
	 * */
	public void getPartialData(String peptideFile, String outFile, int numData) throws IOException{
		FileInputStream fis = new FileInputStream(peptideFile);
		BufferedReader in = new BufferedReader(new InputStreamReader(fis));
		PrintWriter out = new PrintWriter(outFile);
		
		computeNumScan(peptideFile);
		int[] count = new int[numScan+1];
		for (int i=0; i<=numScan; i++)
			count[i] = 0;
		System.out.println(numScan);
		
		String delims = "[, \t]+";
		String line = in.readLine();
		out.println(line);
		line = in.readLine();
		while (line != null){
			String[] re = line.split(delims);
			int scanNum = Integer.parseInt(re[0]);
			if (count[scanNum] < numData){
				out.println(line);
				count[scanNum]++;
			}
			line = in.readLine();
		}
		count = null;
		out.close();
		in.close();
		fis.close();
	}
	
	private void freeTagList(HashSet[] tagSet){
		for (int i=0; i<=numScan; i++)
			tagSet[i].clear();
		
		tagSet = null;
	}
	
	/**
	 * Test function, print tags(patterns) extracted from the peptides in each scan
	 * @param tagSet
	 * 			the tags(patterns) extracted from the peptides in each scan
	 * @param outFile
	 * 			the disk file for holding the output tags
	 * */
	public void printTagSet(HashSet[] tagSet, String outFile) throws IOException{
		PrintWriter writer = new PrintWriter(outFile);
		for (int i=0; i<=numScan; i++){
			if (tagSet[i].size() == 0)
				continue;
			writer.print("Scan "+i+": ");
			Iterator ite = tagSet[i].iterator();
			while (ite.hasNext()){
				Pattern pattern = (Pattern) ite.next();
				writer.print(pattern+", ");
			}
			writer.println();
		}
		
		writer.close();
	}
	
	public static void main(String[] args) throws IOException{
		MassTable mt = new MassTable();
		PeptideFileHandler.gapMode = true;
		PeptideFileHandler pf = new PeptideFileHandler(mt);
		HashSet[] tagSet = pf.extractPatterns("deng_fei_data");
		//pf.printTagSet(tagSet, "tagSet");
		
		System.out.println("hello world!");
	}
}
