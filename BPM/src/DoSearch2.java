
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Vector;


public class DoSearch2 {
	
	private MassTable mt;
	private SuffixTree st;
	private PeptideFileHandler pf;
	private PatternSearch patternSearch;
	private static int numResult = 1000000;
	private ProteinDatabase proteinDatabase;
	private static double massTolerance = 0.05;
	public static int scaledMassTolerance;
	private static int MASS_WATER;
	public static boolean sort = false;
	public static long startTime;
	
	
	public DoSearch2() throws IOException{
		mt = new MassTable();
		mt.generateMassTable();
		pf = new PeptideFileHandler(this.mt);
		MASS_WATER = 18*mt.getScaleFactor();//the mass of H2O
	}
	
	
	/**
	 * compute the mass of a peptide
	 * @param pepString
	 * 			the given peptide
	 * @return the mass of the peptide
	 * */
	public int computePeptideMass(String pepString){
		int mass = 0;
		for (int i=0; i<pepString.length(); i++)
			mass += mt.getScaledMass2(pepString.charAt(i));
		//mass += 18*mt.getScaleFactor();//add the mass of H2O
		mass += MASS_WATER;
		
		return mass;
	}
	
	public static int getScaledMassTolerance(){
		//return (int)Math.round(massTolerance*mt.getScaleFactor());
		return scaledMassTolerance;
	}
	
	public void compScaledMassTolerance(double error){
		scaledMassTolerance = (int)Math.round(massTolerance*mt.getScaleFactor());
	}
	
	
	/**
	 * Test if a given peptide is within the mass tolerance
	 * @param pepString
	 * 			the given peptide string
	 * @param refMass
	 * 			the reference mass value
	 * @return true if yes; false otherwise
	 * */
	public boolean inErrorTolerance(LeafEdge leaf, int refMass){
		int mass = leaf.getPeptideMass();
		if (Math.abs(mass-refMass) <= getScaledMassTolerance())
			return true;
		return false;
	}
	
	
	public void search(String proteinDatabaseFile, String peptideFile, String outputFile) throws IOException{
		System.out.print("Loading database");
		DatabaseFileHandler da = new DatabaseFileHandler();
		proteinDatabase = da.loadDatabase(proteinDatabaseFile);
		computeRunningtime();
		
		System.out.print("Extract peptides from sequence");
		Vector<Vector<Integer>> peptideTable = proteinDatabase.extractPeptideFromSeq2();
		computeRunningtime();
		
		System.out.print("constructing suffix tree");
		this.st = new SuffixTree(proteinDatabase.getSequence(), proteinDatabase);
		computeRunningtime();
		
		System.out.print("Preprocessing");
		Preprocessing pp = new Preprocessing(st, proteinDatabase, this, peptideTable);
		pp.computeMassRangeForNodes(st.getRoot());
		proteinDatabase.freeMem(peptideTable);
		computeRunningtime();
	
		System.out.print("extract patterns from peptides");
		pf.computeNumScan(peptideFile);
		patternSearch = new PatternSearch(this.st, this.mt, this.pf);
		HashSet[] tagSet = pf.extractPatterns(peptideFile);
		computeRunningtime();
		
		System.out.print("searching...");
		search(tagSet, peptideTable, outputFile);
		computeRunningtime();
	}
	
	/**
	 * search peptides against the protein database
	 * @param tagSet
	 * 			extracted tags (patterns) from the peptide file
	 * @param peptideTable
	 * 			extracted peptide table from the protein database
	 * @param outFile
	 * 			output file name
	 * */
	private void search(HashSet[] tagSet, Vector<Vector<Integer>> peptideTable, String outFile) throws IOException{
		PrintWriter writer = new PrintWriter(outFile);
		for (int i=0; i<=pf.getNumScan(); i++){
			patternSearch.clearNodePool();
			HashSet<LeafEdge> searchResult = new HashSet<LeafEdge>();
			Iterator ite = tagSet[i].iterator();
			while (ite.hasNext()){
				Pattern pattern = (Pattern) ite.next();
				ArrayList<LeafEdge> leafSet = patternSearch.searchPattern(pattern, i);
				if (leafSet != null){
					for (int j=0; j<leafSet.size(); j++)
						add2ResultSet(peptideTable, searchResult, leafSet.get(j), i);
				}
			}
			if (searchResult.size() > 0){
				writeResult2File(writer, searchResult, i);
			}
		}
		writer.close();
	}
	
	/**
	 * adds a search result to the result set
	 * @param peptideTable 
	 * 			peptide table extracted from protein sequence
	 * @param searchResult<SuffixPosition, Integer>, where SuffixPosition refers to the result and Integer refers to the time the result occurred 
	 * 			the result set
	 * @param sp
	 * 			the result to be added
	 * */
	public void add2ResultSet(Vector<Vector<Integer>> peptideTable, HashSet<LeafEdge> searchResult, LeafEdge leaf, int scanNum){
		if (searchResult.contains(leaf))
			return;
		searchResult.add(leaf);
		//if (inErrorTolerance(leaf, pf.getMassValue(scanNum)))
			//searchResult.add(leaf);
	}
	
	
	/**
	 * writes the top search results (indicated by the parameter numResult) to file
	 * @param writer
	 * 			file pointer referring to the disk file for holding the output
	 * @param searchResult
	 * 			the obtained search result
	 * */
	public void writeResult2File(PrintWriter writer, HashSet<LeafEdge> searchResult, int scanNum){
		Iterator ite = searchResult.iterator();
		while (ite.hasNext()){
			LeafEdge leaf = (LeafEdge)ite.next();
			//SuffixPosition sp = (SuffixPosition)ite.next();
			writer.print((scanNum-1)+"\t");//scanNum starts from 1, -1 to make it start from 0 in the output file
			writer.print(leaf.getPepString()+"\t");
			writer.print(proteinDatabase.getProteinID(leaf.getSeqNum())+"\t");
			writer.print("<"+leaf.getSeqNum()+","+leaf.getLeafIndex()+">");
			writer.println();
		}
		
		searchResult.clear();
		searchResult = null;
	}
	
	public void setParameters(int minPatternLength, int fixedPatternLength, double searchError, double resultError, int minConfidence){
		Peptide.patternLengthCutoff = minPatternLength;
		Peptide.PATTERN_LENGTH = fixedPatternLength;
		PatternSearch.errorTolerance = searchError;
		DoSearch2.massTolerance = resultError;
		this.compScaledMassTolerance(massTolerance);
		Peptide.MIN_CONFIDENCE = minConfidence;
	}
	
	public void computeRunningtime(){
		
		long end = System.currentTimeMillis();
		long duration = (end-startTime)/1000;
		System.out.println(" (Time elapsed: "+duration+")");
		startTime = end;
	}

	/**
	 * args[0]: protein database
	 * args[1]: input data file
	 * args[2]: minimum pattern length
	 * args[3]: fixed pattern length
	 * args[4]: search error tolerance
	 * args[5]: search result error tolerance
	 * args[6]: confidence cutoff
	 * args[7]: output file name
	 */
	
	public static void main(String[] args) throws IOException{
		startTime = System.currentTimeMillis();
		long start = startTime;
		
		DoSearch2 doSearch = new DoSearch2();
		//set the parameters
		PeptideFileHandler.gapMode = true;
		int minPatternLength = Integer.parseInt(args[2]);
		int fixedPatternLength = Integer.parseInt(args[3]);
		double searchError = Double.parseDouble(args[4]);
		double resultError = Double.parseDouble(args[5]);
		int minConfidence = Integer.parseInt(args[6]);
		doSearch.setParameters(minPatternLength, fixedPatternLength, searchError, resultError, minConfidence);
		System.out.println("List of search parameters:");
		System.out.println("minPatternLength:"+minPatternLength);
		System.out.println("fixedPatternLength:"+fixedPatternLength);
		System.out.println("searchError:"+searchError);
		System.out.println("resultError:"+resultError);
		System.out.println("confidenceCutOff:"+minConfidence);
		System.out.println("gapped-tag:"+PeptideFileHandler.gapMode);
		System.out.println();
		
		doSearch.search(args[0], args[1], args[7]);
		long endTime = System.currentTimeMillis();
		long duration = (endTime - start)/1000;
		System.out.println("Total Running time:"+duration+" seconds");
		System.out.println("hello world!");
	}
	
	/*
	public static void main(String[] args) throws IOException{
		//long startTime = System.currentTimeMillis();
		startTime = System.currentTimeMillis();
		long start = startTime;
		//set the parameters
		PeptideFileHandler.gapMode = true;
		int minPatternLength = 8;
		int fixedPatternLength = 8;
		double searchError = 0.05;
		double resultError = 0.05;
		int minConfidence = 60;
		DoSearch2 doSearch = new DoSearch2();
		doSearch.setParameters(minPatternLength, fixedPatternLength, searchError, resultError, minConfidence);
		
		
		doSearch.search("test.fasta", "deng_fei_data", "output-2-6");
		long endTime = System.currentTimeMillis();
		long duration = (endTime - start)/1000;
		System.out.println("Total Running time:"+duration+" seconds");
		System.out.println("hello world!");
		
		System.out.println("\n=============================");
		System.out.println("List of search parameters:");
		System.out.println("minPatternLength:"+Peptide.patternLengthCutoff);
		System.out.println("fixedPatternLength:"+Peptide.PATTERN_LENGTH);
		System.out.println("searchError:"+PatternSearch.errorTolerance);
		System.out.println("resultError:"+DoSearch2.massTolerance);
		System.out.println("confidenceCutOff:"+Peptide.MIN_CONFIDENCE);
		System.out.println("Use gap mode:"+PeptideFileHandler.gapMode);
		System.out.println();
	}*/
}
