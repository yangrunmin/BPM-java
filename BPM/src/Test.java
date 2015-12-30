import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Vector;


public class Test {
	
	public void generateACommand(int fixedPatternLength,double searchError, double resultError, int confidence, PrintWriter writer) throws IOException{
		int minPatternLength = fixedPatternLength;
		writer.print("java DoSearch2 test.fasta deng_fei_data "+minPatternLength+" "+fixedPatternLength+" ");
		writer.print(searchError+" "+resultError+" "+confidence+" ");
		writer.print("result/output-3-1"+"-p"+fixedPatternLength+"-c"+confidence+"-s"+searchError+"-r"+resultError+" >> ");
		writer.println("result/run-3-1"+"-p"+fixedPatternLength+"-c"+confidence+"-s"+searchError+"-r"+resultError+".log");
		writer.print("java MSGFComparator MSGF-1-29-Q.tsv ");
		writer.print("result/output-3-1"+"-p"+fixedPatternLength+"-c"+confidence+"-s"+searchError+"-r"+resultError+" >> ");
		writer.println("result/run-3-1"+"-p"+fixedPatternLength+"-c"+confidence+"-s"+searchError+"-r"+resultError+".log");
	}
	
	public void generateCommand(String filename) throws IOException{
		PrintWriter writer = new PrintWriter(filename);

		for (int k=4; k<=10; k++){
			for (int i=30; i<=50; i=i+10){
				generateACommand(k, 0, 2.0, i, writer );
				generateACommand(k, 0.05, 2.0, i, writer);
			}
			for (int i=70; i<=90; i=i+10){
				generateACommand(k, 0, 2.0, i, writer );
				generateACommand(k, 0.05, 2.0, i, writer);
			}

		}
		
		
		writer.close();
		
	}
	
	public void t(String peptideFile) throws IOException{
		MassTable mt = new MassTable();
		Peptide.PATTERN_LENGTH = 3;
		PeptideFileHandler pf = new PeptideFileHandler(mt);
		HashSet[] tagSet = pf.extractPatterns(peptideFile);
		pf.getTagSetDetails(tagSet);
		
		//pf.getPartialData(peptideFile, "deng_fei_data_top1", 1);
	
		System.out.println("Hello world");
		
	}
	
	public static void main(String[] args) throws IOException{
		/*MassTable mt = new MassTable();
		mt.generateMassTable();
		DatabaseFileHandler df = new DatabaseFileHandler();
		ProteinDatabase pd = df.loadDatabase("test.fasta");
		SuffixTree st = new SuffixTree(pd.getSequence(), pd);
		PatternSearch ps = new PatternSearch(st, mt);
		
		Pattern pattern3 = new Pattern();
		pattern3.add(11308);
		pattern3.add(11308);
		pattern3.add(7104);
		System.out.println(ps.searchPattern(pattern3));
		
		Pattern pattern = new Pattern();
		pattern.add(7104);
		pattern.add(7104);
		
		ArrayList<SuffixPosition> leafSet = ps.searchPattern(pattern);
		System.out.println(leafSet);
		
		Pattern pattern2 = new Pattern();
		pattern2.add(7104);
		pattern2.add(7104);
		pattern2.add(7104);
		System.out.println(ps.searchPattern(pattern2));*/	
		Test tt = new Test();
		tt.generateCommand("run.sh");
		//tt.t("deng_fei_data");

		
	}

}

class aaa{
	private boolean visited;
	private int num;
	public aaa(boolean hello, int num){
		visited = hello;
		this.num = num;
	}
	public void set(){
		visited = false;
	}
	public String toString(){
		return visited+" "+num;
	}
}
