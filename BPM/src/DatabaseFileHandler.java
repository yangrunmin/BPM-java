import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.Random;
import java.util.Vector;


public class DatabaseFileHandler {
	
	/**
	 * load protein database from file (fasta format, each protein sequence begins with a '>')
	 * @param proteinDatabaseFile 
	 * 			file name for the protein database
	 * @return an instance of ProteinDatabase
	 * */
	public ProteinDatabase loadDatabase(String proteinDatabaseFile) throws IOException{
		ProteinDatabase database = new ProteinDatabase();
		FileInputStream fis = new FileInputStream(proteinDatabaseFile);
		BufferedReader in = new BufferedReader(new InputStreamReader(fis));
		
		boolean isFirstSeq = true;
		boolean isFirstAminoAcid = true;
		StringBuilder sequence = new StringBuilder();
		StringBuilder individualSeq = new StringBuilder();
		String line = in.readLine();
		while (line != null){
			if (line.startsWith(">")){
				database.addProteinID(line);
				if (!isFirstSeq){
					sequence.append('#');
					database.addIndividualSeq(individualSeq.toString());
					individualSeq.setLength(0);
					isFirstAminoAcid = true;
				}
				isFirstSeq = false;
			}
			else{
				if (isFirstAminoAcid){
					line = removeFirstAminoAcid(line);
					isFirstAminoAcid = false;
				}
				line = handleUndefinedCharacter(line);
				individualSeq.append(line);
				sequence.append(line);
			}
			line = in.readLine();
		}
		sequence.append('$');
		String text = sequence.toString();
		text = text.replace('L', 'I');
		database.setSequence(text);
		database.addIndividualSeq(individualSeq.toString());
		
		sequence = null;
		individualSeq = null;
		in.close();
		fis.close();
		return database;
	}
	
	/**
	 * replace unexpected characters other than the 20 amino acids with 'A'
	 * @param text
	 * 			a given amino acid sequence
	 * @return the resulting sequence with all unexpected characters replaced by 'A'
	 * */
	private String handleUndefinedCharacter(String text){
		//replace 'L' with 'I'; they have equal mass
		//text = text.replace('L', 'I');
		//replace unusual symbols with A
		text = text.replace('B', 'A');
		text = text.replace('J', 'A');
		text = text.replace('O', 'A');
		text = text.replace('U', 'A');
		text = text.replace('X', 'A');
		text = text.replace('Z', 'A');
		
		return text;
	}
	
	/**
	 * deletes the first character of a protein if it is 'M'
	 * @param text
	 * 			a text (sequence)
	 * @return the modified sequence
	 * */
	private String removeFirstAminoAcid(String text){
		if (text.startsWith("M"))
			text = text.substring(1);
		
		return text;
	}
	
	/**
	 * TEST function, modify the line beginning with ">" by using the index of each sequence (start from 0)
	 * @param filename input filename for the proteome database
	 * @param filename2 the output file for holding the new format
	 * */
	public void indexDatabase(String proteinDatabaseFile, String outputFile) throws IOException{
		FileInputStream fis = new FileInputStream(proteinDatabaseFile);
		BufferedReader in = new BufferedReader(new InputStreamReader(fis));
		FileOutputStream fos = new FileOutputStream(outputFile);
		BufferedWriter out = new BufferedWriter(new OutputStreamWriter(fos));
		
		int count = 0;
		String line = in.readLine();
		while (line != null){
			if (line.startsWith(">")){
				out.write(">"+count);
				count++;
				out.newLine();
			}
			else {
				out.write(line);
				out.newLine();
			}
			line = in.readLine();
		}
	
		in.close();
		fis.close();
		out.close();
		fos.close();
	}
	
	public static void main(String[] args) throws IOException{
		DatabaseFileHandler df = new DatabaseFileHandler();
		ProteinDatabase pd = df.loadDatabase("test2.fasta");
		Vector<Vector<Integer>> peptideTable = pd.extractPeptideFromSeq();
		pd.writePeptideTable2File(peptideTable, "peptideTable");
		
		System.out.println(pd.getSequence());
		for (int i=0; i<pd.size(); i++)
			System.out.println(pd.getSeqLength(i)+"  "+pd.getIndividualSeq(i));
		
		System.out.println("hello world!");
	}

}
