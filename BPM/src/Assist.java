import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;


public class Assist {

	/**
	 * Test function, given a amino acid sequence, compute the sum of the mass of each amino acid
	 * @param str
	 * 			the given amino acid sequence
	 * @return the sum of the masses
	 * */
	public int computeMass(String str){
		MassTable mt = new MassTable();
		int mass = 0;
		for (int i=0; i<str.length(); i++)
			mass+= mt.getScaledMass(str.charAt(i));
		//System.out.println(mass);
		
		return mass;
	}
	
	/**
	 * Test function, find all the occurrences of a string in the protein database
	 * @param databaseFilename
	 * 			protein database file name
	 * @param queryString 
	 * 			the query string
	 * */
	public void findAll(String databaseFilename, String queryString) throws IOException{
		DatabaseFileHandler df = new DatabaseFileHandler();
		ProteinDatabase pd = df.loadDatabase(databaseFilename);
		pd.findAllOccurrence(queryString);
	}
	
	public static void main(String[] args) throws IOException{
		Assist ass = new Assist();
		System.out.println("mass sum:"+ass.computeMass("NSE"));
		DoSearch2 ds = new DoSearch2();
		System.out.println("Peptide mass:"+ds.computePeptideMass("TVKEEAEKPEREAK"));
		//ass.findAll("test.fasta", "VLPSEQESTK");
		System.out.println("hello world!");
	}
}
