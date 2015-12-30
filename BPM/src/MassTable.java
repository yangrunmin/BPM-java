import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.HashSet;
import java.util.Iterator;



public class MassTable {
	
	private int scaleFactor;//比例系数
	private double maxSize_original;//max mass size before scaling
	private int maxSize;//max mass size after scaling
	private HashSet[] massTable;
	private DecimalFormat df = new DecimalFormat("##0.00");
	
	//all the 19 amino acids, I and L have the same mass, omitted in the table
	public static final char[] aminoAcid = {'A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'};
	
	public static final int[] aminoAcidIndex = {0, -1, 4, 3, 5, 12, 7, 8, 9, -1, 10, 9, 11, 2, -1, 13, 6, 1, 14, 15, -1, 18, 16, -1, 17, -1};
											 //{A, B,  C, D, E, F,  G, H, I, J,  K,  L,  M,  N,  O, P,  Q, R,  S,  T,  U,  V,  W,  X,  Y, Z};

	//the mono mass of the 19 amino acids
	public static final double[] acidMass = {71.03711, 156.10111, 114.04293, 115.02694, 
		160.030654, 129.04259, 128.05858, 57.02146, 137.05891, 113.08406, 128.09496, 
		131.04049, 147.06841, 97.05276, 87.03203, 101.04768, 186.07931, 163.06333, 99.06841};
	public static int[] scaledMass;//序号是0-18
	public static int[] scaledMass2;//序号是ASSIC码
	
	/**
	 * constructor
	 * */
	public MassTable(){
		scaleFactor = 100;
		maxSize_original = 500.0;
		init();
	}
	
	/**
	 * constructor
	 * @param scale
	 * 			scale factor
	 * @param max_size
	 * 			maximum size of the allowed mass (before scaling) 
	 * */
	public MassTable(int scale, double max_size){
		scaleFactor = scale;
		maxSize_original = max_size;
		init();
	}
	
	/**
	 * initialization
	 * */
	public void init(){
		maxSize = (int)Math.round(maxSize_original*scaleFactor);
		scaledMass = new int[19];
		for (int i=0; i<19; i++)
			scaledMass[i] = (int) Math.round(acidMass[i]*scaleFactor);
		//allocate space for massTable
		massTable = new HashSet[maxSize+1];
		for (int i=0; i<=maxSize; i++)
			massTable[i] = new HashSet<String>();
		init2();
	}
	
	private void init2(){
		scaledMass2 = new int[255];
		for (int i=0; i<19; i++)
			scaledMass2[aminoAcid[i]] = this.getScaledMass(aminoAcid[i]);
		scaledMass2['L'] = scaledMass2['I'];
	}
	
	/**
	 * get the "scaleFactor"
	 * @return scaleFactor 
	 * */
	public int getScaleFactor(){
		return scaleFactor;
	}
	
	/**
	 * get the maximum mass (with scaling) allowed
	 * @return maxSize
	 * */
	public int getMaxSize(){
		return maxSize;
	}

	/**
	 * get the scaled mass
	 * @param index
	 * 			the index of the mass
	 * @return
	 * 			the scaledMass at the given pos. 
	 */
	public int getScaledMass(int index){
		return scaledMass[index];
	}
	
	/**
	 * get the scaled mass
	 * @param ch
	 * 			the given amino acid
	 * @return
	 * 			the scaled mass value for the given amino acid
	 * */
	public int getScaledMass(char ch){
		int offset = Character.toUpperCase(ch)-'A';
		//int offset = ch - 'A';
		if (offset<0 || offset>25){
			System.out.println("Unexpected amino acid "+ch);
			return -1;
		}
		int index = aminoAcidIndex[offset];
		if (index>=0 && index<19)
			return scaledMass[index];
		else {
			System.out.println("Unexpected amino acid "+ch);
			return -1;
		}
	}
	
	public int getScaledMass2(char ch){
		return scaledMass2[ch];
	}
	
	/**
	 * get all the mass strings associated with a given mass value
	 * @param mass
	 * 			the given mass value
	 * @return an iterator referring to the mass strings
	 * */
	public Iterator getMassString(int mass){
		Iterator ite = massTable[mass].iterator();
		return ite;
	}
	
	
	/**
	 * get the number of mass strings associated with a given mass value
	 * @param mass
	 * 			the given mass
	 * @return number of mass strings associated with the given mass 
	 */
	public int getNumMassString(int mass){
		return massTable[mass].size();
	}
	
	/**
	 * given a scaled mass, find the corresponding amino acid
	 * @param scaled_mass
	 * 			the given mass value
	 * @return the corresponding amino acid if found; '*' otherwise.
	 */
	public char getAminoAcid(int scaled_mass){
		for (int i=0; i<19; i++)
			if (scaled_mass == scaledMass[i])
				return aminoAcid[i];
		return '*';
	}
	
	/**
	 * generate the mass table, find the mass strings associated with each mass value
	 * */
	public void generateMassTable(){
		for (int i=0; i<=maxSize; i++)
			for (int j=0; j<19; j++)
				extend(i, scaledMass[j]);
	}
	
	/**
	 * find all the strings s such that mass(s)+singleMass=targetMass
	 * @param targetMass
	 * 			target mass value
	 * @param singleMass 
	 * 			a single mass value  
	 * */
	private void extend(int targetMass, int singleMass){
		int sourceMass = targetMass - singleMass;
		if (sourceMass < 0) return;
		else if (sourceMass == 0){
			StringBuilder newElement = new StringBuilder();
			newElement.append(getAminoAcid(singleMass));
			massTable[targetMass].add(newElement.toString());
			return;
		}
		
		int size = massTable[sourceMass].size();
		if (size == 0)
			return;
		String currentElement;
		StringBuilder newElement;
		Iterator ite = massTable[sourceMass].iterator();
		while (ite.hasNext()){
			currentElement = (String)ite.next();
			int elementSize = currentElement.length();
			for (int j=0; j<=elementSize; j++){
				newElement = new StringBuilder(currentElement);
				newElement.insert(j, getAminoAcid(singleMass));
				massTable[targetMass].add(newElement.toString());
			}
		}
	}
	
	/**
	 * get the max length of a mass string over the entire table
	 * @return max length of a mass string over the entire table
	 * */
	public int getMaxPatternNumberLen(){
		int maxPatternNumberLen = 0;
		for (int i=0; i<=maxSize; i++){
			Iterator ite = massTable[i].iterator();
			while (ite.hasNext()){
				String currentElement = (String) ite.next();
				if (maxPatternNumberLen < currentElement.length()){
					maxPatternNumberLen = currentElement.length();
				}
			}
		}
		
		return maxPatternNumberLen;
	}
	
	/**
	 * compute the R value, i.e. (n_ij)^{-1/j} for each row
	 * @param filename
	 * 			the output file for storing the R value for each row
	 * @param maxPatternNumber
	 * 			max length of a mass string  
	 * @return max R value over the entire table
	 * */
	public double computeRValue(String filename, int maxPatternNumberLen) throws FileNotFoundException{
		PrintWriter writer = new PrintWriter(filename);
		double maxRValue = 0.0, currentR;
		int maxIndex = 0, maxLen = 1, maxCount = 0;
		int[] counter = new int[maxPatternNumberLen+1];
		
		for (int i=0; i<=maxSize; i++){
			if (massTable[i].size() == 0)
				continue;
			for (int j=0; j<=maxPatternNumberLen; j++)
				counter[j] = 0;
			Iterator ite = massTable[i].iterator();
			while (ite.hasNext()){
				String currentElement = (String) ite.next();
				counter[currentElement.length()]++;
			}
			
			//write to file
			writer.print(i+":");
			for (int j=1; j<=maxPatternNumberLen; j++){
				if (counter[j]>0){
					currentR = Math.pow(counter[j], 1.0/j);
					if (maxRValue < currentR){
						maxRValue = currentR;
						maxLen = j;
						maxCount = counter[j];
						maxIndex = i;
					}
					writer.print("("+j+","+counter[j]+","+df.format(currentR)+") ");
				}
			}
			writer.println();
		}
		
		System.out.println("maxR value achieves at Entry "+maxIndex+" with "+maxCount+" length-"+maxLen+" pattern numbers");
		
		writer.close();
		return maxRValue;
	}

	/**
	 * write the mass table to disk 
	 * @param filename 
	 * 			file name for storing the table 
	 */
	public void printMassTable2File(String filename) throws FileNotFoundException{
		PrintWriter writer = new PrintWriter(filename);
		for (int i=0; i<=maxSize; i++){
			if (massTable[i].size()==0)
				continue;
			writer.print(i+":");
			Iterator ite = massTable[i].iterator();
			boolean first = true;
			while (ite.hasNext()){
				if (!first) writer.print(", ");
				first = false;
				writer.print(ite.next());
			}
			
			writer.println();
			writer.println();
		}
		writer.close();
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
			mass += getScaledMass2(pepString.charAt(i));
		mass += 18*getScaleFactor();//add the mass of H2O
		
		
		return mass;
	}
	
	
	public static void main(String[] args) throws IOException{
		long startTime = System.currentTimeMillis();
		MassTable mt = new MassTable();
		
		mt.generateMassTable();
		long endTime = System.currentTimeMillis();
		System.out.println("Time:"+(endTime-startTime)/1000);
		mt.printMassTable2File("massTable.txt");
		
		int maxlen=mt.getMaxPatternNumberLen();
		System.out.println("max length of a patter string:"+maxlen);
		System.out.println("maximum R value:"+mt.computeRValue("r.txt", maxlen));
		long endTime2 = System.currentTimeMillis();
		System.out.println("Time:"+(endTime2-endTime)/1000);
		
		System.out.println("hello world!");
	}
}
