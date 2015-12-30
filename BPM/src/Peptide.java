import java.util.HashSet;


public class Peptide {
	private int scanNum = 0;
	private String peptide;
	private int[] confidence;
	private MassTable mt;
	public static int MIN_CONFIDENCE = 60;
	public static int patternLengthCutoff = 8;
	public static int PATTERN_LENGTH = 8;
	
	/**
	 * constructor
	 * */
	public Peptide(){
		peptide = new String();
		mt = null;
		confidence = null;
	}
	
	/**
	 * constructor
	 * */
	public Peptide(int scanNum, String peptide, int[] confidenceValue, MassTable mt){
		this.scanNum = scanNum;
		this.mt = mt;
		this.peptide = peptide;
		int length = peptide.length();
		confidence = new int[length];
		
		for (int i=0; i<length; i++)
			this.confidence[i] = confidenceValue[i];
	}
	
	/**
	 * gets the peptide sequence
	 * @return peptide sequence
	 * */
	public String getPeptide(){
		return peptide;
	}
	
	/**
	 * gets the scanNum to which the peptide belongs
	 * @return the scanNum to which the peptide belongs
	 * */
	public int getScanNum(){
		return scanNum;
	}
	
	/**
	 * constructs a pattern from the peptide
	 * @param tagSet
	 * 			to where the generated tags are stored
	 * @param gapMode
	 * 			true if use gapped-tag, false otherwise
	 * @note all the amino acids are used, 
	 * */
	public void constructPattern2(HashSet[] tagSet, boolean gapMode){
		for (int i=0; i<=peptide.length()-PATTERN_LENGTH; i++){
			if (gapMode)
				constructGapTag(i, tagSet);
			else
				constructOridinaryTag(i, tagSet);
		}
	}
	
	private void constructGapTag(int startPos, HashSet[] tagSet){
		int mass_sum = 0;
		Pattern pattern = new Pattern();
		for (int i=startPos; i<startPos+PATTERN_LENGTH; i++){
			int currentMass = mt.getScaledMass(peptide.charAt(i));
			if (confidence[i] < MIN_CONFIDENCE){
				if (mass_sum+currentMass > mt.getMaxSize()){
					pattern.add(mass_sum);
					mass_sum = currentMass;
				}
				else
					mass_sum += currentMass;
			}
			else {
				if (mass_sum == 0)
					pattern.add(currentMass);
				else{
					pattern.add(mass_sum);
					pattern.add(currentMass);
					mass_sum = 0;
				}
			}
		}
		if (mass_sum > 0)
			pattern.add(mass_sum);
		
		tagSet[scanNum].add(pattern);
	}
	
	private void constructOridinaryTag(int startPos, HashSet[] tagSet){
		Pattern pattern = new Pattern();
		for (int i=startPos; i<startPos+PATTERN_LENGTH; i++){
			int currentMass = mt.getScaledMass(peptide.charAt(i));
			pattern.add(currentMass);
		}
		tagSet[scanNum].add(pattern);
	}
	
	
	
	/**
	 * constructs pattern(mass list) based on the confidence values
	 * @param tagSet
	 * 			the data structure for storing the obtained patterns
	 * @param gapMode
	 * 			true for gapped-tag mode, false for ordinary-tag mode
	 * */
	public void constructPattern(HashSet[] tagSet, boolean gapMode){
		Pattern pattern = new Pattern();
		int length = peptide.length();
		int confidence_sum = 0;
		int mass_sum = 0;
		int patternStartPos = 0, patternEndPos = 0; 
		
		for (int i=0; i<length; i++){
			int currentMass = mt.getScaledMass(peptide.charAt(i));
			if (confidence[i] < MIN_CONFIDENCE){
				confidence_sum += confidence[i];
				mass_sum += currentMass;
				if ((i==length-1) && (mass_sum<=mt.getMaxSize())){
					pattern.add(mass_sum);
					patternEndPos = i;
				}
			}
			else{
				if (mass_sum == 0){
					pattern.add(currentMass);
					patternEndPos = i;
					confidence_sum = 0;
					mass_sum = 0;
				}
				else if (mass_sum <= mt.getMaxSize()){
					pattern.add(mass_sum);
					pattern.add(currentMass);
					patternEndPos = i;
					confidence_sum = 0;
					mass_sum = 0;
				}
				else {
					if (pattern.length() == 0){
						pattern.add(currentMass);
						patternStartPos = i;
						patternEndPos = i;
						confidence_sum = 0;
						mass_sum = 0;
					}
					else {
						if (pattern.length() > 0){
							//tagSet[scanNum].add(pattern);
							//addNewTag(patternStartPos, patternEndPos, tagSet);
							if (gapMode)
								this.constructGapTag(patternStartPos, patternEndPos, tagSet);
							else
								this.constructOrdinaryTag(patternStartPos, patternEndPos, tagSet);
						}
							
						pattern = new Pattern();
						pattern.add(currentMass);
						patternStartPos = i;
						patternEndPos = i;
						confidence_sum = 0;
						mass_sum = 0;
					}
				}
			}	
		}
		if (pattern.length() > 0){
			//tagSet[scanNum].add(pattern);
			//addNewTag(patternStartPos, patternEndPos, tagSet);
			if (gapMode)
				this.constructGapTag(patternStartPos, patternEndPos, tagSet);
			else
				this.constructOrdinaryTag(patternStartPos, patternEndPos, tagSet);
		}
		
		pattern.clear();
		pattern = null;
	}
	
	/**
	 * constructs gapped tags of fixed length from a given part of the peptide
	 * @param startPos
	 * 			a position in the peptide
	 * @param endPos
	 * 			a position in the peptide, [startPos, endPos] indicates a range in the peptide
	 * @param tagSet
	 * 			to which the patterns are stored
	 * */
	private void constructGapTag(int startPos, int endPos, HashSet[] tagSet){
		if (endPos-startPos+1 < PATTERN_LENGTH)
			return;
		for (int i=startPos; i<=endPos-PATTERN_LENGTH+1; i++){
			int mass_sum = 0;
			Pattern pattern = new Pattern();
			for (int j=i; j<i+PATTERN_LENGTH; j++){
				int currentMass = mt.getScaledMass(peptide.charAt(j));
				if (confidence[j] < MIN_CONFIDENCE)
					mass_sum += currentMass;
				else{
					if (mass_sum > 0){
						pattern.add(mass_sum);
						mass_sum = 0;
					}
					pattern.add(currentMass);
				}
			}
			if (mass_sum > 0)
				pattern.add(mass_sum);
			tagSet[scanNum].add(pattern);
		}
	}
	
	/**
	 * constructs ordinary tags of fixed length from a given part of the peptide
	 * @param startPos
	 * 			a position in the peptide
	 * @param endPos
	 * 			a position in the peptide, [startPos, endPos] indicates a range in the peptide
	 * @param tagSet
	 * 			to which the patterns are stored
	 * */
	private void constructOrdinaryTag(int startPos, int endPos, HashSet[] tagSet){
		if (endPos-startPos+1 < PATTERN_LENGTH)
			return;
		for (int i=startPos; i<=endPos-PATTERN_LENGTH+1; i++){
			Pattern pattern = new Pattern();
			for (int j=i; j<i+PATTERN_LENGTH; j++){
				int currentMass = mt.getScaledMass(peptide.charAt(j));
				pattern.add(currentMass);
			}
			tagSet[scanNum].add(pattern);
		}
	}
	
	
	
	/**
	 * adds an ordinary tag (a substring of the peptide) into the tag set
	 * @param patternStartPos
	 * 			the position of the first amino acid of the tag
	 * @param patternEndPos
	 * 			the position of the last amino acid of the tag
	 * @param tagSet
	 * 			to where the tag is added
	 * */
	private void addNewTag(int patternStartPos, int patternEndPos, HashSet[] tagSet){
		Pattern pattern = new Pattern();
		for (int i=patternStartPos; i<=patternEndPos; i++){
			int aminoAcidMass = mt.getScaledMass(peptide.charAt(i));
			pattern.add(aminoAcidMass);
		}
		tagSet[scanNum].add(pattern);
	}
	
	
	
	public static void main(String[] args){
		String str = new String("LLLLIIIPWAIKLLLL");
		int[] confidence = {6,6,6,6,6, 100, 98, 96, 77, 45, 46, 6,6,6,6,6};
		String str2 = new String("IIIIIIIDDDDDDPWA");
		int[] confidence2 = {6,6,6,6,6, 100,60, 98,6,6,6,6,6, 96, 77, 45};
		String str3 = new String("ACDEFG");
		int[] confidence3 = {6, 6, 60, 46, 6,60};
		String str4 = new String("RKGEKLAPTLFK");
		int[] confidence4={7, 9, 9, 11, 11, 33, 91, 86, 14, 11, 6, 6};
		
		HashSet[] tagSet = new HashSet[3];
		for (int i=0; i<3; i++)
			tagSet[i] = new HashSet<Pattern>();
		
		
		MassTable mt = new MassTable();
		Peptide pt = new Peptide(0, str4, confidence4, mt);
		//pt.constructGapTag(13, tagSet);
		pt.constructPattern2(tagSet, true);
		
		
		
		
		
		System.out.println("tagSet[0]:"+tagSet[0]);
		System.out.println("tagSet[1]:"+tagSet[1]);
	}

}
