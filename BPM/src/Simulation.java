import java.util.Random;


public class Simulation {

	public static int SubseqLength = 20;
	private MassTable mt;
	private ProteinDatabase proteinDatabase;
	
	/**
	 * gets a pattern from the protein database
	 * @param patternStyle
	 * 			an integer which is used to indicate how to generate the pattern, right now it can be 1,2,3 or 4
	 * @param random
	 * 			an instance of java.util.Random
	 * */
	public Pattern getPattern(int patternStyle, Random random){
		Pattern pattern = new Pattern();
		String randomSeq = getRandomSubseq(random);
		randomSeq = handleUnknownCharacter(randomSeq);
		
		if (patternStyle == 1){
			for (int i=0; i<SubseqLength; i++){
				int mass = mt.getScaledMass(randomSeq.charAt(i));
				pattern.add(mass);
			}
		}
		
		return pattern;
	}
	
	/**
	 * gets a random subsequence with length <SubseqLength>
	 * @param random
	 * 			an instance of java.util.Random
	 * @return a random subsequence from the protein database
	 * */
	private String getRandomSubseq(Random random){
		return proteinDatabase.getRandomSeq(SubseqLength, random);
	}
	
	public void runSimulation(int patternStyle){
		
	}
	
	/**
	 * replaces unknown amino acids with 'A'
	 * @param str 
	 * 			a given amino acid sequence
	 * @return a string with all unknown characters replaced by 'A'
	 * */
	private String handleUnknownCharacter(String str){
		str = str.replace('B', 'A');
		str = str.replace('J', 'A');
		str = str.replace('O', 'A');
		str = str.replace('U', 'A');
		str = str.replace('X', 'A');
		str = str.replace('Z', 'A');
		return str;
	}
}
