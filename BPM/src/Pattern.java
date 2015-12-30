import java.util.ArrayList;
import java.util.Vector;



/**
 * Defines a pattern which consists of several integers
 * massList: an ArrayList containing the integers
 */
public class Pattern {
	private ArrayList<Integer> massList;
	
	/**
	 * Constructor, creates a new Pattern object, with empty list 
	 */
	public Pattern(){
		massList = new ArrayList<Integer>();
	}
	
	/**
	 * Adds elements into massList
	 * @param value
	 * 			the value of the element to be added 
	 */
	public void add(int value){
		massList.add(value);
	}
	
	/**
	 * Removes specified elements from massList
	 * @param index
	 * 			the position of the element to be removed
	 * */
	public void remove(int index){
		massList.remove(index);
	}
	
	public void clear(){
		massList.clear();
	}
	
	/**
	 * Gets the element at specific position
	 * @param index
	 * 			the given position
	 * @return the required element 
	 */
	public int get(int index){
		return massList.get(index);
	}
	
	/**
	 * get the number of elements in the Pattern
	 * @return #elements in the pattern
	 */
	public int length(){
		return massList.size();
	}
	
	/**
	 * gets a sub pattern containing elements of this pattern in the range [startIndex, endIndex)
	 * @param startIndex
	 * 			the number indicating the beginning range (include)
	 * @param endIndex
	 * 			the number indicating the end of the range(not included)
	 * @return the required sub pattern
	 * */
	public Pattern subPattern(int startIndex, int endIndex){
		Pattern pattern = new Pattern();
		for (int i=startIndex; i<endIndex; i++)
			pattern.add(this.get(i));
		return pattern;
	}
	
	
	/**
	 * splits a pattern into smaller ones with fix length
	 * @param fixedPatternLength
	 * 			a given integer used to control the length of each obtained sub-pattern
	 * @return a set of sub-patterns of fixed length extracted from this pattern
	 * @example Suppose this pattern is [1,2,3,4,5] and fixedPatternLength=3, then the result should
	 * 			be a vector containing 3 elements [[1,2,3],[2,3,4],[3,4,5]] 
	 * */
	public Vector<Pattern> splitPattern(int fixedPatternLength){
		Vector<Pattern> newPatternList = new Vector<Pattern>();
		if (length() <= fixedPatternLength)
			newPatternList.add(this);
		else {
			for (int i=0; i<=length()-fixedPatternLength; i++){
				Pattern newPattern = subPattern(i, i+fixedPatternLength);
				newPatternList.add(newPattern);
			}
		}
		
		return newPatternList;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result
				+ ((massList == null) ? 0 : massList.hashCode());
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		Pattern other = (Pattern) obj;
		if (massList == null) {
			if (other.massList != null)
				return false;
		} else if (!massList.equals(other.massList))
			return false;
		return true;
	}

	
	public String toString(){
		return massList.toString();
	}
}
