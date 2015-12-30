import java.util.Comparator;


public class SuffixPositionComparator implements Comparator<SuffixPosition>{
	
	public int compare(SuffixPosition s1, SuffixPosition s2){
		
		if (s1.getSeqNum() == s2.getSeqNum())
			return s1.getPosInSeq() - s2.getPosInSeq();
		return s1.getSeqNum() - s2.getSeqNum();
	}

}
