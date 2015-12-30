import java.util.Vector;


public class Preprocessing {
	private SuffixTree suffixTree;
	private ProteinDatabase proteinDatabase;
	private Vector<Vector<Integer>> peptideTable;
	private DoSearch2 doSearch;
	
	public Preprocessing(SuffixTree st, ProteinDatabase pd, DoSearch2 doSearch, Vector<Vector<Integer>> peptideTable){
		suffixTree = st;
		proteinDatabase = pd;
		this.peptideTable = peptideTable;
		this.doSearch = doSearch;
	}
	
	public void computeLeafPeptideMass(LeafEdge edge){
		SuffixPosition sp = edge.getSuffixPosition();
		sp.compPeptidePos2(peptideTable, proteinDatabase);
		String pepString = sp.getPeptideSequence(proteinDatabase);
		edge.setPeptideSequence(pepString);
		int mass = doSearch.computePeptideMass(pepString);
		edge.setPeptideMass(mass);
		sp = null;
		pepString = null;
	}
	
	public void computeMassRangeForNodes(Node node){
		for (int i=0; i<21; i++){
			Edge edge = node.getEdge(i);
			if (edge == null) continue;
			Node childNode = edge.getEndNode();
			if (childNode == null)
				computeLeafPeptideMass((LeafEdge)edge);
			else
				computeMassRangeForNodes(childNode);
		}
		
		int min = Integer.MAX_VALUE, max = -1;
		for (int i=0; i<21; i++){
			Edge edge = node.getEdge(i);
			if (edge == null) continue;
			Node childNode = edge.getEndNode();
			if (childNode == null){
				int mass = ((LeafEdge) edge).getPeptideMass();
				if (min > mass)
					min = mass;
				if (max < mass)
					max = mass;
			}
			else{
				if (min > childNode.getMinPeptideMass())
					min = childNode.getMinPeptideMass();
				if (max < childNode.getMaxPeptideMass())
					max = childNode.getMaxPeptideMass();
			}
		}
		
		node.setMinPeptideMass(min-doSearch.getScaledMassTolerance());
		node.setMaxPeptideMass(max+doSearch.getScaledMassTolerance());
	}
}
