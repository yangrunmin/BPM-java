import java.io.IOException;


public class SuffixTreeTest {
	
	public static void main(String[] args) throws IOException{
		MassTable mt = new MassTable();
		mt.generateMassTable();
		DatabaseFileHandler df = new DatabaseFileHandler();
		ProteinDatabase pd = df.loadDatabase("test2.fasta");
		
		SuffixTree st = new SuffixTree(pd.getSequence(), pd);
		System.out.println(Util.printTreeForGraphViz(st));
		st.search("IVQ");
		System.out.println("hello world!");
	}

}
