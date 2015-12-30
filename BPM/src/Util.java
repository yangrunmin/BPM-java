import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Map;


public class Util {

	static String printTreeForGraphViz(SuffixTree tree) {
		return printTreeForGraphViz(tree, true);
	}
	
	
	/**
	 * Generates a .dot format string for visualizing a suffix tree.
	 * 
	 * @param tree
	 *            The tree for which we are generating a dot file.
	 * @return A string containing the contents of a .dot representation of the
	 *         tree.
	 */
	static String printTreeForGraphViz(SuffixTree tree, boolean printSuffixLinks) {
		LinkedList<Node> stack = new LinkedList<Node>();
		stack.add(tree.getRoot());
		Map<Node, Integer> nodeMap = new HashMap<Node, Integer>();
		nodeMap.put(tree.getRoot(), 0);
		int nodeId = 1;

		StringBuilder sb = new StringBuilder(
				"\ndigraph suffixTree{\n node [shape=circle, label=\"\", fixedsize=true, width=0.1, height=0.1]\n");

		while (stack.size() > 0) {
			LinkedList<Node> childNodes = new LinkedList<Node>();
			for (Node node : stack) {
				for (int i=0; i<21; i++){
					Edge edge = node.getEdge(i);
					if (edge == null)
						continue;
					int id = nodeId++;
					if (edge.hasEndNode()){
						childNodes.push(edge.getEndNode());
						nodeMap.put(edge.getEndNode(), id);
					}
					sb.append(nodeMap.get(node)).append(" -> ").append(id).append(" [label=\"");
					sb.append(edge.toString());
					sb.append("\"];");//sb.append("\"];\n");
					sb.append(" ["+edge.getBeginIndex()+","+edge.getEndIndex()+"]");
					if (!edge.hasEndNode())
						sb.append(" suffixPosition:"+((LeafEdge)edge).getSuffixPosition());
					sb.append("\n");
				}
			}
			stack = childNodes;
		}
		if(printSuffixLinks){
			// loop again to find all suffix links.
			sb.append("edge [color=red]\n");
			for (Map.Entry<Node, Integer> entry : nodeMap.entrySet()) {
				Node n1 = entry.getKey();
				int id1 = entry.getValue();

				if (n1.hasSuffixNode()) {
					Node n2 = n1.getSuffixNode();
					Integer id2 = nodeMap.get(n2);
					sb.append(id1).append(" -> ").append(id2).append(" ;\n");
				}
			}
		}
		sb.append("}");
		return (sb.toString());
	}
}
