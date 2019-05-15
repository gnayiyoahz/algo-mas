import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.io.IOException;
import java.io.FileReader;
import java.io.BufferedReader;
import java.lang.Math;
import java.util.ArrayList;

public class Dijkstra {

  // Keep a fast index to nodes in the map
  private Map<String, Vertex> vertexNames;

  /**
   * Construct an empty Dijkstra with a map. The map's key is the name of a vertex
   * and the map's value is the vertex object.
   */
  public Dijkstra() {
    vertexNames = new HashMap<String, Vertex>();
  }

  /**
   * Adds a vertex to the dijkstra. Throws IllegalArgumentException if two vertices
   * with the same name are added.
   * 
   * @param v
   *          (Vertex) vertex to be added to the dijkstra
   */
  public void addVertex(Vertex v) {
    if (vertexNames.containsKey(v.name))
      throw new IllegalArgumentException("Cannot create new vertex with existing name.");
    vertexNames.put(v.name, v);
  }

  /**
   * Gets a collection of all the vertices in the dijkstra
   * 
   * @return (Collection<Vertex>) collection of all the vertices in the dijkstra
   */
  public Collection<Vertex> getVertices() {
    return vertexNames.values();
  }

  /**
   * Gets the vertex object with the given name
   * 
   * @param name
   *          (String) name of the vertex object requested
   * @return (Vertex) vertex object associated with the name
   */
  public Vertex getVertex(String name) {
    return vertexNames.get(name);
  }

  /**
   * Adds a directed edge from vertex u to vertex v
   * 
   * @param nameU
   *          (String) name of vertex u
   * @param nameV
   *          (String) name of vertex v
   * @param cost
   *          (double) cost of the edge between vertex u and v
   */
  public void addEdge(String nameU, String nameV, Double cost) {
    if (!vertexNames.containsKey(nameU))
      throw new IllegalArgumentException(nameU + " does not exist. Cannot create edge.");
    if (!vertexNames.containsKey(nameV))
      throw new IllegalArgumentException(nameV + " does not exist. Cannot create edge.");
    Vertex sourceVertex = vertexNames.get(nameU);
    Vertex targetVertex = vertexNames.get(nameV);
    Edge newEdge = new Edge(sourceVertex, targetVertex, cost);
    sourceVertex.addEdge(newEdge);
  }

  /**
   * Adds an undirected edge between vertex u and vertex v by adding a directed
   * edge from u to v, then a directed edge from v to u
   * 
   * @param nameU
   *          (String) name of vertex u
   * @param nameV
   *          (String) name of vertex v
   * @param cost
   *          (double) cost of the edge between vertex u and v
   */
  public void addUndirectedEdge(String nameU, String nameV, double cost) {
    addEdge(nameU, nameV, cost);
    addEdge(nameV, nameU, cost);
  }

  /**
   * Computes the euclidean distance between two points as described by their
   * coordinates
   * 
   * @param ux
   *          (double) x coordinate of point u
   * @param uy
   *          (double) y coordinate of point u
   * @param vx
   *          (double) x coordinate of point v
   * @param vy
   *          (double) y coordinate of point v
   * @return (double) distance between the two points
   */
  public double computeEuclideanDistance(double ux, double uy, double vx, double vy) {
	  //TODO
	  return Math.sqrt(Math.pow(ux-vx, 2) + Math.pow(uy-vy, 2));
  }

  /**
   * Calculates the euclidean distance for all edges in the map using the
   * computeEuclideanDistance method.
   */

  public void computeAllEuclideanDistances() {
	  // TODO
	  ArrayList<Vertex> allV = new ArrayList<Vertex>(getVertices());
	  int vNum = vertexNames.values().size();
	  for (int i = 0; i < vNum; i++){
		  Vertex currentV = allV.get(i);
		  List<Edge> edges = currentV.adjacentEdges;
		  int d = edges.size();
		  for (int j = 0; j < d; j++ ){
			  Edge currentE = edges.get(j);
			  Vertex v = currentE.source;
			  Vertex u = currentE.target;
			  currentE.distance = computeEuclideanDistance(v.x,v.y,u.x,u.y);
//			  System.out.println(currentE);
		  }
//		  System.out.println(currentV);
	  }
			
  }

  /**
   * Dijkstra's Algorithm. 
   * 
   * @param s
   *          (String) starting city name
   */
  public void doDijkstra(String s) {
	  // TODO
	  ArrayList<Vertex> allV = new ArrayList<Vertex>(getVertices()); 
	  for (Vertex currentV: allV){
		  currentV.distance = Double.POSITIVE_INFINITY;
		  currentV.known = false;
		  currentV.prev = null;
	  }
	  ArrayList<Vertex> pathFound = new ArrayList<Vertex>();
	  Vertex origin = getVertex(s);
	  origin.prev = null;
	  origin.distance = 0;
	  origin.known = true;
	  pathFound.add(origin);
	  pathFound = startWalking(pathFound, origin);  
	  printPF(pathFound);
	  
	  
	  
	  boolean allKnown = check(allV);
//	  printPF(pathFound);
	  while(!allKnown){
		  Vertex x = findMax(pathFound);// find min among known
		  if (x != null){
			pathFound = startWalking(pathFound, x);
		  }
		  printPF(pathFound);
		  allKnown = check(allV);
//		  printPF(pathFound);
	  } 

  }
	
	// to printPF for debugging
  private void printPF(ArrayList<Vertex> pf){
	System.out.println("--------------IN PF------------");
	for(Vertex v :pf){
	  System.out.println(v+"with distance "+v.distance);
	}
	System.out.println("-------------------------------");
  }

	//see if every vertex in an ArrayList is known
  private boolean check(ArrayList<Vertex> vvv){
	  boolean allknown = true;
	  for(Vertex v: vvv){
		  if (v.known == false){
			  allknown = false;
			  return false;
		  }
	  }
	  return allknown;
  }
	
  private ArrayList<Vertex> startWalking(ArrayList<Vertex> pf, Vertex x){
	  Vertex origin = x;
	  List<Edge> edges = origin.adjacentEdges;
	  ArrayList<Vertex> pathFound = pf;  
	  for (Edge currentE: edges){
		  Vertex someNeighbor = currentE.target;
		  if (!someNeighbor.known){
			  double updatedD = x.distance+currentE.distance;
			  boolean seen = pathFound.contains(someNeighbor);
			  if(updatedD < someNeighbor.distance){
//				  System.out.println("--------------------------------");
//				  System.out.println("      original D "+ someNeighbor.distance);
//				  System.out.println("      now gonna be "+ updatedD);
				  someNeighbor.distance = updatedD;
				  someNeighbor.prev = origin;
				  pathFound.add(someNeighbor); 
				  if(seen){
//					  System.out.println("                seen.");
					  pathFound.remove(someNeighbor);
				  }
//				  System.out.println(someNeighbor +"'s distance updated. Prev is "+origin);
			  }
		  }
	  }  
	  Vertex nextMin = findMin(pathFound);
	  if (nextMin != null && !nextMin.known){
		nextMin.known = true;
//	    System.out.println(nextMin+" is now known.");
	  }
	  return pathFound;
  }
	
  
  // find the maximum among all verices that are KNOWN
  // i.e. the vertex recently marked known
  private Vertex findMax(ArrayList<Vertex> pathFound){
	Vertex tempMax = pathFound.get(0);
	// newly connected vertex
	boolean allKnown = check(pathFound);
	if (allKnown)
	  return null;
    for (Vertex v: pathFound){
	  if (v.distance > tempMax.distance && v.known)
	    tempMax = v;
    }
	return tempMax;
  }
	
	
	
  // find minimum distance
  private Vertex findMin(ArrayList<Vertex> pathFound){
	Vertex tempShort = pathFound.get(0);
	for (Vertex v: pathFound){
	  if (!v.known){
	    tempShort = v;
		break;
	  }
    } 
	// newly connected vertex
	boolean allKnown = check(pathFound);
	if (allKnown)
	  return null;
    for (Vertex v: pathFound){
	  if (v.distance < tempShort.distance && !v.known)
	    tempShort = v;
    }
	return tempShort;
  }
	  

  /**
   * Returns a list of edges for a path from city s to city t. This will be the
   * shortest path from s to t as prescribed by Dijkstra's algorithm
   * 
   * @param s
   *          (String) starting city name
   * @param t
   *          (String) ending city name
   * @return (List<Edge>) list of edges from s to t
   */
  public List<Edge> getDijkstraPath(String s, String t) {
	computeAllEuclideanDistances();
	doDijkstra(s);
//	System.out.println("doDijkstra is done!!!!!");
    // TODO
    Vertex start = getVertex(s);
	Vertex to = getVertex(t);
    ArrayList<Edge> shortest = new ArrayList<Edge>();
	Vertex from = to.prev;
	while(to != start){
//	  System.out.println("from "+ from);
//	  System.out.println("to "+to);
	  shortest.add(findEdge(from,to));
	  to = from;
	  from = from.prev;
//		System.out.println("printing path:");
//		for(Edge e:shortest){
//			System.out.println(e);
//		}	
	}
  
    return shortest;
  }
	
  private Edge findEdge(Vertex from, Vertex to){
	List<Edge> neighborhood = to.adjacentEdges;
//	  System.out.println("printing neighborhood:");
//	  for(Edge e:neighborhood){
//			System.out.println(e);
//		}
	  
	for (Edge e : neighborhood){
	  if(e.source == to && e.target == from){
		return e;
	  }
	}
	return null;
  }

  /**
   * Prints out the adjacency list of the dijkstra for debugging
   */
  public void printAdjacencyList() {
    for (String u : vertexNames.keySet()) {
      StringBuilder sb = new StringBuilder();
      sb.append(u);
      sb.append(" -> [ ");
      for (Edge e : vertexNames.get(u).adjacentEdges) {
        sb.append(e.target.name);
        sb.append("(");
        sb.append(e.distance);
        sb.append(") ");
      }
      sb.append("]");
      System.out.println(sb.toString());
    }
  }


  /** 
   * A main method that illustrates how the GUI uses Dijkstra.java to 
   * read a map and represent it as a graph. 
   * You can modify this method to test your code on the command line. 
   */
  public static void main(String[] argv) throws IOException {
    String vertexFile = "cityxy.txt"; 
    String edgeFile = "citypairs.txt";

    Dijkstra dijkstra = new Dijkstra();
    String line;

    // Read in the vertices
    BufferedReader vertexFileBr = new BufferedReader(new FileReader(vertexFile));
    while ((line = vertexFileBr.readLine()) != null) {
      String[] parts = line.split(",");
      if (parts.length != 3) {
        vertexFileBr.close();
        throw new IOException("Invalid line in vertex file " + line);
      }
      String cityname = parts[0];
      int x = Integer.valueOf(parts[1]);
      int y = Integer.valueOf(parts[2]);
      Vertex vertex = new Vertex(cityname, x, y);
      dijkstra.addVertex(vertex);
    }
    vertexFileBr.close();

    BufferedReader edgeFileBr = new BufferedReader(new FileReader(edgeFile));
    while ((line = edgeFileBr.readLine()) != null) {
      String[] parts = line.split(",");
      if (parts.length != 3) {
        edgeFileBr.close();
        throw new IOException("Invalid line in edge file " + line);
      }
      dijkstra.addUndirectedEdge(parts[0], parts[1], Double.parseDouble(parts[2]));
    }
    edgeFileBr.close();

    // Compute distances. 
    // This is what happens when you click on the "Compute All Euclidean Distances" button.
    dijkstra.computeAllEuclideanDistances();
    
    // print out an adjacency list representation of the graph
    dijkstra.printAdjacencyList();

    // This is what happens when you click on the "Draw Dijkstra's Path" button.

    // In the GUI, these are set through the drop-down menus.
    String startCity = "SanFrancisco";
    String endCity = "Boston";

    // Get weighted shortest path between start and end city. 
    List<Edge> path = dijkstra.getDijkstraPath(startCity, endCity);
    
    System.out.print("Shortest path between "+startCity+" and "+endCity+": ");
    System.out.println(path);
  }

}
