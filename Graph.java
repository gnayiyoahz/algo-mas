import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.stream.Collectors;
import java.lang.Math;

public class Graph {

  // Keep a fast index to nodes in the map
  private Map<Integer, Vertex> vertexNames;

  /**
   * Construct an empty Graph with a map. The map's key is the name of a vertex
   * and the map's value is the vertex object.
   */
  public Graph() {
    vertexNames = new HashMap<>();
  }

  /**
   * adds a vertex to the graph. Throws IllegalArgumentException if two vertices
   * with the same name are added.
   * 
   * @param v
   *          (Vertex) vertex to be added to the graph
   */
  public void addVertex(Vertex v) {
    if (vertexNames.containsKey(v.name))
      throw new IllegalArgumentException("Cannot create new vertex with existing name.");
    vertexNames.put(v.name, v);
  }

  /**
   * Gets a collection of all the vertices in the graph
   * 
   * @return (Collection<Vertex>) collection of all the vertices in the graph
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
  public Vertex getVertex(int name) {
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
  public void addEdge(int nameU, int nameV, Double cost) {
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
   * @param name
   *          (String) name of vertex u
   * @param name2
   *          (String) name of vertex v
   * @param cost
   *          (double) cost of the edge between vertex u and v
   */
  public void addUndirectedEdge(int name, int name2, double cost) {
    addEdge(name, name2, cost);
    addEdge(name2, name, cost);
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
    return Math.sqrt(Math.pow(ux - vx, 2) + Math.pow(uy - vy, 2));
  }

  /**
   * Computes euclidean distance between two vertices as described by their
   * coordinates
   * 
   * @param u
   *          (Vertex) vertex u
   * @param v
   *          (Vertex) vertex v
   * @return (double) distance between two vertices
   */
  public double computeEuclideanDistance(Vertex u, Vertex v) {
    return computeEuclideanDistance(u.x, u.y, v.x, v.y);
  }

  /**
   * Calculates the euclidean distance for all edges in the map using the
   * computeEuclideanCost method.
   */
  public void computeAllEuclideanDistances() {
    for (Vertex u : getVertices())
      for (Edge uv : u.adjacentEdges) {
        Vertex v = uv.target;
        uv.distance = computeEuclideanDistance(u.x, u.y, v.x, v.y);
      }
  }



  // STUDENT CODE STARTS HERE

  public void generateRandomVertices(int n) {
    vertexNames = new HashMap<>(); // reset the vertex hashmap
    
    // Your code here...
    for (int i = 0; i < n; i++){
	  int x = (int)(Math.random()*101);
	  int y = (int)(Math.random()*101);
	  Vertex myVertex = new Vertex(i,x,y);
      addVertex(myVertex);
	}
	for (int i = 0; i < n; i++){
	  for (int j = i+1; j < n; j++){
		addUndirectedEdge(i,j,1.0);
	  }
	}
    computeAllEuclideanDistances(); // compute distances
/*
	ArrayList<Vertex> allV = new ArrayList<Vertex>(getVertices()); 
	for (Vertex v:allV){
	  System.out.println(v);
	}
*/
  }

  public List<Edge> nearestNeighborTsp() {
	// initialization  
//	ArrayList<Vertex> allV = new ArrayList<Vertex>(getVertices()); // save the vertices 
//	for (Vertex currentV: allV){
	for (Vertex currentV:getVertices()){
	  currentV.distance = Double.POSITIVE_INFINITY;
	  currentV.known = false;
	  currentV.prev = null;
    }
	  
//	int random = (int) (Math.random()*allV.size());
	int r = vertexNames.size();
	int random = (int) (Math.random()*r);
	Vertex origin = getVertex(random); 
	System.out.println("Origin is "+origin);
	origin.prev = null;
	origin.distance = 0;
	origin.known = true;

	ArrayList<Edge> edgesAdded = new ArrayList<Edge>(); 
	for(int i = 0; i < r-1; i++){
	  edgesAdded = startWalking(edgesAdded,origin);
	  origin = edgesAdded.get(edgesAdded.size()-1).target; // recently added
//	  printE(edgesAdded);
	}
	edgesAdded.add(findEdge(origin,edgesAdded.get(0).source));
	  
	  
//	printE(edgesAdded); 
    return edgesAdded;
  }

	  
  public List<Edge> bruteForceTsp() {
//	ArrayList<Vertex> allV = new ArrayList<Vertex>(getVertices()); // save the vertices
	// permutation
	int size = vertexNames.size(); // number of vertices
	String intStr = "";
	for (int r = 0; r < size; r++){
	  intStr+=r; // intStr = 012...(n-1) for n vertices
	}
	ArrayList<String> sequences = new ArrayList<String>();
	sequences = permutation(sequences, intStr, 1, size-1);
	// start from position 1: permutation around a table -> 0 is fixed as the starting point 

/*	Tried but not used!
	for (int i = 0; i < sequences.size(); i++){	
	  sequences.set(i,sequences.get(i)+0);
	}	  
	sequences = removePalindrome(sequences);
	// palindromes like 012340 is calculated again in 043210, remove duplicates (half) 
	// removing them takes a long time, not used but I'll leave the code here
	Lessons learned: try not to go through a large list!
*/
	ArrayList<Edge> bfEdges = new ArrayList<Edge>();
	double tempMinD = Double.POSITIVE_INFINITY; 
	for (String s : sequences){ 
	  s += 0; // go back to starting point
	  double pathD = 0;
	  ArrayList<Edge> currentP = new ArrayList<Edge>();
	  for (int i = 0; i < size; i++){ //s.length() = size + 1
	    int j = i+1;
		int u = Character.getNumericValue(s.charAt(i));
		int v = Character.getNumericValue(s.charAt(j));
		Edge currentE = findEdge(getVertex(u),getVertex(v));
		currentP.add(currentE);  
		pathD += currentE.distance;
	  }
	  if (pathD < tempMinD){
		tempMinD = pathD;
		bfEdges = currentP;  
	  }
	}
  
//	printE(bfEdges);  
    return bfEdges;
  }

	
// general helper methods	
  // return edge from a to b	
  private Edge findEdge(Vertex from, Vertex to){
	List<Edge> neighborhood = from.adjacentEdges;
	for (Edge e : neighborhood){
	  if(e.source == from && e.target == to){
		return e;
	  }
	}
	return null;
  }	
		
  // to print out an arraylist for debugging
  private void printE(ArrayList<Edge> eee){
	System.out.println("------------Edges--------------");
	for(Edge e: eee){
	  if (e == null)
		  System.out.println("Some error occured.");
	  System.out.println(e+"with distance "+e.distance);
	}
	System.out.println("-------------------------------");
  }	
	
// helper method for nearestNeighborTsp	
  // start walking from some vertex	
  private ArrayList<Edge> startWalking(ArrayList<Edge> ea, Vertex x){
	Vertex origin = x;
	List<Edge> edges = origin.adjacentEdges;
	ArrayList<Edge> edgesAdded = ea;
	// some initialization
	Edge shortest = edges.get(0);
	double tempShort = shortest.distance;
	Vertex nextToGo = shortest.target;    
	for (Edge currentE: edges){
	  if (!currentE.target.known){
		shortest = currentE;
		tempShort = currentE.distance;
		nextToGo = currentE.target;
		break;
	  }
	}
	// decide on the next edge to add  
    for (Edge currentE: edges){
	  if (currentE.distance < tempShort && !currentE.target.known){
		tempShort = currentE.distance;
		shortest = currentE;
		nextToGo = currentE.target;
	  }
    } 
	  
	edgesAdded.add(shortest);
	nextToGo.known = true;

	return edgesAdded;
  }

// helper methods for bruteForceTsp
  private ArrayList<String> permutation(ArrayList<String> str, String s, int start, int end){
	ArrayList<String> strings = str;
	if (start == end){
	  strings.add(s);
	} else {
	  for (int i = start; i < end+1 ; i++){ 
		s = swap(s,start,i); // which one be the first
		strings = permutation(strings, s, start+1, end); // sort the substring with the first one fixed
		s = swap(s,start,i); // change back for the next i
	  } 
	}
	return strings;
  } 
	
  private String swap(String input, int i, int j){
	char temp;
	char[] charArray = input.toCharArray(); 
    temp = charArray[i] ; 
    charArray[i] = charArray[j]; 
    charArray[j] = temp; 
    return String.valueOf(charArray); 
  } 

	
/*	Tried but not used!!
  private ArrayList<String> removePalindrome(ArrayList<String> str){
	ArrayList<String> duplicate = new ArrayList<String>();
	ArrayList<String> strings = new ArrayList<>(str);
	for (String s:str){
	  if (duplicate.contains(s)){ 
		strings.remove(s);
	  } else {
		duplicate.add(getPalindrome(s)); 
	  }
	}
	return strings;  
  }
	
  private String getPalindrome(String s){
	String sp = "";
	int len = s.length();
	for (int i = len-1; i > -1; i--){
	  sp+=s.charAt(i);
	}  
	return sp;
  }
*/	

	
  // STUDENT CODE ENDS HERE



  /**
   * Prints out the adjacency list of the graph for debugging
   */
  public void printAdjacencyList() {
    for (int u : vertexNames.keySet()) {
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
}
