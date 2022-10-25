/*
* AUTHOR: Kevin Nisterenko
* FILE: DGraph.java
* ASSIGNMENT: Programming Assignment 11 - The Traveling Salesman Problem
* COURSE: CSc 210; Fall 2021
* PURPOSE: This program defines the DGraph class, which represents a directed, weighted
* graph, where the vertices are integers and the weight of the edges are doubles. To 
* store this information the underlying structure is that of an adjacency list, where the 
* vertices are the slots of an arraylist and each slot contains a linkedlist with the edges
* to other vertices, which is a class itself. The graph is constructed with a fixed 
* size/number of vertices given to the constructor. There are methods to add new edges, get
* the cost of a path in the graph and also check the weight of edges and get the minimum
* connection between vertices. These are all used to solve the TSP problem and help with three
* implementations of solutions a heuristic, a brute force backtrack and a solution I defined,
* which is a mix of both heuristics and backtracking to decrease the number of permutations
* checked in the backtracking and get a good enough path in a faster time. 
*
* There are no inputs for this specific file. 
*/

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

public class DGraph {
	private int numVertices;
	private ArrayList<LinkedList<Edge>> adjList = new ArrayList<>();
	
	/*
	 * Provides the Edge class, which represents a path between two vertices, 
	 * since it will be added to the adjacency list of the origin vertex, it 
	 * contains the information on the weight of the edge and the destination 
	 * vertex.
	 */
	private class Edge {
		private int destinationVertex;
		private double weight; 
		
		/*
		 * Constructor for an Edge object, it takes a destination vertex and also 
		 * the weight of the path and sets the attributes to these parameters.
		 * 
		 * @param vertex, integer representing destination vertex in the edge
		 * @param weight, double representing the weight of the edge
		 */
		public Edge(int vertex, double weight) {
			this.destinationVertex = vertex;
			this.weight = weight;
		}
		
	}
	
	/*
	 * Constructor for a DGraph object, it takes the number of vertices and 
	 * creates the empty adjacency list so that edges between the vertices can
	 * be added later. 
	 * 
	 * @param numVertices, integer representing the number of vertices in the 
	 * graph
	 */
	public DGraph(int numVertices) { 
		this.numVertices = numVertices; 
		for (int i = 0; i < numVertices; i++)
			this.adjList.add(new LinkedList<Edge>());
	} 
	
	/*
	 * Public method used to add edges between vertices in the graph, 
	 * it takes from an to vertices and the weight for the edge between 
	 * them. A new Edge object is created and then added to the correct 
	 * slot of the adjacency list. 
	 * 
	 * @param from, integer representing the index of the origin/from vertex
	 * of the edge, also the slot in the list where the edge is to be added
	 * @param to, integer representing the index of the destination vertex
	 * @param weight, double representing the weight of the edge between the 
	 * two given vertices in the specific direction
	 */
	public void addEdge(int from, int to, double weight) { 
		this.adjList.get(from).add(new Edge(to, weight)); 
	}
	
	/*
	 * Private helper method used to calculate the distance between two
	 * given vertices. It iterates over the inner linked list of the 
	 * origin/from vertex and find the edge containing the destination 
	 * vertex. 
	 * 
	 * @param vertex, integer representing the index of the origin/from vertex
	 * of the edge, also the slot in the list where the edge can be found
	 * @param destination, integer representing the index of the 
	 * destination vertex
	 * @return closeWeight, double representing the weight of the 
	 * edge between two vertices in the specific given direction
	 */
	private double edgeWeight(int vertex, int destination) {
		List<Edge> connections = this.adjList.get(vertex);
		double closeWeight = 0;
		// Iterates over all edges, if the edge connects to the destination,
		// we can set the weight to it
		for (int i = 0; i < connections.size(); i++) {
			if (connections.get(i).destinationVertex == destination) {
				closeWeight = connections.get(i).weight;
			}
		}
		return closeWeight;
	}
	
	/*
	 * Private helper method that takes a path of travel of vertices in the graph 
	 * and calculates the total cost to travel said path. It uses a for loop to 
	 * iterate over the given path and it grabs the weight of the path between 
	 * every two vertices by calling the edgeWeight method. Lastly, fence-posting 
	 * is used to connect back the path to the beginning in a cycle manner. 
	 * 
	 * @param path, list of integer where every element is a vertex in the graph
	 * @return cost, double representing the total cost of travel in a given cycle
	 * on the graph
	 */
	private double getPathCost(List<Integer> path) {
		double cost = 0.0;
		// Iterate over the path, for every two vertices get the cost between 
		// them
		for (int i = 0; i < path.size()-1; i++) {
			cost += edgeWeight(path.get(i)-1, path.get(i+1)-1);
		}
		// Connect the last back to the first to close the cycle
		int last = path.get(path.size()-1)-1;
		cost += edgeWeight(last, 0);
		return cost;
	}
	
	/*
	 * Private helper method that iterates over the edges of a given vertex and 
	 * returns an unvisited edge with the minimum distance between the vertices. 
	 * A for loop is used to calculate and set the return value to this minimum edge
	 * and a boolean array of the visited vertices is also given. 
	 * 
	 * @param connections, linked list of edge objects that represent all the connected
	 * vertex of a specific vertex in the adjacency list
	 * @param visited, boolean array representing the visited vertices in the path, 
	 * true if visited, false otherwise, and since the vertices are represented  with 
	 * integers, we use the index to identify the vertices
	 * @return retVal, Edge object representing the unvisited edge with the minimum 
	 * weight/distance in the list of edges of a vertex
	 */
	private Edge getMinConnection(List<Edge> connections, boolean[] visited) {
		double min = Double.MAX_VALUE;
		Edge retVal = connections.get(0);
		// Iterate over every edge for the current vertex
    	for (int i = 0; i < connections.size(); i++) {
    		// If vertex has still not been visited
    		if (!visited[connections.get(i).destinationVertex]) {
    			// If the path is the smallest one
    			if (min > connections.get(i).weight) {
    				// Save this path and destination vertex
    				min = connections.get(i).weight;
    				retVal = connections.get(i);
    			}
    		}
    	}
    	return retVal;
	}
	
	/*
	 * Private helper method used to round a double to one decimal place. It 
	 * takes a double and returns another double, which is the rounded version 
	 * of the given number.
	 * 
	 * @param val, double representing the cost of a path
	 * @return double, given value rounded down to one decimal place
	 */
	private double roundNum(double val) {
	    int weight = (int) Math.pow(10, 1); // Power of 10, and 1, means 1 decimal
	    return (double) Math.round(val * weight) / weight;
	}
	
	/*
	 * Heuristic method to solve the TSP problem by using the Nearest Neighbor
	 * Algorithm idea. We start from the starting vertex, find the minimum 
	 * possible path to another vertex and move to that one, we do this until we 
	 * hit all vertices. This gives a possible solution in quadratic time, although 
	 * since it only performs a local optimum check, it may not be the best possible 
	 * solution to the graph, but a close one nonetheless. The cost of the found path 
	 * is then returned and the path may be printed if indicated by the given boolean 
	 * for printing. 
	 * 
	 * @param start, starting vertex of the path
	 * @param print, boolean to indicate whether or not we want to print the 
	 * cost and path for this run of the algorithm
	 * @return cost, total cost of the found path
	 */
	public double heuristicTSP(int start, boolean print) {
		boolean[] visited = new boolean[numVertices];
		int currVertex = start;
		double cost = 0;
		List<Integer> path = new ArrayList<Integer>();
		path.add(start+1);
		// Move until we have every vertex in our path
		for (int i = 0; i < this.numVertices - 1; i++) {
			// Visit this vertex
			visited[currVertex] = true;
			LinkedList<Edge> connections = this.adjList.get(currVertex);
			// Get the smallest edge/connection
			Edge minEdge = getMinConnection(connections, visited);
			double min = minEdge.weight;
			currVertex = minEdge.destinationVertex;
	    	// Move in the smallest direction
	    	path.add(currVertex+1); // +1 to go back to file number system
	    	// Visit the new minimum vertex and increment the cost
	    	visited[currVertex] = true;
	    	cost += min;
		}
		// Connect back to start to close the cycle 
		List<Edge> connections = this.adjList.get(currVertex);
		for (int i = 0; i < connections.size(); i++) {
			if (connections.get(i).destinationVertex == start) {
				cost += connections.get(i).weight;
			}
		}
		cost = roundNum(cost);
		if (print) System.out.println("cost = " + cost + ", visitOrder = " + path);
		return cost;
	}
	
	/*
	 * Brute force backtracking method to solve the traveling salesman problem. It 
	 * uses recursive backtracking to find all permutations with the given vertex at 
	 * the start of the path, and then keep the best possible one in the whole 
	 * graph, however at a factorial cost. It calls the backtrackHelper method to
	 * actually perform the recursion and path finding.
	 * 
	 * @param start, starting vertex of the path
	 * @param print, boolean to indicate whether or not we want to print the 
	 * cost and path for this run of the algorithm
	 * @return cost, total cost of the found path
	 */
	public double backtrackTSP(int start, boolean print) {
		// Keeps track of *which* vertices we have visited
		boolean[] visited = new boolean[numVertices];
		// Initial cost is 0
		double cost = 0;
		// Maximum possible cost so we can compare and move down
		double result = Double.MAX_VALUE;
		// So we can keep tracking how our recursive function is getting called
		LinkedList<Integer> path = new LinkedList<>();
		path.add(start+1);
		List<Integer> savedPath = new LinkedList<>();
		cost = backtrackHelper(start, visited, cost, result, path, savedPath);
		cost = roundNum(cost);
		// Since we have been adding our possible paths, now we can just get the last cycle
		List<Integer> minPath = savedPath.subList(savedPath.size()-this.numVertices,
				savedPath.size());
		if (print) System.out.println("cost = " + cost + ", visitOrder = " + minPath);
		return cost;
	}
	
	/*
	 * Private recursive method to solve the TSP problem with brute-force backtracking. It 
	 * first checks if the given path is already a complete cycle, and tests if it is the 
	 * minimum possible cycle and returns the minimum cost. If this base case is not hit, 
	 * then we must visit vertices to complete the path. To do so, a for loop is used to 
	 * iterate over every vertex and check if said vertex has not been visited. Then we visit 
	 * it and do the recursive call to test the paths with this vertex added. After the 
	 * recursive call, we unvisit the vertex so we can "go back", and thus achieve backtracking.
	 * 
	 * @param vertex, integer representing the current vertex in the call
	 * @param visited, boolean array representing the visited vertices in the path, 
	 * true if visited, false otherwise, and since the vertices are represented  with 
	 * integers, we use the index to identify the vertices 
	 * @param cost, double representing the total cost of travel in a given cycle
	 * on the graph
	 * @param minCost, double representing the current minimum cost of travel in 
	 * on the graph  
	 * @param path, list of integer where every element is a vertex in the graph
	 * @param resPath, list of integer used to record the traveled paths, the last/optimum
	 * path will always be the last one in it (of size numVertices)
	 * @return cost, double representing the total minimum cost of travel in the graph
	 */
	private double backtrackHelper(int vertex, boolean[] visited, double cost,
			double minCost, List<Integer> path, List<Integer> savedPath) {
		// Grab closing distance (last vertex to start)
		double closeWeight = edgeWeight(vertex, 0);
		// If we hit every vertex and we can connect back to the start, 
		// we hit the base case and we can return the smallest cost
		if (path.size() == this.numVertices && closeWeight > 0) {
			if (minCost > cost + closeWeight) {
				minCost = cost + closeWeight;
				savedPath.addAll(path);
			} 
			return minCost;
		}
		// This is where the backtracking happens, we iterate over every vertex 
		for (int i = 0; i < this.adjList.size(); i++) {
			// If we haven't yet visited the vertex and we have a connection to it
			if (!visited[i] && edgeWeight(vertex, i) > 0) {
				// Mark it as visited, do the recursive call to find the path
				visited[i] = true;
				path.add(i+1);
				minCost = backtrackHelper(i, visited, cost + edgeWeight(vertex, i),
						minCost, path, savedPath);
				// Reset the visited status so we can backtrack
				visited[i] = false;
				path.remove(path.size()-1);
			}
		}
		return minCost;
	} 
    
    /*
	 * My solution to the TSP, it uses ideas from both the backtracking and heuristic
	 * methods to solve the traveling salesman problem. It uses recursive backtracking
	 * to quickly find permutations with the given vertex at the start of the path. However
	 * it will not find a factorial number of permutations since there is an added constraint 
	 * to the paths explored, it will look for the minimum weight neighbor similar to heuristics
	 * and then keep the best possible path with this. It calls the mineHelper method to
	 * actually perform the recursion and path finding.
	 * 
	 * @param start, starting vertex of the path
	 * @param print, boolean to indicate whether or not we want to print the 
	 * cost and path for this run of the algorithm
	 * @return cost, total cost of the found path
	 */
    public double mineTSP(int start, boolean print) {
    	boolean[] visited = new boolean[numVertices];
    	// Set the resulting path and the total path formed in the recursion
    	List<Integer> resPath = new ArrayList<Integer>();
    	List<Integer> path = new ArrayList<>();
    	path.add(start+1);
    	// Set the starting to the maximum possible double to then go down and find the minimum
    	double cost = Double.MAX_VALUE;
    	cost = mineHelper(start, visited, path, cost, resPath);
    	cost = roundNum(cost);
    	// Gets the last path of size numVertices (which will be the optimum one)
    	resPath = resPath.subList(resPath.size()-this.numVertices, resPath.size());
    	if (print) if (print) System.out.println("cost = " + cost + ", visitOrder = " + resPath);
    	return cost;
    }
	
    /*
     * Private recursive method to solve the TSP problem with a mix of heuristics and
     * backtracking. It first checks if the given path is already a complete cycle, and
     * tests if it is the minimum possible cycle and returns the minimum cost. If this
     * base case is not hit, then we must visit vertices to complete the path. To do so,
     * a for loop is used to iterate over every vertex and check if said vertex has not been
     * visited and if that vertex is the minimum distance neighbor (added constraint that 
     * decreases the number of permutations/paths we visit). Then we visit it and do the 
     * recursive call to test the paths with this vertex added. After the recursive
	 * call, we unvisit the vertex so we can "go back", and thus achieve backtracking.
	 * 
	 * @param u, integer representing the current vertex in the call
	 * @param visited, boolean array representing the visited vertices in the path, 
	 * true if visited, false otherwise, and since the vertices are represented  with 
	 * integers, we use the index to identify the vertices 
	 * @param path, list of integer where every element is a vertex in the graph
	 * @param cost, double representing the total cost of travel in a given cycle
	 * on the graph 
	 * @param resPath, list of integer used to record the traveled paths, the last/optimum
	 * path will always be the last one in it (of size numVertices)
	 * @return cost, double representing the total cost of travel in the resulting path of
	 * the graph
     */
    private double mineHelper(int u, boolean[] visited, List<Integer> path, double cost,
    		List<Integer> resPath) {
		visited[u] = true;
    	// if all the vertices are visited, then the Hamiltonian cycle exists
    	if (path.size() == this.numVertices) {
    		if (this.getPathCost(path) < cost) {
    			cost = this.getPathCost(path);
    			resPath.addAll(path);
    		}
    		return cost;
    	}

    	// Check if every edge starting from vertex `u` leads to a solution or not and if 
    	// that solution is small enough (min edge)
    	for (int i = 0; i < adjList.get(u).size(); i++) {
    		Edge minEdge = this.getMinConnection(adjList.get(u), visited);
    		Integer v = adjList.get(u).get(i).destinationVertex;
    		if (!visited[v] && v == minEdge.destinationVertex) {
    			path.add(v+1);
    			cost = mineHelper(v, visited, path, cost, resPath);
    			// backtrack for the path
    			visited[v] = false;  // so v could be used in another path
    			path.remove(path.size() - 1);
    		}
    	}
    	return cost;
    }
}
