/*
* AUTHOR: Kevin Nisterenko
* FILE: PA11Main.java
* ASSIGNMENT: Programming Assignment 11 - The Traveling Salesman Problem
* COURSE: CSc 210; Fall 2021
* PURPOSE: This program defines the PA11Main class, which takes a .mtx file as input 
* builds the graph represented by it by calling the appropriate methods from DGraph and 
* then takes another command line input to call the appropriate solution to the TSP 
* problem on the built graph. The appropriate message is printed on regular output
* for each input, showcasing the cost of the path, and depending on the command, the
* path itself or the time it took in milliseconds for the solution to be found. 
*
* This file takes command line inputs: 
* -The first input is a string representing a .mtx file path which will be used to 
*  construct the graph. 
* 
* -The second input is the command to be executed. It can be one of the following:
*   HEURISTIC - an approximate algorithm that uses a heuristic,
*   BACKTRACK - an exact algorithm that uses brute-force search,
*   MINE - a variation on the heuristic and backtracking approaches that
*          improves performance,
*   TIME - a measuring of the running time of the above approaches
*/

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;

public class PA11Main {
	public static void main(String[] args) {
		// Open file
		Scanner fileObj = null;
	    try {fileObj = new Scanner(new File(args[0]));} 
	    catch (FileNotFoundException e) {e.printStackTrace();}
	    
	    DGraph graph = readFile(fileObj);
	    
	    runCommand(args[1], graph);
	 }
	
	/* 
	 * Public method used to build the graph from the input file, it uses a while loop to 
	 * go over each line and initiate a new DGraph object and also build the edges between 
	 * vertices. The graph is then returned by this method.
	 * 
	 * @param fileObj, scanner object representing the input file to build the graph from
	 * @return graph, DGraph object, representing a directed and weighted graph 
	 */
	public static DGraph readFile(Scanner fileObj) {
		boolean flag = true;
		DGraph graph = null;
		while (fileObj.hasNext()) {
	    	  String line = fileObj.nextLine();
	    	  // Skip the header/comment lines on the file
	    	  if (line.charAt(0) != '%') {
	    		  String[] lineArr = line.trim().split("( )+");
    			  // Flag the first line to build empty (no connections) 
    			  // graph of appropriate size
	    		  if (flag) {
	    			  graph = new DGraph(Integer.valueOf(lineArr[0]));
	    			  flag = false;
	    		  // Add edges to 'build' the connections to the graph
	    		  } else {
	    			  int from = Integer.valueOf(lineArr[0]) - 1;
	    			  int to = Integer.valueOf(lineArr[1]) - 1;
	    			  double weight = Double.valueOf(lineArr[2]);
	    			  graph.addEdge(from, to, weight);
	    		  }
	    	  }
	      }
		return graph;
	}
	
	/* 
	 * Public method used to call the appropriate graph method or function (if input 
	 * is 'time') on the given graph. 
	 * 
	 * @param command, string representing the command specified in the command line for which 
	 * TSP operation to perform.
	 * @param graph, DGraph object, representing a directed and weighted graph 
	 */
	public static void runCommand(String command, DGraph graph) {
		if (command.toLowerCase().equals("heuristic")) {
			graph.heuristicTSP(0, true);
		} else if (command.toLowerCase().equals("backtrack")) {
			graph.backtrackTSP(0, true);
		} else if (command.toLowerCase().equals("mine")) {
			graph.mineTSP(0, true);
		} else if (command.toLowerCase().equals("time")) {
			timeCommand(graph);
		}
	}
	
	/* 
	 * Public method used when the specified input command is 'time'. It calls each of the 
	 * three TSP traversals and times them in milliseconds. This information is printed in 
	 * the output.
	 * 
	 * @param graph, DGraph object, representing a directed and weighted graph 
	 */
	public static void timeCommand(DGraph graph) {
		long startH = System.currentTimeMillis();
		double costH = graph.heuristicTSP(0, false);
		long endH = System.currentTimeMillis();
		
		long elapsedH = endH - startH;
		System.out.println("heuristic: cost = " + costH + ", " + elapsedH + " milliseconds");
		
		long startM = System.currentTimeMillis();
		double costM = graph.mineTSP(0, false);
		long endM = System.currentTimeMillis();
		long elapsedM = endM - startM;
		System.out.println("mine: cost = " + costM + ", " + elapsedM + " milliseconds");
		
		long startB = System.currentTimeMillis();
		double costB = graph.backtrackTSP(0, false);
		long endB = System.currentTimeMillis();
		long elapsedB = endB - startB;
		System.out.println("backtrack: cost = " + costB + ", " + elapsedB + " milliseconds");
	}
}

