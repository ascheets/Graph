#include <iostream>
#include <string>
#include "Graph.h"

using namespace std;

int main()
{
	//graph that holds locations (strings)
	Graph < string > graph;

	//locations
	string chicago = "Chicago, IL";
	string newYork = "New York, NY";
	string seattle = "Seattle, WA";
	string sanFrancisco = "San Francisco, CA";
	string tempe = "Tempe, AZ";
	string boulder = "Boulder, CO";
	string stLouis = "St. Louis, MO";
	string minneapolis = "Minneapolis, MN";
	string orlando = "Orlando, FL";
	string newHaven = "New Haven, CT";

	//add the locations to the graph
	graph.addVertex(chicago);
	graph.addVertex(newYork);
	graph.addVertex(seattle);
	graph.addVertex(sanFrancisco);
	graph.addVertex(tempe);
	graph.addVertex(boulder);
	graph.addVertex(stLouis);
	graph.addVertex(minneapolis);
	graph.addVertex(orlando);
	graph.addVertex(newHaven);

	//create routes between the cities
	//chicago destinations
	graph.setEdge(chicago, minneapolis, 354);
	graph.setEdge(chicago, seattle, 1736);
	graph.setEdge(chicago, boulder, 927);
	graph.setEdge(chicago, stLouis, 259);
	graph.setEdge(chicago, orlando, 984);

	//east coast routes
	graph.setEdge(minneapolis, newHaven, 1056);
	graph.setEdge(newHaven, newYork, 67);
	graph.setEdge(orlando, newYork, 940);

	//west coast routes
	graph.setEdge(seattle, sanFrancisco, 681);
	graph.setEdge(boulder, sanFrancisco, 938);
	graph.setEdge(stLouis, tempe, 1264);
	graph.setEdge(tempe, sanFrancisco, 663);

	//print all the vertex info
	graph.print();

	//find a path through the graph
	graph.topSort();

	// //find the routes through the graph
	// graph.printRoutingTable();

	return 0;
}
