#pragma once
#include <iostream>
#include "GraphExceptions.h"
#include "SafeArray/SafeArray.h"

using namespace std;

template <class T>
class Graph
{
 public:
    Graph();
    ~Graph();

    void addVertex(T& v);
    void setEdge(T& from, T& to, int weight);
    void topSort();
    int size();
    bool isEmpty();
    bool getEdge(const T& fromVertex, const T& toVertex);
    int getEdgeWeight(const T& fromVertex, const T& toVertex);
    SafeArray <T*> getIncomingNeighbors(const T& vertex);
    SafeArray <T*> getOutgoingNeighbors(const T& vertex);
    void print();
    void myPrint();
    void printRoutingTable();

 private:
    //methods
    bool find(const T& v);
    void resize();
    int findIndex(const T& v);
    void printStars();
 
    //data members
    
    //non-allocatable
    SafeArray <T> vertices;
    SafeArray <T*> sortedVertices;

    //number of vertices
    int numElements;
    //number of spaces allocated
    int capacity;

    //allocatable
    int** edges;
};

//ctor
template <class T>
Graph <T> :: Graph()
{
    //intialize numElements
    numElements = 0;

    //by default start with capacity of 10
    capacity = 10;

    //initialize edges
    //create an array of capacity int pointers
    edges = new int*[capacity];

    //for each row
    for(int i = 0; i < capacity; i++)
    {
	//create single dimension array of ints
	//store in corresponding int pointer
	edges[i] = new int[capacity];
    }

    //all rows and cols
    for(int i = 0; i < capacity; i++)
    {
	for(int j = 0; j < capacity; j++)
	{
	    //initialize all weights to be zero
	    edges[i][j] = 0;
	}
    }
}

//dtor
template <class T>
Graph <T> :: ~Graph()
{
    //make sure to delete allocatable members
    
    //delete edges
    //for each row
    for(int i = 0; i < capacity; i++)
    {
	//delete each row in the 2d array
	delete [] edges[i];
    }
}

//addVertex
template <class T>
void Graph <T> :: addVertex(T& v)
{
    //is element already a vertex?
    if(find(v))
    {
	cout << "Exception...element already in graph" << endl;
	//GraphException error;
	//throw error;
    }
    else
    {
	//update numElements
	numElements++;
	
	//set data in vertices    
	vertices.push_back(v);
	
	//update adjacency matrix if needed
	if(numElements > capacity)
	{
	    resize();
	}	
    }
}

//setEdge
template <class T>
void Graph <T> :: setEdge(T& from, T& to, int weight)
{
    //find index of from and to in vertices
    int fromIndex = findIndex(from);
    int toIndex = findIndex(to);

    //set edges[from][to] equal to weight
    edges[fromIndex][toIndex] = weight;
}

//topSort
template <class T>
void Graph <T> :: topSort()
{
    printStars();

    cout << "Finding a path through the graph, calling topSort" << endl;

    //first of all, clear previous sortedVertices array if not empty
    while(!(sortedVertices.isEmpty()))
    {
	sortedVertices.pop_back();
    }

    //make a copy of the edges matrix
    int** edgesCopy = new int*[numElements];

    //for each row
    for(int i = 0; i < numElements; i++)
    {
	edgesCopy[i] = new int[numElements];
    }	

    //for all rows and cols
    for(int i = 0; i <numElements; i++)
    {
	for(int j = 0; j < numElements; j++)
	{
	    //set edgesCopy[i][j] to edges[i][j]
	    edgesCopy[i][j] = edges[i][j];
	}
    }

    //use this array to keep track of where have been
    SafeArray <int> traversed;

    //want to repeat code numElements times
    for(int q = 0; q < numElements; q++)
    {
	bool hasIncoming;

	//find vertex with no incoming edges

	//want to find a column of zeros
	//for all columns
	for(int j = 0; j < numElements; j++)
	{
	    bool inTraversed = false;

	    //run through traversed to see if column j...
	    //has been added
	    for(int r = 0; r < traversed.size(); r++)
	    {
		if(traversed.at(r) == j)
		{
		    inTraversed = true;
		}
	    }
	    
	    //if this column hasn't been added to traversed yet
	    if(!inTraversed)
	    {
		hasIncoming = false;

		//for all rows
		for(int i = 0; i < numElements; i++)
		{
		    //if the entry is not zero, this destination...
		    //has a prerequisite destination
		    if(edgesCopy[i][j] != 0)
		    {
			hasIncoming = true;
			//have found out that this node has incoming
			//no need to continue here...
			break;
		    }
		}
	    
		//was node found to have incoming?
		if(!hasIncoming)
		{
		    //process node
		    cout << vertices.at(j) << endl;

		    //remove outgoing edges
		    for(int i = 0; i < numElements; i++)
		    {
			edgesCopy[j][i] = 0;
		    }

		    //add this column to the traversed list
		    traversed.push_back(j);
		    //add this column to the list of sorted vertices
		    sortedVertices.push_back(&(vertices.at(j)));

		    //want to restart loop as adjacency matrix...
		    //has been modified
		    break;
		}
	    }
	}
    }

    //make sure to delete edgesCopy here

    //for each row
    for(int i = 0; i < numElements; i++)
    {
	//delete each row
	delete [] edgesCopy[i];
    }

    printStars();
}

//size
template <class T>
int Graph <T> :: size()
{
    return numElements;
}

//isEmpty
template <class T>
bool Graph <T> :: isEmpty()
{
    return (numElements == 0);
}

//getEdge
template <class T>
bool Graph <T> :: getEdge(const T& fromVertex, const T& toVertex)
{
    bool retVal = false;

    int fromIndex = findIndex(fromVertex);
    int toIndex = findIndex(toVertex);

    if(edges[fromIndex][toIndex] != 0)
    {
	retVal = true;
    }

    return retVal;

}

//getEdgeWeight
template <class T>
int Graph <T> :: getEdgeWeight(const T& fromVertex, const T& toVertex)
{
    int weight;

    int fromIndex = findIndex(fromVertex);
    int toIndex = findIndex(toVertex);

    weight = edges[fromIndex][toIndex];

    return weight;    
}

//getIncomingNeighbors
template <class T>
SafeArray <T*> Graph <T> :: getIncomingNeighbors(const T& vertex)
{
    //initialize return array
    SafeArray <T*> retArray;

    //find index of vertex
    int index = findIndex(vertex);

    //want to find all rows with nonzero values...
    //in the column corresponding to index 

    //for all rows
    for(int i = 0; i < numElements; i++)
    {
	if(edges[i][index] != 0)
	{
	    //get address of incoming vertex...
	    //use to create pointer
	    T* elemAddress = &(vertices.at(i));
	    //add pointer to retArray
	    retArray.push_back(elemAddress);
	}
    }

    return retArray;
}

//getOutgoingNeighbors
template <class T>
SafeArray <T*> Graph <T> :: getOutgoingNeighbors(const T& vertex)
{
    //initialize return array
    SafeArray <T*> retArray;

    //find index of vertex
    int index = findIndex(vertex);

    //want to find all columns with nonzero values...
    //in the row corresponding to index 

    //for all cols
    for(int j = 0; j < numElements; j++)
    {
	if(edges[index][j] != 0)
	{
	    //get address of incoming vertex...
	    //use to create pointer
	    T* elemAddress = &(vertices.at(j));
	    //add pointer to retArray
	    retArray.push_back(elemAddress);
	}
    }

    return retArray;

}

//find
template <class T>
bool Graph <T> :: find(const T& v)
{
    bool retVal = false;

    for(int i = 0; i < vertices.size(); i++)
    {
	if(vertices[i] == v)
	{
	    retVal = true;
	    break;
	}
    }

    return retVal;

}

//findIndex
template <class T>
int Graph <T> :: findIndex(const T& v)
{
    //if retVal is ever -1, we know the vertex wasn't found
    int retVal = -1;
    
    //search vertices for the element
    for(int i = 0; i < numElements; i++)
    {
	//if found element, set retVal to current i
	if(vertices.at(i) == v)
	{
	    retVal = i;
	    break;
	}
    }

    return retVal;

}

//resize
template <class T>
void Graph <T> :: resize()
{
    //want to double capacity in this function...

    //temp Matrix for copying purposes
    int** temp;
    //new capacity of matrix
    int newCap = capacity*2;

    //initialize tempMatrix
    temp = new int*[newCap];

    //for each row...
    for(int i = 0; i < newCap; i++)
    {
	temp[i] = new int[newCap];
    }

    //all rows and cols
    for(int i = 0; i < newCap; i++)
    {
	for(int j = 0; j < newCap; j++)
	{
	    //initialize all weights to be zero
	    temp[i][j] = 0;
	}
    }

    //now, copy old info into temp
    for(int i = 0; i < capacity; i++)
    {
	for(int j = 0; j < capacity; j++)
	{
	    temp[i][j] = edges[i][j];	    
	}
    }

    //the deed is done, handle memories
    
    //toDelete points to old edges
    int** toDelete = edges;

    //edges points to new edges
    edges = temp;

    //delete oldEdges
    //for each row
    for(int i = 0; i < capacity; i++)
    {
	//delete each row in the 2d array
	delete [] toDelete[i];
    }

    //update capacity
    capacity = newCap;

}

//print
template <class T>
void Graph <T> :: print()
{
    //print each vertex and all of...
    //its incoming/outgoing neighbors

    //for each vertex in vertices
    for(int p = 0; p < vertices.size(); p++)
    {
	printStars();

	cout << vertices.at(p) << endl;

	cout << "Incoming neighbors: " << endl;

	const SafeArray <T*> & incoming = getIncomingNeighbors(vertices.at(p));

	//print incoming neighbors

	for(int q = 0; q < incoming.size(); q++)
	{
	    for(int r = 0; r < 10; r++)
	    {
		cout << " ";
	    }

	    cout << *(incoming.at(q)) << endl;
	}				 

	cout << "Outgoing neighbors: " << endl;

	const SafeArray <T*> & outgoing = getOutgoingNeighbors(vertices.at(p));

	//print outgoing neighbors

	for(int q = 0; q < outgoing.size(); q++)
	{
	    for(int r = 0; r < 10; r++)
	    {
		cout << " ";
	    }

	    cout << *(outgoing.at(q)) << endl;
	}

	printStars();
    }

}

//myPrint
template <class T>
void Graph <T> :: myPrint()
{
    //print vertices (calling SafeArray.print())
    vertices.print();

    //all rows and cols
    for(int i = 0; i < capacity; i++)
    {
	for(int j = 0; j < capacity; j++)
	{
	    cout << edges[i][j] << " ";

	}	
	
	cout << endl;

    }

}

//printStars
template <class T>
void Graph <T> :: printStars()
{
    for(int i = 0; i < 30; i++)
    {
	cout << "*";
    }

    cout << endl;
}

//printRoutingTable
template <class T>
void Graph <T> :: printRoutingTable()
{
    //topologically sort array of the vertices
    topSort();

    //topSort refreshes sortedVertices
    //sortedVertices.print();

    //for each node, treating each node as source
    //and incrementing source each iteration...
    for(int source = 0; source < sortedVertices.size(); source++)
    {
	//create a SafeArray for this particular...
	//node's predecessor nodes
	T* predecessors[numElements];

	//initialize array of predessors
	for(int i = 0; i < numElements; i++)
	{
	    predecessors[0] = 0;
	}	

	//initialize a dValues SafeArray
	//this array holds the distance from the current source...
	//to the node in question
	SafeArray <int> dValues(numElements);

	//initialize all but source node dValues to be infinite
	for(int r = 0; r < dValues.size(); r++)
	{
	    if(r == source)
		dValues.at(r) = 0;
	    else
		dValues.at(r) = 1000000;
	}

	//dValues.print();

	//starting from source and continuing to end of sortedGraph
	//relax all of the edges
	for(int currentNode = source; currentNode < sortedVertices.size(); currentNode++)
	{
	    //get outgoing neighbors
	    const SafeArray <T*> & outgoing = getOutgoingNeighbors(*(sortedVertices.at(currentNode)));

	    //for all of the outgoing neighbors, relax edges
	    for(int q = 0; q < outgoing.size(); q++)
	    {
		//need to find index of outgoing in sorted array
		int outgoingIndex = -1;
		
		for(int s = 0; s < sortedVertices.size(); s++)
		{
		    if(sortedVertices.at(s) == outgoing.at(q))
		    {
			outgoingIndex = s;
			break;
		    }
		}

		//find the original index in edges matrix of
		//currentNode and outgoingNode
		int currentOriginalIndex = findIndex(*(sortedVertices.at(currentNode)));
		int outgoingOriginalIndex = findIndex(*(sortedVertices.at(outgoingIndex)));

		//if distance from source to outgoing is larger than
		//distance from source to currentNode plus
		//distance from currentNode to outgoing
		if(dValues.at(outgoingIndex) > dValues.at(currentNode) + 
		   edges[currentOriginalIndex][outgoingOriginalIndex])
		{
		    //set dValue at outgoingIndex to smaller distance
		    dValues.at(outgoingIndex) = dValues.at(currentNode) + 
			edges[currentOriginalIndex][outgoingOriginalIndex];
		    
		    //set the predecessor of outgoing to the currentNode
		    predecessors[outgoingIndex] = sortedVertices.at(currentNode);
		}
	    }
	}

	//for this node as source, print out predecessors
	cout << "SOURCE NODE: " << *(sortedVertices.at(source)) << " to: " << endl;
			
	for(int i = 0; i < numElements; i++)
	{
	    if(dValues.at(i) < 1000000 &&
	       dValues.at(i) != 0)
	    {
	    	cout << *(sortedVertices.at(i)) << ": " << 
		    " total weight: " << dValues.at(i) << " ";

		int searchInt = i;
		T* predecessor;

		//trace predecessors up to source node
		while(searchInt != source)
		{
		    predecessor = predecessors[searchInt];

		    //cout << *(predecessors[searchInt]) << endl; 
		    
		    if(predecessor == sortedVertices.at(source))
			break;

		    //find the index of predecessor in the sorted array
		    for(int j = 0; j < sortedVertices.size(); j++)
		    {
			if(sortedVertices.at(j) == predecessor)
			{
			    searchInt = j;
			    break;
			}
		    }
		}

		cout << "Next hop: " << *(sortedVertices.at(searchInt)) << endl;
	    }
	}

	cout << endl << endl;
	
    }

}
