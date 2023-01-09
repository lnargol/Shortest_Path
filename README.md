---
title: "Calculating Shortest Paths in Graphs"
subtitle : "CSE 701 - Final Project"
author: "Nargol Lotfizadeh - 400486891"
date: "2022-12"
output: pdf_document
---


# Documentation

# Summary

* In discrete mathematics, a graph is a structure amounting to a set of objects in which some pairs of the objects are in some sense ”related”. In graph theory the shortest path problem is the problem of finding a path between two vertices (or nodes) in a graph such that the sum of the weights of its constituent edges is minimized.Two vertices are adjacent when they are both incident to a common edge. A path in an undirected graph is a sequence of vertices
$P = v_{1}, v_{2}, v_{3}, ..., v_{n} \in V \times V \times ... \times V$
such that $v_{i}$ is adjacent to $v_{i+1}$ for $1\leq i < n$. Such a path $P$ is called a path of length $n-1$ from $v_{1} to v_{n}$.

* Let $e_{i,j}$ be the edge incident to both $v_{i}$ and $v_{j}$. Given a real-valued weight function $f : E \longmapsto R$ and a graph G, the shortest path from $v$ to $v'$ is the path $P = (v_{1}, v_{2},...,v_{n})$ (where $v_{1} = v$ and $v_{n} = v'$)that over all possible n minimizes the sum $\sum_{i = 1}^{n-1} f(e_{i,i+1}).$

# Graph Storing Format
 * In order to store our graph, we need to store its vertices and edges, and in order to do so, we are going to define classes for edges and graphs. we store out graphs in an adjacency list, which will be explained in detail below. We will use a vector of vectors from class instances of edge, *vector< vector>edge >*.


# Classes and Functions


* Class `edge`: 
  This class consists of 3 member attributes; start, end, and weight. The first is the node that the edge originates from, the middle is the destination of the node, and the last is the weight of the edge, which will be used in calculating the cost of paths between nodes.
Getters and Setters have been written for this class for the attributes to be accessible.


* Class `graph`: This class is the main class to store our graph. It has four attributes. Setters and Getters have been generated for the private attributes. 
    + `uint64_t vertices`: This is an integer type variable and stores the number of vertices present in the graph.
    + `vector<<vector>edge> adj`: This is the adjacency list that represents the whole graph. For each vertex, we store a vector of edges which are originated from that vertex, and then put them all together in a vector.
    + `vector<int64_t> weights`: All weights used in the graph will be added to this vector at the time of construction. We will then use this vector to assess whether there is a negative value among weights or not. 
    + `vector<edge> edges`: All edges will be added to this vector and later used in the algorithms.

    This class also has functions. These functions are:
    + `void add_edge` : It takes the start, end, and the weight of an edge and constructs it. Then, it adds the new edge to the vector of edges for the graph. 
    + `bool non_negative_weights`: As we will see further in this file, the Dijkstra algorithm does not accept negative weights in the graph. As such, this function will be used to detect if the graph contains negative weights. 
    This class has only one constructor, which generates the graph by reading a file. The file should have a specific format. 
    + `graph(string file_name)`: This constructor takes the file name and constructs the graph.The first line consists of two integers, namely the number of vertices and the number of edges. Using `for`, we then proceed to read each line, which consists of three elements. The first element is the starting node, the second element is the destination node, and the third element is the weight. After reading each line and separating the elements using white space, the program constructs an instance of the type `edge`, adds it to the vector of edges corresponding with the starting node. For example, if the program reads the line `2 3 1`, it constructs the new edge with starting node 2, ending node 3, weight 1, and adds it to the vector of edges for vertex 2.
    

* `Class timer`: This class basically is used to measure the amount of time each algorithm takes to be executed. This class uses the `chrono` library. 

* `cast`: This function has been implemented to safely store unsigned integers into signed ones and vice versa. This is only viable in case the integer in non negative, of course.

# Algorithms

## Bellman Ford

* Bellman Ford is a single source algorithm that computes the shortest path from a specified vertex or node, to all other vertices. This algorithm accepts negative weights in the graph, and is able to detect a negative cycle, if there exists one. In this case, no shortest path can be found, as moving on the negative cycle can be done unlimited times, and every time we can find a shorter path. Bellman Ford can detect this, and then it throws an exception, stating that a negative cycle exists. 

* Bellman Ford advantages from a concept named **relaxation**, which basically means substituting better (cheaper) paths for the previous (more expensive) ones. The algorithm starts by defining distance and predecessor vector of integers, which hold the distance of the source vertex from all other vertices. For example, the integer in the second index of the distance vector shows the current shortest path from source to vertex two. In initialization part, all distances are set to infinity, and all predecessors are set to -1, which indicates that no predecessor has been found yet. The Bellman Ford algorithms simply *relaxes* all edges, and does this process `vertices - 1` times. In each execution of the relaxation, the distance between more nodes get calculated correctly until eventually all the distances are the shortest distances possible. It basically updates the current shortest paths found in each step to see if the costs of the paths can be reduced.

* A it is a single source algorithm, we use a for loop to calculate all distances for all vertices in the graph.



## Dijkstra

* Similar to Bellman Ford, Dijkstra is also a single source shortest path, and probably the most famous one. The initialization part of this algorithm is similar to that of Bellman Ford, and it also uses relaxation and updates the found distances in each step. However, instead of going over edges `vertices - 1` times (which is Bellman Ford's approach), Dijkstra select the closest vertex that has still yet to be visited, and performs relaxation of all of that vertex's outgoing edges. Each vertex is removed from the unvisited vector once it has been visited. 

* If the user only wants the shortest path between two specific nodes, both Dijkstra and Bellman Ford can easily be modified in a way that they stop once the algorithm found the shortest path from source to sink (the target node).

* In more advanced version of this algorithm, a min heap is used as a priority queue for finding the vertex with minimum distance and this enhances the algorithm's performance. 

<p align="center">
    <img src="dijkstra.png" alt="drawing" width="200" height="600"/>
</p>
## Floyd Warshall

* Unlike Bellman Ford and Dijkstra, Floyd Warshall is not a single source shortest path algorithm and is designed to find all shortest paths between all the vertices. The algorithm does this by comparing all possible paths through the graph between each pair of vertices. Floyd Warshall is an example of dynamic programming.

* The principle of this algorithm is as follows:

    + For each of the pairs between vertices, the shortest path from $i$ to $j$ using only the first $k$ (let's show this as $shortestPath(i,j,k)$) vertices as the intermediate nodes along the path, two scenarios exist: 
    
      First, the shortest path does not go through k (and therefore only uses vertices in the set of ${1,...,k-1}$).
    
      Second, the shortest path does go through k (from $i$ to $k$ and then from $k$ to $j$, both only using intermediate vertices in ${1,...,k-1}$).
    
      The above logic leads us to the following formula which helps us calculates all shortest paths in the graph:
    
      $shortestPath(i,j,k) = min(shortestPath(i,j,k-1), shortestPath(i,k,k-1) + shortestPath(k,j,k-1))$.

# Algorithm Comparison

* First of all, Dijkstra algorithm can not deal with negative weights, which is not the case for Bellman Ford.

*Bellman Ford and Dijkstra are both single source shortest path algorithms, while Floyd Warshall is not. 

* Dynamic programming approach is taken to implement the Bellman Ford Algorithm while Greedy approach is followed for Dijkstra.

* Floyd Warshall can be implemented in a distributed system which makes it suitable for data structures such as graphs of graphs. These kind of data structures can be used in maps.

# Operator Overloading

* Multiple overloads have been included in the code. As we have a self-defined format of storing the graph, in many situations it is useful to print out the adjacency list, the list of vectors (unsigned integers), the list of weights (signed integers), and so on. As such, five overloads have been written for printing vectors to facilitate debugging and observing the flow of the program.  



# Flow of the Program

* The command line can consist of either 3 or 4 arguments. The proper format is 
    + `[EXE File Name][Input File Name][Output File Name][(Optional) Desired Algorithm Number]`
    + If no argument is mentioned for the optional desired algorithm number, the program executes all three algorithms and output the result.

* A separate timer has been put for measuring the time each algorithm takes, and its then printed out to both the terminal and the file for the user to compare.


# Error Handling 

* Several Error handling has been embedded in the program. First of all, the format of the command line inputs is checked in *main*, whether there are the right number of arguments specified in the command line, if the fourth argument is really a number, and whether the files have been opened correctly. In each case if there is an error, the program prints out the error and terminates the program.

* Another issue to have in mind is that Dijkstra Algorithms does not accepts graphs with negative weights. Hence, this matter is checked in the Dijkstra function to ensure the weights in the graph are completely non-negative. If this is not the case, the function throws an exception, which would be caught in main.



# Input and Output

* As explained the graph is given to the program via an input file. This file should consist of 2 numbers, let's say m (number of vertices), and n (number of edges) respectively. It should be followed by n lines, each line corresponding to an edge. The input file name is given to the program as the second argument of the command line.

* Output file name is the third argument of the command line specified by the user. All algorithms' outputs will be written to the file, alongside with their execution time. The execution time, and the completion of the algorithms will also be printed out in terminal.

* Sample compilation command line : `g++ -Wall -Wextra -Wconversion -Wsign-conversion -Wshadow -Wpedantic -std=c++20 shortest_path.cpp -o shortest_path`

* Sample execution : `./shortest_path graph.txt out.txt 2`   



* **Sample Input File**:
    +``` 4 3
    
      0 2 2
      
      1 3 1
      
      2 1 1
      
      1 0 3```
* **Sample Command Line**: [./shortest_path] [graph.txt] [output.txt] 1
* **Sample Output File**:
    
   ``` ** OUTPUT FOR THE GIVEN GRAPH IS AS FOLLOWS: ** 
    
    ----------------------------------------------------
    BELLMAN FORD Algorithm has calculated the shortest paths .
    The paths and their corresponding cost are as follows : 
    
    Shortest path found from source 0 to sink 1 with cost 3 --> (0, 2, 1)
    
    Shortest path found from source 0 to sink 2 with cost 2 --> (0, 2)
    
    Shortest path found from source 0 to sink 3 with cost 4 --> (0, 2, 1, 3)
    
    Shortest path found from source 1 to sink 0 with cost 3 --> (1, 0)
    
    Shortest path found from source 1 to sink 2 with cost 5 --> (1, 0, 2)
    
    Shortest path found from source 1 to sink 3 with cost 1 --> (1, 3)
    
    Shortest path found from source 2 to sink 0 with cost 4 --> (2, 1, 0)
    
    Shortest path found from source 2 to sink 1 with cost 1 --> (2, 1)
    
    Shortest path found from source 2 to sink 3 with cost 2 --> (2, 1, 3)
    
    No path from 3 to 0 has been found!
    
    No path from 3 to 1 has been found!
    
    No path from 3 to 2 has been found!
    ----------------------------------------------------
    [BELLMAN FORD] Elapsed time: 0.174833 milli seconds.
    ```

Thank you for reading! :)