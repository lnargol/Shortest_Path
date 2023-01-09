/**
 * @file shortest_path.cpp
 * @author Nargol Lotfizadeh (lotfizan@mcmaster.ca)
 * @brief Final Project for CSE 701 - McMaster University
 *        This program is designed to implement three algorithms proposed for finding shortest paths in graphs.
 * @version 0.1
 * @date 2022-12-27
 *
 * @copyright Copyright (c) 2022
 *
 */

#include <iostream>
#include <vector>
#include <fstream>
#include <limits>
#include <algorithm>
#include <chrono>

using namespace std;

/**
 * @brief This class has been included to measure execution time of functions.
 * The code for this class has been borrow from "https://baraksh.com/CSE701/".
 */

class timer
{
public:
    /**
     * @brief Stores the start time of the timer
     *
     */
    void start()
    {
        start_time = chrono::steady_clock::now();
    }

    /**
     * @brief Stores the end time of the timer
     *
     */
    void end()
    {
        elapsed_time = chrono::steady_clock::now() - start_time;
    }

    /**
     * @brief Calculates the time passed by seconds
     *
     * @return double
     */
    double seconds() const
    {
        return elapsed_time.count();
    }
    // chrono::duration<double, milli> elapsed_time_milli = end_time - start_time;
private:
    chrono::time_point<chrono::steady_clock> start_time = chrono::steady_clock::now(); /** < Stores the current time as the start */
    chrono::duration<double, milli> elapsed_time = chrono::duration<double>::zero();   /** <Calculates the time passed by milliseconds */
};

/**
 * @brief This class is used to generate edges in order to construct graphs later.
 *
 */
class edge
{
public:
    /**
     * @brief Construct a new edge object
     *
     * @param u starting node
     * @param v destination node or ending node
     * @param q the weight associated with the edge
     */
    edge(uint64_t u, uint64_t v, int64_t q)
    {
        start = u;  /**< Stores the origin of the edge */
        end = v;    /**< Stores the destination of the edge */
        weight = q; /**< Stores the weight (cost) of the edge */
    }

    /**
     * @brief Get the start object
     *
     * @return starting node
     */
    uint64_t get_start()
    {
        return start;
    }

    /**
     * @brief Get the end object
     *
     * @return ending node
     */
    uint64_t get_end()
    {
        return end;
    }

    /**
     * @brief Get the weight object
     *
     * @return weight of the edge
     */
    int64_t get_weight()
    {
        return weight;
    }

private:
    uint64_t start; /**< The origin of the edge */
    uint64_t end;   /**< The destination of the edge */
    int64_t weight; /**< The weight of the edge (cost) */
};

/**
 * @brief operator << overloaded for a vector of type uint64_t
 *
 * @param out the outstream
 * @param v the vector
 * @return ostream& printed output
 */
ostream &operator<<(ostream &out, const vector<uint64_t> &v)
{
    out << '(';
    for (uint64_t i = 0; i < v.size() - 1; i++)
        out << v[i] << ", ";
    out << v[v.size() - 1] << ')';
    return out;
}

/**
 * @brief operator << overloaded for a vector of type int64_t (used for printing predecessor)
 *
 * @param out  the outstream
 * @param v the vector
 * @return ostream&  printed output
 */
ostream &operator<<(ostream &out, const vector<int64_t> &v)
{
    out << '(';
    for (uint64_t i = 0; i < v.size() - 1; i++)
        out << v[i] << ", ";
    out << v[v.size() - 1] << ')';
    return out;
}

/**
 * @brief operator << overloaded for a vector of type edge
 *
 * @param out the outstream
 * @param v the vector
 * @return ostream&
 */
ostream &operator<<(ostream &out, vector<edge> &v)
{
    for (uint64_t i = 0; i < v.size(); i++)
    {
        out << v.at(i).get_start() << " " << v.at(i).get_end() << " " << v.at(i).get_weight() << endl;
    }
    // out << ')';
    return out;
}

/**
 * @brief operator << overloaded for a vector of vectors of type edge (an example is adjacency list)
 *
 * @param out the outstream
 * @param v the vector
 * @return ostream&
 */
ostream &operator<<(ostream &out, vector<vector<edge>> &v)
{
    out << "{\n";
    for (uint64_t i = 0; i < v.size(); i++)
    {
        out << v.at(i);
    }
    out << '}';
    return out;
}

/**
 * @brief operator << overloaded for an edge
 *
 * @param out the outstream
 * @param e the vector
 * @return ostream&
 */
ostream &operator<<(ostream &out, edge &e)
{
    out << '(' << e.get_start() << " " << e.get_end() << " " << e.get_weight() << ')' << " \n";
    return out;
}

/**
 * @brief During the course of the program we store signed integers in an unsigned integer.
 *        Of course this is only valid when the integer is non negative.
 *        This function helps not to get implicit conversion warnings.
 *
 * @param i The signed integer
 * @return uint64_t The value of the signed integer
 */
uint64_t cast(int64_t &i)
{
    if (i < 0)
        throw invalid_argument("Number is negative! Can't cast it to uint64_t!");
    else
        return static_cast<uint64_t>(i);
}

/**
 * @brief During the course of the program we store unsigned integers in a signed integer.
 *        This function helps not to get implicit conversion warnings.
 *
 * @param i The unsigned integer
 * @return The signed integer containing the unsigned integer
 */
int64_t cast(uint64_t &i)
{
    return static_cast<int64_t>(i);
}

/*! @class graph
    @brief handles generating the graph using reading input from a file

*/
class graph
{
public:
    /**
     * @brief Construct a new graph object
     *
     * @param file_name  input file name
     */
    graph(string file_name)
    {
        /** The input file is being read line by line,
         * Splitted by white space,
         * and converted into instances of class edge.
         */

        string line;
        ifstream my_file;
        // bool read_first_line = false;
        string delimiter = " ";
        size_t pos = 0;
        string token;
        uint64_t vertices_no = 0;
        uint64_t num_of_edges = 0;
        vector<edge> all_edges;

        my_file.open(file_name);

        if (!my_file.is_open())
            throw invalid_argument("\n\nError in opening file!\n\n");

        getline(my_file, line); // read first line for vertices and number of edges
        if ((pos = line.find(delimiter)) != std::string::npos)
        {
            vertices_no = stoull(line.substr(0, pos));
            line.erase(0, pos + delimiter.length());
            num_of_edges = stoull(line);
            // cout << "num of edges : " << num_of_edges << endl;
        }
        vector<vector<edge>> adj_vec(vertices_no);
        for (uint64_t i = 0; i < num_of_edges; i++)
        {
            // cout << "i is " << i << endl;
            getline(my_file, line);
            vector<uint64_t> edge_elements;
            while (((pos = line.find(delimiter)) != string::npos) || (pos = line.find("\n") != string::npos))
            {
                token = line.substr(0, pos);
                line.erase(0, pos + delimiter.length());
                edge_elements.push_back(stoull(token));
            }
            edge_elements.push_back(stoull(line));
            if (edge_elements.size() < 3)
            {
                throw invalid_argument("Invalid arguments entered for edges! ");
            }

            edge new_edge = edge(edge_elements.at(0), edge_elements.at(1), cast(edge_elements.at(2)));
            all_edges.push_back(new_edge);
            adj_vec.at(edge_elements.at(0)).push_back(new_edge);
        }
        // Constructed elements are assigned to the graph.
        vertices = vertices_no;
        adj = adj_vec;
        edges = all_edges;
    }

    /**
     * @brief function for detecting negative values among the weights
     *
     * @return true if there exists a negative weight
     * @return false if all weights are non negative
     */
    bool non_negative_weights() const
    {
        for (int64_t w : weights)
        {
            if (w < 0)
                return false;
            else
                continue;
        }
        return true;
    }

    /**
     * @brief Get the adj object
     *
     * @return vector<vector<edge> > adjacency list
     */
    vector<vector<edge>> get_adj() const
    {
        return adj;
    }

    /**
     * @brief Get the edges object
     *
     * @return vector<edge>  edges
     */
    vector<edge> get_edges() const
    {
        return edges;
    }

    /**
     * @brief Get the vertices object
     *
     * @return uint64_t  vertices
     */
    uint64_t get_vertices() const
    {
        return vertices;
    }

private:
    uint64_t vertices; /**< Detailed description after the member */
    vector<vector<edge>> adj;
    vector<int64_t> weights;
    vector<edge> edges;
    /**
     * @brief adding edge to the graph
     *
     * @param u the starting node
     * @param v the destination node
     * @param weight the weight of the edge
     */
    void add_edge(uint64_t u, uint64_t v, int64_t weight)
    {
        adj.at(u).push_back(edge(u, v, weight));
        weights.push_back(weight);
    }
};

/**
 * @brief function to discover whether a path exists according to the path that has been constructed for it
 *
 * @param vec the vector
 * @return true if a path exists
 * @return false if a path does not exist
 */
bool found_path(const vector<uint64_t> &vec)
{
    if (vec.size() == 1) // No path were found
        return false;

    else
        return true;
}

/**
 * @brief Using the vector of vectors of type int64_t, which is the distances and the predecessors of the vertices, it constructs the path from the source to sink, if there exists one.
 *
 * @param vec the vector constructed by either Dijkstra or Bellman Ford
 * @param source the source vertex
 * @param sink the sink vertex
 * @return vector <uint64_t> contains the elements pf the path in order
 */
vector<uint64_t> find_path(const vector<vector<int64_t>> &vec, const uint64_t &source, const uint64_t &sink)
{
    // cout << "source and sink are " << source << " " << sink << endl;

    vector<int64_t> dist = vec.at(0);
    vector<int64_t> pred = vec.at(1);

    uint64_t iterator = sink;
    vector<uint64_t> path;
    // bool path_found = true;
    while (iterator != source)
    {
        // cout << "iterator " << iterator << endl;
        if (pred.at(iterator) == -1)
        {
            // path_found = false;
            //  cout << "PATH NOT FOUND! \n";
            break;
        }

        path.push_back(iterator);
        // Follow the path one by one by getting their predecessor
        iterator = cast(pred.at(iterator));
    }

    path.push_back(source);

    reverse(path.begin(), path.end());
    // cout << path << endl;
    return path;
}

/**
 * @brief This algorithm is a single source shortest path algorithm for graphs. Detailed Description is available in Readme.
 *        Bellman Ford reaches the shortest path via relaxing all edges (vertices - 1) times.
 *         Relaxing is the procedure in which the distances are replaced with better (cheaper) paths by considering other edges and their weights.
 *
 * @param g the graph
 * @param source the source node
 * @return vector < vector < uint64_t> > This vector contains of 2 vectors of type uint64_t. The first vector (index 0) is the vector that has the distances for the source to all other node,
 * and the second vector (index 1) contains the predecessors of the vertices.
 */
vector<vector<int64_t>> bellman_ford(const graph &g, const uint64_t &source)
{

    vector<edge> edges = g.get_edges();
    uint64_t vertices = g.get_vertices();
    vector<int64_t> dist(vertices);
    vector<int64_t> pred(vertices);

    // Stage 1: Initialization
    for (uint64_t i = 0; i < vertices; i++)
    {
        dist.at(i) = numeric_limits<int>::max();
        pred.at(i) = -1;
    }
    dist.at(source) = 0;

    uint64_t u, v;
    int64_t w;

    // Stage 2: Relax all edges
    for (uint64_t i = 0; i < vertices - 1; i++)
    {
        for (edge temp_edge : edges)
        {
            u = temp_edge.get_start();
            v = temp_edge.get_end();
            w = temp_edge.get_weight();
            if (dist.at(u) + w < dist.at(v))
            {
                dist.at(v) = dist.at(u) + w;
                pred.at(v) = cast(u);
            }
        }
    }

    // Stage 3 : check for negative cycles; can not find shortest path if there exists a negative cycle
    for (edge temp_edge : edges)
    {
        u = temp_edge.get_start();
        v = temp_edge.get_end();
        w = temp_edge.get_weight();
        if (dist.at(u) + w < dist.at(v))
        {

            // At this point we know that a negative cycle exists;
            // No shortest paths can be found when we have a negative cycle
            throw invalid_argument("Negative Cycle in the Graph!");
        }
        vector<vector<int64_t>> answer;
        answer.push_back(dist);
        answer.push_back(pred);
    }
    vector<vector<int64_t>> ans;
    ans.push_back(dist);
    ans.push_back(pred);
    return ans;
}

/**
 * @brief As Bellman Ford is a single source algorithm, we use this function to apply Bellman Ford to all vertices,
 *        hence having all shortest paths within the whole graph.
 *
 * @param g the graph
 */
void bellman_ford_all_paths(const graph &g, ofstream &out)
{
    vector<edge> edges = g.get_edges();
    uint64_t vertices = g.get_vertices();
    vector<uint64_t> dist(vertices);
    vector<uint64_t> pred(vertices);
    vector<vector<int64_t>> ans;
    bool path_found;
    vector<uint64_t> path;
    out << "\n----------------------------------------------------";
    out << "\nBELLMAN FORD Algorithm has calculated the shortest paths .\n"
        << "The paths and their corresponding cost are as follows : \n";
    for (uint64_t i = 0; i < vertices; i++)
    {
        ans = bellman_ford(g, i);
        vector<int64_t> i_dist = ans.at(0);

        for (uint64_t j = 0; j < vertices; j++)
        {
            if (i != j)
            {

                path = find_path(ans, i, j);

                path_found = found_path(path);
                if (path_found)
                {
                    out << "\nShortest path found from source " << i << " to sink " << j << " with cost " << i_dist.at(j) << " --> "
                        << path << endl;
                }
                else
                {
                    out << "\nNo path from " << i << " to " << j << " has been found!" << endl;
                }
            }
        }
    }
    out << "----------------------------------------------------";
}

/**
 * @brief Dijkstra Algorithm is a single source shortest path algorithm. The initialization is similar to Bellman Ford, but instead of relaxing all edges,
 * it finds the vertex with minimum distance and checks the neighbors of that vertex.
 *
 * @param g he graph
 * @param source source vertex
 * @return vector < vector <uint64_t > >   This vector contains of 2 vectors of type uint64_t. The first vector (index 0) is the vector that has the distances for the source to all other node,
 * and the second vector (index 1) contains the predecessors of the vertices.
 */
vector<vector<int64_t>> dijkstra(const graph &g, const uint64_t &source)
{
    vector<edge> edges = g.get_edges();
    uint64_t vertices = g.get_vertices();
    vector<vector<edge>> adj = g.get_adj();
    vector<int64_t> dist(vertices);
    vector<int64_t> pred(vertices);
    vector<uint64_t> unvisited;

    // Stage 1
    for (uint64_t i = 0; i < vertices; i++)
    {
        dist.at(i) = numeric_limits<int>::max();
        pred.at(i) = -1;
        unvisited.push_back(i);
    }
    dist.at(source) = 0;

    uint64_t min_dist_vertex;
    // uint64_t min_ind;
    uint64_t current_dist;

    while (!unvisited.empty())
    {
        current_dist = numeric_limits<int>::max();

        // find the vertex with minimum distance
        for (uint64_t i = 0; i < unvisited.size(); i++)
        {
            if (dist.at(unvisited.at(i)) < cast(current_dist))
            {
                min_dist_vertex = unvisited.at(i);
                // min_ind = i;
            }
        }

        unvisited.erase(std::remove(unvisited.begin(), unvisited.end(), min_dist_vertex), unvisited.end());

        uint64_t current_end;
        uint64_t temp_dist;

        // check the min_dist_vertex's neighbors
        for (edge e : adj.at(min_dist_vertex))
        {
            current_end = e.get_end();
            if (find(unvisited.begin(), unvisited.end(), current_end) != unvisited.end()) // if it was found
            {
                temp_dist = dist.at(min_dist_vertex) + e.get_weight();
                if (temp_dist < cast(dist.at(current_end)))
                {
                    dist.at(current_end) = cast(temp_dist);
                    pred.at(current_end) = cast(min_dist_vertex);
                }
            }
        }
        uint64_t infinity_count = 0;

        for (uint64_t i = 0; i < unvisited.size(); i++)
        {
            if (dist.at(unvisited.at(i)) == numeric_limits<int>::max())
            {
                infinity_count++;
            }
        }
        // if none of the unvisited vertices are not updated, we should break the loop, as continuing it would not change the results.
        if (infinity_count == unvisited.size())
        {
            break;
        }
    }

    vector<vector<int64_t>> ans;
    ans.push_back(dist);
    ans.push_back(pred);
    return ans;
}

/**
 * @brief As Dijkstra is a single source algorithm, we use this function to apply Bellman Ford to all vertices,
 *        hence having all shortest paths within the whole graph.
 *
 * @param g
 */
void dijkstra_all_paths(const graph &g, ofstream &out)
{
    // In Dijkstra Algorithms weights can not be negative.

    if (!g.non_negative_weights())
    {
        throw invalid_argument("\nGraph used can not have negative weights for Dijkstra Algorithm! Please try again!\n");
    }
    vector<edge> edges = g.get_edges();
    uint64_t vertices = g.get_vertices();
    vector<uint64_t> dist(vertices);
    vector<uint64_t> pred(vertices);
    vector<vector<int64_t>> ans;
    bool path_found;
    vector<uint64_t> path;
    out << "\n----------------------------------------------------\n";
    out << "[DIJKSTRA] has calculated the shortest paths .\n"
        << "The paths and their corresponding cost are as follows : \n";
    for (uint64_t i = 0; i < vertices; i++)
    {
        ans = dijkstra(g, i);
        vector<int64_t> i_dist = ans.at(0);

        for (uint64_t j = 0; j < vertices; j++)
        {
            if (i != j)
            {

                path = find_path(ans, i, j);

                path_found = found_path(path);
                if (path_found)
                {
                    out << "\nShortest path found from source " << i << " to sink " << j << " with cost " << i_dist.at(j) << " --> "
                        << path << endl;
                }
                else
                {
                    out << "\nNo path from " << i << " to " << j << " has been found!" << endl;
                }
            }
        }
    }
    out << "\n----------------------------------------------------";
}

/**
 * @brief Floyd Warshall is an algorithm to find all shortest paths between vertices using dynamic programming.
 *
 * @param g the graph
 * @param out ostream
 */
void floyd_warshall(const graph &g, ofstream &out)
{

    vector<edge> edges = g.get_edges();
    uint64_t vertices = g.get_vertices();
    vector<vector<int64_t>> dist(vertices, vector<int64_t>(vertices, 0));
    vector<vector<int64_t>> next(vertices, vector<int64_t>(vertices, 0));

    // Stage 1
    // Initialize the dist and pred arrays

    for (uint64_t i = 0; i < vertices; i++)
    {
        for (uint64_t j = 0; j < vertices; j++)
        {
            dist[i][j] = numeric_limits<int>::max();
        }
    }

    for (uint64_t i = 0; i < vertices; i++)
    {
        for (uint64_t j = 0; j < vertices; j++)
        {
            next[i][j] = -1;
        }
    }

    for (edge e : edges)
    {
        dist[e.get_start()][e.get_end()] = e.get_weight();
        next[e.get_start()][e.get_end()] = e.get_end();
    }

    for (uint64_t i = 0; i < vertices; i++)
    {
        dist[i][i] = 0;

        next[i][i] = cast(i);
    }

    for (uint64_t k = 0; k < vertices; k++)
    {
        for (uint64_t i = 0; i < vertices; i++)
        {
            for (uint64_t j = 0; j < vertices; j++)
            {
                if (dist[i][j] > dist[i][k] + dist[k][j])
                {
                    // this is the main formula of the algorithm
                    dist[i][j] = dist[i][k] + dist[k][j];
                    next[i][j] = next[i][k];
                }
            }
        }
    }

    out << "\n----------------------------------------------------"
        << "\n[FLOYD WARSHALL]" << endl;
    for (uint64_t i = 0; i < vertices; i++)
    {
        for (uint64_t j = 0; j < vertices; j++)
        {
            if (i != j)
            {
                vector<uint64_t> path;
                if (next[i][j] == -1)
                {
                    out << "\nNO PATH FOUND from " << i << " to " << j << endl;
                    continue;
                    ;
                }
                path.push_back(i);

                uint64_t node = i;
                while (node != j)
                {
                    node = cast(next[node][j]);
                    path.push_back(node);
                }

                out << "\n Shortest path from " << i << " to " << j << " with cost " << dist[i][j] << " -->";
                out << path << endl;
            }
        }
    }
    out << "\n----------------------------------------------------";
}

/**
 * @brief main program
 *
 * @return int
 */
int main(int argc, char *argv[])
{
    cout << "\n\n-------------------------------------------------------------------\n"
         << "** Welcome to calculate the shortest paths in weighted directed graphs!\n"
         << "** This program calculates the shortest paths in weighted directed graphs using 3 main algorithms.\n"
         << "** 1. BELLMAN FORD     2. DIJKSTRA     3. FLOYD WARSHALL\n"
         << "** Please use the following format for command line --> [./executable_file_name][input_filename][output_filename][Algorithm Number (1,2,3)].\n"
         << "** If a number is not specified all three algorithms will be executed.\n\n"
         << "------------------------------------------------------------------------\n\n";

    /** If user does not specify which algorithm to use there should be 3 arguments : exe file, input file name, and output file name.
     * If user does specify the algorithm, there should be 4 arduments: exe file, input file name, output file name, algorithm number.
     */
    if (!(argc == 3 || argc == 4)) // Wrong number of arguments
    {
        cout << "argc is " << argc << endl;
        cout << "ERROR! You entered the wrong number of inputs in command line.\n"
             << "Please use the following format for command line --> [./executable_file_name][input_filename][output_filename]\n\n";
        return -1;
    }

    string input_filename = argv[1];
    string output_filename = argv[2];
    int64_t algo_num = -1;
    if (argc == 4)
    {
        algo_num = atoi(argv[3]);
        if (!(algo_num == 1 || algo_num == 2 || algo_num == 3))
        {
            cout << "ERROR! The number specified for the algorithm is not valid! Please try again using only 1, 2, or 3. \n";
            return -1;
        }
    }

    /*< Checking if input and output files are opened correctly>*/
    ifstream in(input_filename);

    if (!in.is_open())
    {
        cout << "Error opening input file!";
        return -1;
    }

    ofstream out(output_filename);

    if (!out.is_open())
    {
        cout << "Error opening output file!";
        return -1;
    }

    try
    {

        graph g = graph(input_filename);
        out << "\n** OUTPUT FOR THE GIVEN GRAPH IS AS FOLLOWS: ** \n";
        if (algo_num == 1)
        {
            timer t_bellman_ford;
            bellman_ford_all_paths(g, out);
            t_bellman_ford.end();
            cout << "\nBELLMAN FORD has successfully executed."
                 << "\n[BELLMAN FORD] Elapsed time: " << t_bellman_ford.seconds() << " milli seconds.\n\n";
            out << "\n[BELLMAN FORD] Elapsed time: " << t_bellman_ford.seconds() << " milli seconds.\n\n";
        }

        else if (algo_num == 2)
        {
            timer t_dijkstra;
            dijkstra_all_paths(g, out);
            t_dijkstra.end();
            cout << "\nDIJKSTRA has successfully executed."
                 << "\n[DIJKSTRA] Elapsed time: " << t_dijkstra.seconds() << " milli seconds.\n\n";
            out << "\n[DIJKSTRA] Elapsed time: " << t_dijkstra.seconds() << " milli seconds.\n\n";
        }

        else if (algo_num == 3)
        {
            timer t_floyd_warshall;
            floyd_warshall(g, out);
            t_floyd_warshall.end();
            cout << "\nFLOYD WARSHALL has successfully executed."
                 << "\n[FLOYD WARSHALL] --> Elapsed time: " << t_floyd_warshall.seconds() << " milli seconds.\n\n";
            out << "\n[FLOYD WARSHALL] --> Elapsed time: " << t_floyd_warshall.seconds() << " milli seconds.\n\n";
        }

        else if (algo_num == -1)
        { // No algorithm number wa specified by user; All algorithms will be generated.
            timer t_bellman_ford;
            bellman_ford_all_paths(g, out);
            t_bellman_ford.end();
            cout << "\nBELLMAN FORD has successfully executed."
                 << "\nOutput is written to file."
                 << "\n[BELLMAN FORD] Elapsed time: " << t_bellman_ford.seconds() << " milli seconds.\n\n";
            out << "\n[BELLMAN FORD] Elapsed time: " << t_bellman_ford.seconds() << " milli seconds.\n\n";

            timer t_dijkstra;
            dijkstra_all_paths(g, out);
            t_dijkstra.end();
            cout << "\nDIJKSTRA has successfully executed."
                 << "\nOutput is written to file."
                 << "\n[DIJKSTRA] Elapsed time: " << t_dijkstra.seconds() << " milli seconds.\n\n";
            out << "\n[DIJKSTRA] Elapsed time: " << t_dijkstra.seconds() << " milli seconds.\n\n";

            timer t_floyd_warshall;
            floyd_warshall(g, out);
            t_floyd_warshall.end();
            cout << "\nFLOYD WARSHALL has successfully executed."
                 << "\nOutput is written to file."
                 << "\n[FLOYD WARSHALL] --> Elapsed time: " << t_floyd_warshall.seconds() << " milli seconds.\n\n";
            out << "\n[FLOYD WARSHALL] --> Elapsed time: " << t_floyd_warshall.seconds() << " milli seconds.\n\n";
        }
    }

    catch (exception &e)
    {
        cout << "ERROR! " << e.what() << "\n";
    }
}