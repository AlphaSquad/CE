#ifndef BUILDGRAPHS_H
#define BUILDGRAPHS_H

#include <iostream>
#include <string>
#include <set>
#include <random> 		// for levels of graph
#include <utility>                   // for std::pair
#include <algorithm>                 // for std::for_each
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/subgraph.hpp>
#include <boost/config.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/copy.hpp>
#include <boost/graph/connected_components.hpp>
#include "Commons.h" 

struct root_t {
  typedef boost::vertex_property_tag kind;
};

typedef boost::adjacency_list<boost::listS, boost::listS, boost::bidirectionalS, boost::property<boost::vertex_name_t, std::string, boost::property<root_t ,bool> >, boost::vecS> Cotree;
typedef typename boost::graph_traits<Cotree>::edge_descriptor CEdge;
typedef typename boost::graph_traits<Cotree>::vertex_descriptor CVertex;
typedef boost::graph_traits<Cotree>::vertex_iterator cvertex_iter;
typedef boost::graph_traits<Cotree>::edge_iterator cedge_iter;
typedef boost::graph_traits<Cotree>::out_edge_iterator c_out_edge_iter; // OUT! JEEZ!
typedef boost::graph_traits<Cotree>::adjacency_iterator caiter;


class BuildGraphs {
public:
  BuildGraphs(unsigned int seed = 0, bool print = true);
  Graph buildCograph (unsigned int vertices);
  void graphProperties(std::vector<float>& prop, const Tree& t, const Graph& g);
  Graph addEdges(Graph& g, float pc);
  Graph removeEdges(Graph& graph, float pc);
  Graph editEdges(Graph& graph, float pc);
  void writeGraphs (const std::vector<Graph>& graphs, std::string filename, std::vector<std::pair<int,int>> props = {});
  void writeGraph (const Graph& g, std::string filename);
  void ReadEdgeList(std::string sFileName, Graph& g, bool& FileExists);
  std::string getNewick(const Tree& t);
  std::vector<Graph> readGraphs(std::string filename);
  int calculateEditDistance(const Graph& g1, const Graph& g2);
private:
  Cotree buildCotree (unsigned int vertices);
  std::string getSubtree(const Tree& t, const TVertex& v);
  Graph graphFromTree(const Cotree& cotree);
  int rootedTree(int toAdd, Cotree& t, CVertex& v, bool series, int num);
  
  std::mt19937 generator;
};

#endif
