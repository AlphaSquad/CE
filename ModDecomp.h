#ifndef MODDECOMP_H
#define MODDECOMP_H

#include <iostream>
#include <string>
#include <set>
#include <utility>                   // for std::pair
#include <algorithm>                 // for std::for_each
#include <limits> // for max_int
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/subgraph.hpp>
#include <boost/config.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/copy.hpp>
#include <boost/graph/connected_components.hpp>
#include "DyckWord.h"

class ModDecomp {
public:
  Tree decompose(Graph& g);
  Tree decompose_components(Graph& g);
  Tree calcModDecomp(const std::vector<Vertex>& factorization, Graph& graph);
  std::vector<Vertex> calcFacPerm(const Graph& graph);
  //stepwise? 
private:
  DyckWord parenthesizing(const std::vector<Vertex>& fac, const Graph& graph, std::map<Vertex, int>& lcutters, std::map<Vertex, int>& rcutters);
  Tree buildFractureTree(const DyckWord& dyck, Graph& g, const std::map<Vertex, int>& lcutters, const std::map<Vertex, int>& rcutters); 
  void moduleDetDelMerge(); //Detect modules, delete dummies, merge modules [one step?]  
  void buildModDecomp(const Graph& graph, const std::vector<Vertex>& fac, const std::map<Vertex, int>& lcutters, const std::map<Vertex, int>& rcutters, Tree& fractureTree); // merging moved to here  
};

#endif //MODDECOMP_H