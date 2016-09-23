#ifndef COGRAPH_EDIT
#define COGRAPH_EDIT

#include <iostream>
#include <string>
#include <set>
#include <utility>                   // for std::pair
#include <algorithm>                 // for std::for_each
#include <limits> // for max_int
#include <time.h> 
#include <math.h>   
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/subgraph.hpp>
#include <boost/config.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/copy.hpp>
#include "ModDecomp.h"

class Cograph_Editing {
public:
  Graph cograph_editing (Graph& g);
private:
  Tree editTree(const Tree& t, const TVertex& v, Graph& g);
  Tree editTreeRefined(const Tree& t, const TVertex& v, Graph& g);
};

#endif