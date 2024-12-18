#include "ModDecomp.h"

// These functions provide begin and end implementations for the BGL iterator
// pairs. Thus we can simply use for(auto e : edges(g)) {} instead of messy
// iterator syntax.
namespace std
{
	template <typename A> A begin(const pair<A, A>& s)
	{
		return s.first;
	}

	template <typename A> A end(const pair<A, A>& s)
	{
		return s.second;
	}
}


//DEBUG
void printTree(const Tree& fractureTree, bool toDelete = false, bool name = false, bool corrV = false)
{
  typedef std::map<TVertex, int> IndexMap;
  IndexMap mapIndex;
  boost::associative_property_map<IndexMap> propmapIndex(mapIndex);
  tvertex_iter vi, vi_end;
  int i = 0; // index
  for (boost::tie(vi,vi_end) = boost::vertices(fractureTree); vi != vi_end; ++vi)
    boost::put(propmapIndex,*vi,i++);
  
  std::cout << boost::num_vertices(fractureTree) << std::endl;
  std::cout << boost::num_edges(fractureTree) << std::endl;
  
  if (toDelete)
    boost::write_graphviz(std::cout, fractureTree, boost::make_label_writer(boost::get(&vertex_prop::toDelete,fractureTree)), boost::default_writer(), boost::default_writer(), propmapIndex);
  else if (name)
    boost::write_graphviz(std::cout, fractureTree, boost::make_label_writer(boost::get(&vertex_prop::name,fractureTree)), boost::default_writer(), boost::default_writer(), propmapIndex);
  else if (corrV)
    boost::write_graphviz(std::cout, fractureTree, boost::make_label_writer(boost::get(&vertex_prop::corrV,fractureTree)), boost::default_writer(), boost::default_writer(), propmapIndex);
  else
    boost::write_graphviz(std::cout, fractureTree, boost::make_label_writer(boost::get(&vertex_prop::type,fractureTree)), boost::default_writer(), boost::default_writer(), propmapIndex);
}

//DEBUG
void printGraph(const Graph& g)
{
  boost::write_graphviz(std::cout, g);
}

template<typename T>
int findFractures(T iter, T siter, const std::vector<std::string>& cutset, const Graph& graph, std::map<Vertex, int>& lbrackets, std::map<Vertex, int>& rbrackets) {
    
    const auto& name = boost::get(boost::vertex_name, graph);
    const auto& facPos = boost::get(facPos_t(), graph);
    
    for (; iter != siter; ++iter) // iter != siter: we cannot be a left fracture if right of pair
    {
      if (std::binary_search(cutset.begin(), cutset.end(), boost::get(name,*iter)))
      {
	//left fracture found
	//insert closing bracket right of siter
	++rbrackets[*siter];
	++lbrackets[*iter];

	return boost::get(facPos,*iter);
      }
    }
    
    return -1;
}

// Given a factorizing permutation constructs the dyck word (parenthesized factorization) depending on the graph
DyckWord ModDecomp::parenthesizing(const std::vector<Vertex>& fac, const Graph& graph, std::map<Vertex, int>& lcutters, std::map<Vertex, int>& rcutters) // make sure the graph names are the names appearing in the string
{
  const auto& name = boost::get(boost::vertex_name, graph);

  // we count for every string opening and closing brackets before/after vertices to build the dyck word from it later on
  std::map<Vertex, int> lbrackets;
  std::map<Vertex, int> rbrackets;
  
  for(size_t i = 0; i < fac.size() - 1; ++i) {
    rcutters[fac[i]] = -1;
    lcutters[fac[i]] = fac.size() + 1;
    
    std::vector<std::string> nh1;
    for(const auto& v : boost::adjacent_vertices(fac[i], graph)) {
      nh1.push_back(boost::get(name, v));
    }
    
    std::vector<std::string> nh2;
    for(const auto& v : boost::adjacent_vertices(fac[i + 1], graph)) {
      nh2.push_back(boost::get(name, v));
    }
    
    std::sort(nh1.begin(), nh1.end());
    std::sort(nh2.begin(), nh2.end());

    std::vector<std::string> cutset;
    
    std::set_symmetric_difference(
      nh1.begin(), nh1.end(),
      nh2.begin(), nh2.end(),
      std::back_inserter(cutset)
    );

    int r = findFractures(fac.begin(), fac.begin() + i, cutset, graph, lbrackets, rbrackets);

    // small hack so we can use the same function for l and r brackets
    if(r != -1) {
      lcutters[fac[i]] = r;
    }

    rcutters[fac[i]] = findFractures(fac.rbegin(), fac.rbegin() + fac.size() - i - 2, cutset, graph, rbrackets, lbrackets);
  }
  
  rcutters[fac[fac.size() - 1]] = -1;
  lcutters[fac[fac.size() - 1]] = fac.size() + 1;
  // The last vertex does not have a pair and thus does not have any cutters whatsoever

   //Now build the dyck word
   return DyckWord(graph, lbrackets, rbrackets, fac);
}

// add child to node curr with lfrac size
TVertex addChild(const TVertex& curr, int size, Tree& fractureTree)
{
  TVertex v = boost::add_vertex(fractureTree);
  fractureTree[v].name = "Internal"; // These will later on become the module nodes [series/parallel/prime]
  fractureTree[v].lfrac = size;
  fractureTree[v].rfrac = -1;
  fractureTree[v].toDelete = false; // we need to initalize this, might be set to "true" in "leaveChild"
  fractureTree[v].root = false;
  boost::add_edge(curr,v,fractureTree);
  return v;
}

// does internal vertex curr has to be deleted
bool deleteable(TVertex curr, Vertex lastChild, const Graph& g, const Tree& fractureTree)
{
  const auto& facPos = boost::get(facPos_t(), g);
  return !fractureTree[curr].toDelete and (fractureTree[curr].rfrac > boost::get(facPos,lastChild) or fractureTree[curr].lfrac < boost::get(facPos,fractureTree[curr].containedNodes[0]));
}

// recurr back to father node of curr as curr has been finished
TVertex leaveChild(const TVertex& curr, const DyckWord& dyck, int i, const Graph& g, const std::map<Vertex, int>& lcutters, const std::map<Vertex, int>& rcutters, Tree& fractureTree)
{
  const auto& facPos = boost::get(facPos_t(), g);
  
  int pos = fractureTree[curr].containedNodes.size();
  Vertex lastChild = fractureTree[curr].containedNodes[pos - 1];
  if (boost::out_degree(curr,fractureTree) == 1) // we only have one child: dummy node 
    fractureTree[curr].toDelete = true;
  
  const auto& name = boost::get(boost::vertex_name, g);
  fractureTree[curr].corrV = boost::get(name, fractureTree[curr].containedNodes[0]); // initialise
  fractureTree[curr].rfrac = std::max(boost::get(facPos,lastChild), fractureTree[curr].rfrac); // we can only set this now

  if (deleteable(curr, lastChild, g, fractureTree))
    fractureTree[curr].toDelete = true; // this node is a dummy node
    
  int ltemp = fractureTree[curr].lfrac;
  int rtemp = fractureTree[curr].rfrac;
  std::vector<Vertex> nodes = fractureTree[curr].containedNodes;
  
  TVertex newCurr = boost::source(*(boost::in_edges(curr,fractureTree).first),fractureTree); // Assume dyck word is well-formed

  assert(curr != newCurr);

  // we leave this node for good, give its cutters to parent node
  fractureTree[newCurr].containedNodes.insert(fractureTree[newCurr].containedNodes.end(), nodes.begin(), nodes.end()); // save all the contained Nodes, even if they're not direct children
  if (i + 1 < dyck.size() and dyck.getToken(i + 1).getType() != Type::RBracket) // If we do not leave the node in the next step  we need to consider the cutters of the last child
  {
    fractureTree[newCurr].lfrac = std::min(std::min(lcutters.at(lastChild), fractureTree[newCurr].lfrac), ltemp); // leftmost cutter
    fractureTree[newCurr].rfrac = std::max(std::max(rcutters.at(lastChild), fractureTree[newCurr].rfrac), rtemp); // rightmost cutter
  }
  else
  { // this was really the last node, we must not consider the cutters of the last child
    fractureTree[newCurr].lfrac = std::min(fractureTree[newCurr].lfrac, ltemp); // leftmost cutter
    fractureTree[newCurr].rfrac = std::max(fractureTree[newCurr].rfrac, rtemp); // rightmost cutter
  }
  
  return newCurr;
}

// a vertex corresponding to the original graph has been found, add as leaf
void addLeaf(const TVertex& curr, const Token& token, const DyckWord& dyck, int i, Graph& g, const std::map<Vertex, int>& lcutters, const std::map<Vertex, int>& rcutters, Tree& fractureTree)
{ 
//   assert(token.getType() == Token::Node);

  const auto& cv = boost::get(corrV_t(), g); // this is not really the color but it does work [stores the corresponding vertex in the graph]
  const auto& facPos = boost::get(facPos_t(), g);
  const auto& name = boost::get(boost::vertex_name, g);
  
  TVertex leaf = boost::add_vertex(fractureTree);
  
  Vertex v = token.getVertex();
  
  fractureTree[leaf].name = token.toString(g);
  fractureTree[leaf].type = token.toString(g); // type in Modular Decomposition = name in original graph
  fractureTree[leaf].lfrac = dyck.size() + 1;
  fractureTree[leaf].rfrac = -1;
  fractureTree[leaf].toDelete = false; // if this is uninitialized, strange things may happen
  fractureTree[leaf].root = false;
  fractureTree[leaf].corrV = boost::get(name,v); // go over name
  fractureTree[leaf].containedNodes = std::vector<Vertex>{v};
  boost::put(cv, v, leaf);
   
  fractureTree[curr].children.push_back(v); // direct children
  fractureTree[curr].containedNodes.push_back(v); //we need this for module detection
  
  boost::add_edge(curr,leaf,fractureTree);
  //
  if (boost::out_degree(curr,fractureTree) == 1) // we're the first node added
    fractureTree[curr].lfrac = std::min(boost::get(facPos,v),std::min(fractureTree[curr].lfrac, lcutters.at(v))); //is firstChild the minimum?
  
  if (i + 1 < dyck.size() && dyck.getToken(i + 1).getType() != Type::RBracket) //iter + 1 always exists as the last node is always ")". This means, we are not the last vertex added
  {
    fractureTree[curr].lfrac = std::min(fractureTree[curr].lfrac, lcutters.at(v)); // leftmost cutter, we now can leave firstChild out
    fractureTree[curr].rfrac = std::max(fractureTree[curr].rfrac, rcutters.at(v)); // rightmost cutter
  }
  
  // If we are the last vertex added: we must not consider fractures going "out" of the module
}

// Given the Dyck word builds the fracture tree
Tree ModDecomp::buildFractureTree(const DyckWord& dyck, Graph& g, const std::map<Vertex, int>& lcutters, const std::map<Vertex, int>& rcutters)
{
  // create root
  Tree fractureTree;
  TVertex root = boost::add_vertex(fractureTree);
  fractureTree[root].name = "Internal";
  fractureTree[root].root = true;
  fractureTree[root].toDelete = false;
  TVertex curr = root;
  for (int i = 0; i < dyck.size(); i++){
    Token t = dyck.getToken(i);
    switch (t.getType())
    {
      case Type::LBracket : curr = addChild(curr, dyck.size() + 1, fractureTree); break;
      case Type::RBracket : curr = leaveChild(curr, dyck, i, g, lcutters, rcutters, fractureTree); break;
      case Type::Node : addLeaf(curr, t, dyck, i, g, lcutters, rcutters, fractureTree); break;
    }
  }
  if (boost::out_degree(root,fractureTree) == 1)
    fractureTree[root].toDelete = true; 
  return fractureTree;
}

// delete the marked vertex v
void deleteVertex(TVertex& v, Tree& fractureTree)
{
  auto parent = *boost::inv_adjacent_vertices(v, fractureTree).first;

  std::vector<TEdge> toDelete;
  for (auto ei : boost::out_edges(v, fractureTree)) // remove all edges from node
  {
    boost::add_edge(parent, boost::target(ei,fractureTree), fractureTree); // link children to parent of vertex to be removed
    toDelete.push_back(ei);
  }
  
  for (const auto& e : toDelete)
    boost::remove_edge(e, fractureTree);
  
  fractureTree[parent].children.insert(
    fractureTree[parent].children.end(),
    fractureTree[v].children.begin(),
    fractureTree[v].children.end()); // add direct children to father

  // The children should be in "contained Nodes" anyways
  boost::remove_edge(parent,v,fractureTree); // remove edge to node  
  boost::remove_vertex(v,fractureTree); // all edges are already gone
}

// Detects all modules and deletes dummy nodes which do not represent a module 
void moduleDetDel(Tree& fractureTree) //Detect modules, delete dummies 
{ // This craps on my name mapping
  // Delete all nodes which flags have been set in the previous step
  tvertex_iter vi, vi_end, next;
  boost::tie(vi,vi_end) = boost::vertices(fractureTree);
  for (next = vi; vi != vi_end; vi = next)
  {
    ++next;
    if (fractureTree[*vi].toDelete && !fractureTree[*vi].root) // we cant do this if we are root (link to parents)
    {    
      deleteVertex(*vi, fractureTree);
    }
    else if (fractureTree[*vi].toDelete && fractureTree[*vi].root) // this is only the case iff root has exactly one child. Otherwise the root is always a module: The one containing all nodes
    {
      auto child = *boost::adjacent_vertices(*vi, fractureTree).first;
      
      fractureTree[child].root = true; // we can do this
      boost::remove_edge(*(boost::out_edges(*vi,fractureTree).first),fractureTree); // remove _the_ out edge (can only be one)
      // we dont link child to parent because root does not have parent
      boost::remove_vertex(*vi,fractureTree); // there are no in-edges so we can delete the old root now
    }
  }
}

// not really LCA, get direct child of v which is an ancestor of curr
TVertex getLCA(const TVertex& curr, const TVertex& v, const Tree& fractureTree)
{
  TVertex tmp = curr;
  while (*boost::inv_adjacent_vertices(tmp,fractureTree).first != v)
    tmp = *boost::inv_adjacent_vertices(tmp,fractureTree).first;
  return tmp;
}

// merging DOES occur, find the modules to merge
bool getOriginalModule(const SubGraph& rootG, const std::vector<Vertex>& fac, const TVertex& v, const Vertex& v1, const Vertex& v2, TVertex& tA1, TVertex& tA2, const std::map<Vertex, int>& lcutters, const std::map<Vertex, int>& rcutters, Tree& fractureTree)
{
  const auto& cv = boost::get(corrV_t(), rootG);
  const auto& facPos = boost::get(facPos_t(), rootG);
  
  // real children of *ti?
  bool child1 = std::find(fractureTree[v].children.begin(), fractureTree[v].children.end(),v1) != fractureTree[v].children.end();
  bool child2 = std::find(fractureTree[v].children.begin(), fractureTree[v].children.end(),v2) != fractureTree[v].children.end();
  
  /* This should work as well
   * 
  bool child1 = std::find(fractureTree[v].children.begin(), fractureTree[v].children.end(),v1) != fractureTree[v].children.end();
  bool child2 = std::find(fractureTree[v].children.begin(), fractureTree[v].children.end(),v2) != fractureTree[v].children.end();
  
  TVertex curr1 = child1 ? boost::get(cv,v1) : getLCA(curr1,v,fractureTree);
  TVertex curr2 = child2 ? boost::get(cv,v2) : getLCA(curr2,v,fractureTree);
  
  Vertex lastVertex = *(fractureTree[curr1].containedNodes.rbegin());

  int posF = child1 ? lcutters.at(v1) : lcutters.at(lastVertex);
  int posL = child1 ? rcutters.at(v1) : rcutters.at(lastVertex);  
  
  bool rep1 = posF >= (child1 ? boost::get(facPos,v2) : boost::get(facPos,*(fractureTree[curr1].containedNodes.begin())) );
  bool rep2 = posL <= (child2 ? boost::get(facPos,v2) : boost::get(facPos,*(fractureTree[curr2].containedNodes.rbegin())) );
  
  bool rep = rep1 and rep2;
   * 
   */
  TVertex curr1 = boost::get(cv,v1);
  TVertex curr2 = boost::get(cv,v2);
  bool rep = false;
  int posF;
  int posL;  
  
  // which of the two vertices is a "real" child of the module = leaf as direct child
  if (child1 and !child2)
  {
    curr2 = getLCA(curr2,v,fractureTree);
    posF = lcutters.at(v1); 
    posL = rcutters.at(v1);// only v1 is a real child should not have any cutters with v2
    rep = (posF >= boost::get(facPos,v2) and posL <= boost::get(facPos,*(fractureTree[curr2].containedNodes.rbegin())));
  }
  else if (child2 and !child1)
  {
    curr1 = getLCA(curr1,v,fractureTree);
    Vertex lastVertex = *(fractureTree[curr1].containedNodes.rbegin());
    posF = lcutters.at(lastVertex); 
    posL = rcutters.at(lastVertex);
    rep = (posF >= boost::get(facPos,*(fractureTree[curr1].containedNodes.begin())) and posL <= boost::get(facPos,v2));
  }
  else if (!child1 and !child2)
  {
    curr1 = getLCA(curr1,v,fractureTree);
    curr2 = getLCA(curr2,v,fractureTree);
    Vertex lastVertex = *(fractureTree[curr1].containedNodes.rbegin()); 
    posF = lcutters.at(lastVertex); 
    posL = rcutters.at(lastVertex);   
    rep = (posF >= boost::get(facPos,*(fractureTree[curr1].containedNodes.begin())) and posL <= boost::get(facPos,*(fractureTree[curr2].containedNodes.rbegin())));
  }
  else
  {
    posF = lcutters.at(v1);
    posL = rcutters.at(v1); // both are real children, left and right cutter between v1 and v2 must not exist
    rep = false; // just to be sure that posF = fac.size() +1 and posL = -1 has to be true
  }
  
  bool fc = (posF == fac.size() + 1);
  bool lc = (posL == -1);
  
  // set vertices which have to be added/merged
  if (curr1 != curr2 and ((fc and lc) or rep))
  {
    tA1 = curr1;
    tA2 = curr2;
    return true;
  }
  
  return false;
}

// we have detected a case where merging might occur/module type has to be assigned
void detectMerging(const SubGraph& rootG, const Graph& graph, const TVertex& v, const std::vector<Vertex>& fac, int minPos, int maxPos, const std::map<Vertex, int>& lcutters, const std::map<Vertex, int>& rcutters, Tree& fractureTree)
{
  // temporaries
  const auto& cv = boost::get(corrV_t(), rootG);
  bool added = false;
  TVertex merged;
  std::vector<TVertex> newModules;
  std::vector<std::pair<TVertex,TVertex> > toDelete;
  std::vector<std::pair<TVertex,TVertex> > toAdd;
  
  
  for (std::vector<Vertex>::const_iterator it = fac.begin() + minPos; it != fac.begin() + maxPos; ++it) // my beloved paris [sic!]
  {
    bool foundNew = false;
    
    TVertex tA1;
    TVertex tA2;
    
    foundNew = getOriginalModule(rootG, fac, v, *it, *(it+1), tA1, tA2, lcutters, rcutters, fractureTree); // true if merging occurs

    if (foundNew)
    {
      if (!added)
      {
	merged = boost::add_vertex(fractureTree);
	fractureTree[merged].name = "Internal"; // we do not necessarily need to initalize the other values
	fractureTree[merged].root = false;
	if (boost::edge(*it,*(it + 1),graph).second)
	  fractureTree[merged].type = "Series";
	else
	  fractureTree[merged].type = "Parallel";
	//Add contained nodes
	if (fractureTree[tA1].name == "Internal")
	  fractureTree[merged].containedNodes.insert(fractureTree[merged].containedNodes.end(), fractureTree[tA1].containedNodes.begin(), fractureTree[tA1].containedNodes.end());
	else
	{
	  fractureTree[merged].children.insert(fractureTree[merged].children.begin(), fractureTree[tA1].children.begin(), fractureTree[tA1].children.end());
	  fractureTree[merged].containedNodes.insert(fractureTree[merged].containedNodes.begin(), fractureTree[tA1].containedNodes.begin(), fractureTree[tA1].containedNodes.end());
	}
	toAdd.push_back(std::make_pair(merged,tA1));
	toDelete.push_back(std::make_pair(v,tA1));
	newModules.push_back(merged);
	added = true;
      }
      //boost::add_edge(merged,t2,fractureTree);
      if (fractureTree[tA2].name == "Internal")
	fractureTree[merged].containedNodes.insert(fractureTree[merged].containedNodes.end(), fractureTree[tA2].containedNodes.begin(), fractureTree[tA2].containedNodes.end());
      else
      {
	fractureTree[merged].children.insert(fractureTree[merged].children.begin(), fractureTree[tA2].children.begin(), fractureTree[tA2].children.end());
	fractureTree[merged].containedNodes.insert(fractureTree[merged].containedNodes.begin(), fractureTree[tA2].containedNodes.begin(), fractureTree[tA2].containedNodes.end());
      }
      toAdd.push_back(std::make_pair(merged,tA2));
      toDelete.push_back(std::make_pair(v,tA2));
      foundNew = false;
    }
    else if (getLCA(boost::get(cv,*it),v, fractureTree) != getLCA(boost::get(cv,*(it+1)),v,fractureTree))
      added = false;
  }
  // add merged module
  for (std::vector<std::pair<TVertex, TVertex> >::iterator it = toAdd.begin(); it != toAdd.end(); ++it)  
    boost::add_edge(it->first, it->second, fractureTree);
  
  for (std::vector<std::pair<TVertex, TVertex> >::iterator it = toDelete.begin(); it != toDelete.end(); ++it)
    boost::remove_edge(it->first,it->second,fractureTree); 

  for (const auto& vert : newModules)
    boost::add_edge(v, vert, fractureTree);
}

// Checks whether a degree sequence corresponds to a thin spider
bool checkThinSpider(unsigned int s, const std::vector<unsigned int>& degrees)
{
  if (degrees[0] != 0)
    return false;
  for (unsigned int k = 2; k < s; k++)
    if (degrees[k] != 0)
      return false;
  for (unsigned int k = degrees.size() - s + 1; k < degrees.size(); k++)
    if (degrees[k] != 0)
      return false;
  return true;
}

bool checkThickSpider(unsigned int s, const std::vector<unsigned int>& degrees)
{
  if (degrees[degrees.size() - 1] != 0)
    return false;
  for (unsigned int k = 0; k < s - 1; k++)
    if (degrees[k] != 0)
      return false; 
  for (unsigned int k = s; k < degrees.size() - s + 1; k++)
    if (degrees[k] != 0)
      return false;
  return true;
}

// does the calculation for module assigning
void calculateModuleType(const Graph& graph, SubGraph& tmp, const TVertex& v, const std::vector<Vertex>& fac, int minPos, int maxPos, const std::map<Vertex, int>& lcutters, const std::map<Vertex, int>& rcutters, Tree& fractureTree)
{
  Tree repGraph; // The representative graph, relevant for prime nodes
  int out_d = -1;
  bool foundChild = false;
  unsigned int max_degree = fractureTree[v].containedNodes.size(); // in a spider the head and legs have to be real children
  std::vector<unsigned int> degrees(max_degree,0); // store degrees for spider recognition (initialise all degrees with 0)
  svertex_iter si, si_end;
  for (boost::tie(si, si_end) = boost::vertices(tmp); si != si_end; ++si)
  {
    const auto& name = boost::get(boost::vertex_name, tmp);
    
    if((size_t)boost::degree(*si,tmp) >= degrees.size()) {
       std::cout << boost::get(name,*si) << " " << max_degree << " " << boost::degree(*si,tmp) << std::endl;
       for (const auto& v : boost::adjacent_vertices(*si,tmp))
       {
	 std::cout << boost::get(name, v) << " ";
       }
       std::cout << std::endl;
       for (const auto& v : boost::vertices(tmp))
       {
	 std::cout << boost::get(name,v) << " ";
       }
       std::cout << std::endl;
    }
    
    assert(boost::degree(*si,tmp) < degrees.size());
    
    degrees[boost::degree(*si,tmp)]++; // degree sequence
    
    Vertex t = tmp.local_to_global(*si); // we need the global name for searching in the fractureTree
    std::vector<Vertex>::iterator child = std::find(fractureTree[v].children.begin(), fractureTree[v].children.end(),t);
    if (child != fractureTree[v].children.end()) // This means, we are a "real" child of the module
    {     
      out_d = boost::degree(*si,tmp);
      foundChild = true;
      if (out_d != 0 and out_d != fractureTree[v].containedNodes.size() - 1) // out degree is neither 0 nor k-1 => we are prime
      {
	fractureTree[v].type = "Prime"; // As the node is prime we need to detect merging
	//break; //we cannot break because of spider recognition
      }
    }
    // If we are no real child we will be considered later on (or already have been)
  }
  if (!foundChild) // We moved out_d outwards. If there is no "real" child it should be -1
    out_d = -1;
  
  bool spider = false;
  // Spider detection. First if case: "Direct" spider: body and legs are direct children, head is arbitrary 
  if (max_degree >= 3 and fractureTree[v].children.size() >= boost::out_degree(v,fractureTree) - 1) // we can only be a spider if we have >= 4 nodes => max_degree is >= 3 and all but one child has to be leaves 
  {
    unsigned int s = degrees[1]; 
    // check for thin spider
    if (s < degrees.size() and s >= 2)
    {
      unsigned int v_s = degrees[degrees.size() - s];
      if (s == v_s) // we possibly are a spider, continue checking
	spider = checkThinSpider(s, degrees);
    }
    if (!spider)
    {
      s = degrees[degrees.size() - 2]; // |V| - 2, to all vertices but 1 leg
      if (s > 0)
      {
	unsigned int v_s = degrees[s - 1]; // s - 1 cannot be larger than degrees.size() - 1
	if (s == v_s) // possibly a thick spider
	  spider = checkThickSpider(s, degrees);
      }
    }
  }
  if (spider)
    fractureTree[v].type = "Spider";
  else if (out_d == 0)
    fractureTree[v].type = "Parallel";
  else if (out_d == fractureTree[v].containedNodes.size() - 1)
    fractureTree[v].type = "Series";
  else if (out_d == -1) // This means we dont have any real children! Assign module type depending on the edges between modules
  {
    TVertex child = *boost::adjacent_vertices(v, fractureTree).first;
    int deg_global = 0;
    for (const auto& w : boost::adjacent_vertices(v,fractureTree))
      if (w != child and boost::edge(fractureTree[child].containedNodes[0],fractureTree[w].containedNodes[0], graph).second) // If they are in a module all or no edges are existent
	deg_global += fractureTree[w].containedNodes.size();
      
    if (deg_global == fractureTree[v].containedNodes.size() - fractureTree[child].containedNodes.size()) // neighboured to all other modules
      fractureTree[v].type = "Series";
    else if (deg_global == 0) // neighboured to no other module
      fractureTree[v].type = "Parallel";
    else // well, something else => can only by prime. neat.
    {
      fractureTree[v].type = "Prime";
      detectMerging(tmp.root(), graph, v, fac, minPos, maxPos, lcutters, rcutters, fractureTree);
    }
  }
  else // detect merging
  {
    fractureTree[v].type = "Prime";
    detectMerging(tmp.root(), graph, v, fac, minPos, maxPos, lcutters, rcutters, fractureTree); // this also assigns the module type if still missing (this is a lie)
  }
  // weak modules cannot exist in undirected graphs - we dont need to detect them!
}

// assigns the type to an internal module
void assignModuleType(SubGraph& rootG, const Graph& graph, TVertex v, const std::vector<Vertex>& fac, const std::map<Vertex, int>& lcutters, const std::map<Vertex, int>& rcutters, Tree& fractureTree)
{
  const auto& facPos = boost::get(facPos_t(), graph);
  SubGraph tmp = rootG.create_subgraph();
  int minPos = boost::num_vertices(rootG); //for module merging if prime
  int maxPos = 0;
  for (const auto& vert : fractureTree[v].containedNodes){
    boost::add_vertex(vert,tmp);  // They are already part!
    minPos = std::min(minPos, boost::get(facPos,vert));
    maxPos = std::max(maxPos, boost::get(facPos,vert));
  }

  assert(boost::num_vertices(tmp) == fractureTree[v].containedNodes.size());
  
  if (boost::num_edges(tmp) == 0) // No edges between children nodes => parallel
    fractureTree[v].type = "Parallel";      
  else if (boost::num_edges(tmp) == boost::num_vertices(tmp) * (boost::num_vertices(tmp) - 1)/2) // all edges between children => series
    fractureTree[v].type = "Series";
  else // Then its more complex
    calculateModuleType(graph, tmp, v, fac, minPos, maxPos, lcutters, rcutters, fractureTree);  
}

// Finally builds the modular decomposition with the reduced fracture tree and the graph structure 
void ModDecomp::buildModDecomp(const Graph& graph, const std::vector<Vertex>& fac, const std::map<Vertex, int>& lcutters, const std::map<Vertex, int>& rcutters, Tree& fractureTree) //we need representative graphs for module merging, factorization for module merging
{ 
  SubGraph rootG;
  boost::copy_graph(graph, rootG);
  
  tvertex_iter ti,ti_end;
  for (boost::tie(ti,ti_end) = boost::vertices(fractureTree); ti != ti_end; ++ti)
  {
    if (fractureTree[*ti].name == "Internal" and fractureTree[*ti].type != "Series" and fractureTree[*ti].type != "Parallel" and fractureTree[*ti].type != "Prime")  
    { //We are a module - need to assign type (if we are not a merged node whichs type already has been assigned in creation
      assignModuleType(rootG, graph, *ti, fac, lcutters, rcutters, fractureTree);     
    }
    //printTree(fractureTree);
  }
}

// Given the factorizing permutation calculates the modular decomposition
Tree ModDecomp::calcModDecomp(const std::vector<Vertex>& factorization, Graph& graph)
{ 
  const auto& facPos = boost::get(facPos_t(), graph);
  int ct = 0;
  for (const auto& v : factorization)
    boost::put(facPos, v, ct++);

  std::map<Vertex, int> lcutters;
  std::map<Vertex, int> rcutters; // rcutter of vertex
  // Check for connected components first
  DyckWord dyck = parenthesizing(factorization, graph, lcutters, rcutters); // cutters get set here
//   std::cout << dyck.toString(graph) << std::endl;
  Tree fractureTree = buildFractureTree(dyck, graph, lcutters, rcutters);
//   printTree(fractureTree,true);
  moduleDetDel(fractureTree);
//   printTree(fractureTree);
  buildModDecomp(graph, factorization, lcutters, rcutters, fractureTree);
//   printTree(fractureTree);
  return fractureTree;
}

//DEBUG
void printPartition(const std::list<std::vector<Vertex> >& partition, const Graph& graph)
{
  const auto& name = boost::get(boost::vertex_name, graph);
  for (const auto& vec : partition)
  {
    for (const auto& v : vec)
    {
      std::cout << boost::get(name, v) << " " << std::ends;
    }
    std::cout << "/ " << std::ends;
  }
  std::cout << std::endl;
}

// function initPartition from paper
void initPartition(const Graph& graph, std::list<std::vector<Vertex> >& partition, std::list<std::vector<Vertex> >& modules, std::list<std::vector<Vertex> >& pivots, Vertex& center, std::map<std::vector<Vertex>, Vertex>& fPivots)
{
  bool singleton = true;
  std::list<std::vector<Vertex> >::iterator classPos;
  for (std::list<std::vector<Vertex> >::iterator it = partition.begin(); it != partition.end(); ++it)
  {
    if (it->size() > 1)
    {
      singleton = false;
      classPos = it;
      break;
    }
  }

  if (singleton) // propagate this outwards [we are finished]
    throw "Solution found";
  
  Vertex fPivot;
  if (modules.empty())
  {
    try
    {
      fPivot = fPivots.at(*classPos);
    } catch (...) // no pivots set
    {
      fPivot = (*classPos)[0]; // choose pivot as first vertex
    }
    std::vector<Vertex> nonNeighbours;
    std::vector<Vertex> pivot{fPivot};
    std::vector<Vertex> neighbours;
    aiter ai, ai_end;
    boost::tie(ai,ai_end) = boost::adjacent_vertices(fPivot,graph);
    std::vector<Vertex> adj;
    for (; ai != ai_end; ++ai)
      adj.push_back(*ai);
    for (std::vector<Vertex>::iterator it = classPos->begin(); it != classPos->end(); ++it)
    {
      if (*it != fPivot)
      {
	if (std::find(adj.begin(), adj.end(), *it) != adj.end())
	  neighbours.push_back(*it);
	else
	  nonNeighbours.push_back(*it);
      }
    }
    
    // if pivot is in a "real" module, one of neighbours or nonNeighbours will definitely be empty, do not add empty sets to the partition please.
    if (!nonNeighbours.empty())
      partition.insert(classPos,nonNeighbours);
    partition.insert(classPos,pivot); // center
    if (!neighbours.empty())
      partition.insert(classPos,neighbours);
    partition.erase(classPos); // classPos is gone now
    center = fPivot;
    
    // And no empty sets to pivots nor to modules
    if (neighbours.size() < nonNeighbours.size())
    {
      if (!neighbours.empty())
      {
	/*std::cout << "added neighbors as pivot: " << std::ends;
	for (const auto& v : neighbours )
	  std::cout << v << " " << std::ends;
	std::cout << std::endl;*/
	pivots.push_back(neighbours);
      }
      modules.insert(modules.begin(), nonNeighbours); // begin, so modPiv is consistent in any case
    }
    else if (neighbours.size() >= nonNeighbours.size() && !neighbours.empty())
    {
      if (!nonNeighbours.empty())
      {
	/*std::cout << "added non-neighbors as pivot: " << std::ends;
	for (const auto& v : nonNeighbours )
	  std::cout << v << " " << std::ends;
	std::cout << std::endl;*/
	pivots.push_back(nonNeighbours);
      }
      modules.insert(modules.begin(), neighbours);
    }
    // else case: original class was singleton -> contradiction
  }
  else
  {
    std::vector<Vertex> x = *(modules.begin());
//     pivots.push_back(std::vector<Vertex>{x[0]});
    /*std::cout << "added first module as pivot: " << std::ends;
    for (const auto& v : x )
      std::cout << v << " " << std::ends;
    std::cout << std::endl;*/
    pivots.push_back(x);
    fPivots.insert(std::make_pair(*(modules.begin()),x[0]));
    modules.erase(modules.begin());
  }
}

// create the pivot set for vertex y from set x
std::vector<Vertex> pivotSet(Vertex y, const std::vector<Vertex>& x, const Graph& g)
{
  //return N(y) without x
  std::vector<Vertex> v;
  aiter ai, ai_end;
  for (boost::tie(ai,ai_end) = boost::adjacent_vertices(y,g); ai != ai_end; ++ai)
  {
    if (std::find(x.begin(), x.end(), *ai) == x.end())
      v.push_back(*ai);
  }
  return v;  
}

template <typename T>
bool isSubset(std::vector<T> B, std::vector<T> A)
{
    std::sort(A.begin(), A.end());
    std::sort(B.begin(), B.end());
    return std::includes(A.begin(), A.end(), B.begin(), B.end());
}

bool insertRight(const std::vector<Vertex>& x, const std::vector<Vertex>& y, const std::list<std::vector<Vertex> >& partition, const Vertex& center)
{
  bool foundX = false;
  bool foundY = false;
  bool foundCenter = false;
  for (const auto& vec : partition)
  {
    if (vec[0] == center)
    {
      foundCenter = true;
      if (foundX and foundY)
	return false;
      else if (foundX or foundY)
	return true;
    }
    else // check whether (*it) is Y or X
    {
      /*for (const auto& v : vec)
	std::cout << v << " ";
      std::cout << "is " << (isSubset(vec,x) ? "subset" : "not subset") << " of ";
      for (const auto& v : x)
	std::cout << v << " ";
      std::cout << std::endl;
      std::cout << "and " << std::ends;
      std::cout << "is " << (isSubset(vec,y) ? "subset" : "not subset") << " of ";
      for (const auto& v : y)
	std::cout << v << " ";
      std::cout << std::endl;*/
      if (isSubset(vec,x))
	foundX = true;
      else if (isSubset(vec,y))
      {
	foundY = true;
	if (foundCenter and foundX)
	  return false;
	else if (foundCenter or foundX)
	  return true;
      }
    }     
  }
  return true; //Y has not been found
}

// after the set x has been split to x and x_a, add pivots accordingly
void addPivot(const std::vector<Vertex>& x, const std::vector<Vertex>& x_a, std::list<std::vector<Vertex> >& modules, std::list<std::vector<Vertex> >& pivots)
{
  if (std::find(pivots.begin(), pivots.end(), x) != pivots.end() or x_a.size() < x.size())
  {
    /*std::cout << "added x_a as pivot: " << std::ends;
    for (const auto& v : x_a )
      std::cout << v << " " << std::ends;
    std::cout << std::endl;*/
    pivots.push_back(x_a);
    std::list<std::vector<Vertex> >::iterator it = std::find(modules.begin(), modules.end(), x);
    if (it == modules.end())
      modules.push_back(x);
  }
  else
  {
    /*std::cout << "added x as pivot: " << std::ends;
    for (const auto& v : x )
      std::cout << v << " " << std::ends;
    std::cout << std::endl;*/
    pivots.push_back(x);
    std::list<std::vector<Vertex> >::iterator it = std::find(modules.begin(), modules.end(), x);
    if (it != modules.end())
    {
      modules.insert(it, x_a);
      modules.erase(it);
    }
    else
      modules.push_back(x_a);
  }
}

// as described in the paper
void refine(const std::vector<Vertex>& pivotSet, const std::vector<Vertex>& Y, std::list<std::vector<Vertex> >& partition, std::list<std::vector<Vertex> >& modules, std::list<std::vector<Vertex> >& pivots, Vertex& center, const Graph& graph)
{  
  std::map<std::vector<Vertex>,std::vector<Vertex> > overlap;
  // get the "properly overlapping modules"
  for (const auto& v : pivotSet)
  {
    for (const auto& vec : partition)
    {
      if (std::find(vec.begin(), vec.end(), v) != vec.end())
      {
	overlap[vec].push_back(v);
	break; // we can only be in one part of the partition [partition - duh]
      }
    }
  }
  for (std::list<std::vector<Vertex> >::iterator vec = partition.begin(); vec != partition.end(); ++vec)
  {
    std::vector<Vertex> toMap = *vec; // save the position we're working on
    
    std::vector<Vertex>::iterator pos = std::remove_if(vec->begin(), vec->end(), [&pivotSet](Vertex elem) 
    {
      return (std::find(pivotSet.begin(), pivotSet.end(), elem) != pivotSet.end());   
    }); //moves the elements to delete after pos
    unsigned int numRemove = std::distance(pos, vec->end());
    if (numRemove > 0 && numRemove < vec->size()) // more than 1 element would be removed and proper subset (not all elements are removed)
    {
      vec->erase(pos, vec->end()); // delete the previous calculated elements of x_a
      std::vector<Vertex> ins = *vec;
      // add x_a to partition
      if (insertRight(*vec, Y, partition, center))
      {
	/*std::cout << "inserted ";
	for (const auto& v : overlap[toMap])
	  std::cout << v << " ";
	std::cout << "right of ";
	for (const auto& v : *vec)
	  std::cout << v << " ";
	std::cout << std::endl;*/
	partition.insert(++vec, overlap[toMap]); // toMap because we erased from vec such that the map for vec is now empty
	--vec; // or do we
      }
	// it is okay to just use "++vec" because vec now points to the newly added element, for which numRemove always would be = vec->size (they obv. fully overlap)
      else
      {
	/*std::cout << "inserted ";
	for (const auto& v : overlap[toMap])
	  std::cout << v << " ";
	std::cout << "left of ";
	for (const auto& v : *vec)
	  std::cout << v << " ";
	std::cout << std::endl;*/
	partition.insert(vec, overlap[toMap]);
      }
      addPivot(ins, overlap[toMap], modules, pivots);
    }
  }
}

// core function
void partitionRefinement(std::list<std::vector<Vertex> >& partition, const Graph& graph)
{
  std::list<std::vector<Vertex> > pivots;
  std::list<std::vector<Vertex> > modules;
  Vertex center;
  std::map<std::vector<Vertex>,Vertex> fPivots;
  while (true) // repeat until initPartition throws - this means partition of singletons
  {
    try
    {
//       std::cout << "init" << std::endl;
//       printPartition(partition, graph);
      initPartition(graph, partition, modules, pivots, center, fPivots);
//       printPartition(partition,graph);
    } catch (...) // catch specific?
    {
//       printPartition(partition, graph);
      return; // if the error occured, partition is already a set of singletons
    }
    while (!pivots.empty())
    {
      std::vector<Vertex> curr = pivots.front(); // get first pivot from pivots
      pivots.pop_front();
      for (const auto& v : curr)
      {
	std::vector<Vertex> p = pivotSet(v,curr,graph);
// 	std::cout << "refine with " << v << std::endl;/*<< " of pivot set: " << std::ends*/;
// 	for (const auto& v_ : curr)
// 	  std::cout << v_ << " " << std::ends;
// 	std::cout << "having the neighbours: " << std::ends;
// 	for (const auto& v_ : p)
// 	  std::cout << v_ << " " << std::ends;
// 	std::cout << std::endl;
// 	printPartition(partition, graph);
	refine(p, curr, partition, modules, pivots, center, graph);
// 	printPartition(partition, graph);
      }
    }
  }
}

// given a graph calculates ONE possible factorizing permutation (there might be as many as N! different ones)
std::vector<Vertex> ModDecomp::calcFacPerm(const Graph& graph)
{
  vertex_iter vi, vi_end;
  boost::tie(vi,vi_end) = boost::vertices(graph);
  std::vector<Vertex> vertices;
  for (; vi != vi_end; ++vi)
    vertices.push_back(*vi);
  std::list<std::vector<Vertex> > partition;
  partition.push_back(vertices);
  partitionRefinement(partition, graph);
  std::vector<Vertex> facPerm;
  for (const auto& vec : partition)
    facPerm.push_back(vec[0]); // should be singletons
  return facPerm;
}

Tree ModDecomp::decompose(Graph& g)
{
  Graph h;
  const auto& name = boost::get(boost::vertex_name, g);
  
  vertex_iter vi, vi_end;
  int ct = 0;
  for (boost::tie(vi,vi_end) = boost::vertices(g); vi != vi_end; ++vi, ct++)
  {
    std::string vname = "v" + std::to_string(ct);
    boost::put(name, *vi, vname);
  }
  std::vector<Vertex> facPerm = calcFacPerm(g);
  Tree t = calcModDecomp(facPerm, g);
  
  return t;
}
