#include "ModDecomp.h"

std::map<std::string, TVertex> namelist_c;
std::map<std::string, Vertex> namelist_g;
std::map<Vertex, Vertex> lcutters; // lcutter of vertex
std::map<Vertex, Vertex> rcutters; // rcutter of vertex
std::map<Vertex, unsigned int> facPos; // position in factorizing permutation
Tree fractureTree;

void setName(vertex_name& node, std::string name, Vertex vertex)
{
	if(namelist_g.insert(std::pair<std::string, Vertex>(name,vertex)).second)
	{
		node.name = name;
	}
	else
	{
		std::cerr << "Duplicate name" << std::endl;
	}
}

void initNode(vertex_prop& node, std::string name, TVertex& vertex)
{
	node.lfrac = std::numeric_limits<int>::max();
	node.rfrac = -1;
	node.toDelete = false;
	node.leaf = false;
	node.root = false;
	if(namelist_c.insert(std::pair<std::string, TVertex>(name,vertex)).second)
	{
		node.name = name;
		if (name != "Internal")
		  node.leaf = true;
	}
	else // Internal nodes in fracture Tree are label as you please.
	{
		node.name = name;
	}
}

/* Given a factorizing permutation constructs the dyck word (parenthesized factorization) depending on the graph */
std::vector<std::string> parenthesizing(std::vector<std::string> fac, Graph graph) // make sure the graph names are the names appearing in the string
{
  // we count for every string opening and closing brackets before/after vertices to build the dyck word from it later on
  std::map<Vertex, unsigned int> lbrackets;
  std::map<Vertex, unsigned int> rbrackets;
  
  std::vector<std::string> dyck_word;
  int pos = 0;
  
  
  for (std::vector<std::string>::iterator siter = fac.begin(); siter != fac.end() - 1; ++siter) // -1 because we will be looking for pairs
  {
    // We should be able to read the pos from the factorizing permutation
    /*facPos.insert(std::make_pair(namelist_g[*siter],pos++));
    if (siter == fac.end() - 1)
      facPos.insert(std::make_pair(namelist_g[*(siter + 1)],pos)); // before finishing: set the pos of the last vertex in the factorizing permutation*/
      
      
      
    // The two adjacency lists as set
    std::vector<std::string> nh1;
    for (boost::tie(ai, ai_end) = boost::adjacent_vertices(namelist_g[*siter], graph); ai != ai_end; ++ai)
      nh1.push_back(graph[*ai].name);
      
    std::vector<std::string> nh2;
    for (boost::tie(ai, ai_end) = boost::adjacent_vertices(namelist_g[*(siter + 1)], graph); ai != ai_end; ++ai)
      nh2.push_back(graph[*ai].name);
    
    std::vector<std::string> cutset; 
    /*std::set_symmetric_difference(nh1.begin(), nh1.end(), 
				  nh2.begin(), nh2.end(), 
				  cutset.begin(), std::strcmp);*/ 
    // store xor in cutset - these are the nodes for which the pair does not agree upon
    // Okay, this is possibly inefficient:
    for (std::vector<std::string>::iterator iter = nh1.begin(); iter != nh1.end(); ++iter)
      if (*iter != *(siter + 1)) // this is the other node of the pair
	if (std::find(nh2.begin(), nh2.end(), *iter) == nh2.end()) // not contained in its neighbours
	  cutset.push_back(*iter);
      
    for (std::vector<std::string>::iterator iter = nh2.begin(); iter != nh2.end(); ++iter)
      if (*iter != *(siter))
	if (std::find(nh1.begin(), nh1.end(), *iter) == nh1.end())
	  cutset.push_back(*iter);
	
    //cutset is not necessarily sorted now
      
    // left fractures
    for (std::vector<std::string>::iterator iter = fac.begin(); iter != siter; ++iter) // iter != siter: we cannot be a left fracture if right of pair
    {
      if (std::find(cutset.begin(), cutset.end(), *iter) != cutset.end())
      {
	//left fracture found
	//insert closing bracket right of siter
	if (rbrackets.find(namelist_g[*siter]) != rbrackets.end())
	  rbrackets[namelist_g[*siter]]++;
	else
	  rbrackets[namelist_g[*siter]] = 1;
	//insert opening bracket left of the vertex corresponding to iter
	if (lbrackets.find(namelist_g[*iter]) != lbrackets.end())
	  lbrackets[namelist_g[*iter]]++;
	else
	  lbrackets[namelist_g[*iter]] = 1;
	lcutters.insert(std::make_pair(namelist_g[*siter], namelist_g[*(std::find(cutset.begin(), cutset.end(), *iter))]));
	break; // left fracture found, cutters right from this are ignored
      }
    }
    
    // right fractures
    for (std::vector<std::string>::iterator iter = fac.end() - 1; iter != siter; iter--)
    {
      if (std::find(cutset.begin(), cutset.end(), *iter) !=  cutset.end())
      {
	//right fracture found
	//insert opening bracket left of (siter + 1)
	if (lbrackets.find(namelist_g[*(siter+1)]) != lbrackets.end())
	  lbrackets[namelist_g[*(siter + 1)]]++;
	else
	  lbrackets[namelist_g[*(siter + 1)]] = 1;
	//insert closing bracket right of the vertex corresponding to iter
	if(rbrackets.find(namelist_g[*iter]) != rbrackets.end())
	  rbrackets[namelist_g[*iter]]++;
	else
	  rbrackets[namelist_g[*iter]] = 1;
	rcutters.insert(std::make_pair(namelist_g[*siter], namelist_g[*(std::find(cutset.begin(), cutset.end(), *iter))]));
	break; // right fracture found, cutter left from this are ignored
      }
    }
  }
   //Now build the dyck word
  for (std::vector<std::string>::iterator iter = fac.begin(); iter != fac.end(); ++iter)
  {
    unsigned int brackets = 0;
    if (lbrackets.find(namelist_g[*iter]) != lbrackets.end()) // add all opening brackets before vertex
      brackets = lbrackets[namelist_g[*iter]];
    for (unsigned int i = 0; i < brackets; i++)
      dyck_word.push_back("(");
    
    dyck_word.push_back(*iter); // vertex
    
    brackets = 0;
    if (rbrackets.find(namelist_g[*iter]) != rbrackets.end()) // add all closing brackets after vertex
      brackets = rbrackets[namelist_g[*iter]];
    for (unsigned int i = 0; i < brackets; i++)
      dyck_word.push_back(")");
  }
  
  return dyck_word;
}

/* Given the Dyck word builds the fracture tree */
void buildFractureTree(std::vector<std::string> pfac)
{
  TVertex root = boost::add_vertex(fractureTree);
  initNode(fractureTree[root], "Internal", root); // we could also call this "Internal"? (root?)
  fractureTree[root].root = true;
  TVertex curr = root;
  
  Vertex lastChild;
  
  for (std::vector<std::string>::iterator iter = pfac.begin(); iter != pfac.end(); ++iter){
    if (*iter == "(") // opening bracket: add vertex as child of curr
    {
      TVertex v = boost::add_vertex(fractureTree);
      initNode(fractureTree[v], "Internal", v); // These will later on become the module nodes [series/parallel/prime]
      boost::add_edge(curr,v,fractureTree);
      curr = v; // new vertex becomes curr
    }
    else if (*iter == ")") // closing bracket: backtrack to father of current node
    {     
      if (boost::out_degree(curr,fractureTree) == 1) // we only have one child: dummy node 
	fractureTree[curr].toDelete = true;
      unsigned int ltemp = fractureTree[curr].lfrac;
      unsigned int rtemp = fractureTree[curr].rfrac;
      std::vector<std::string> nodes = fractureTree[curr].containedNodes;
      curr = boost::source(*(boost::in_edges(curr,fractureTree).first),fractureTree); // no null-check. Assume dyck word is well-formed
      // we leave this node for good, give its cutters to parent node
      fractureTree[curr].containedNodes.insert(fractureTree[curr].containedNodes.begin(), nodes.begin(), nodes.end());
      fractureTree[curr].lfrac = std::min(fractureTree[curr].lfrac, ltemp);
      fractureTree[curr].rfrac = std::min(fractureTree[curr].rfrac, rtemp);
      if (fractureTree[curr].rfrac > facPos[lastChild] && fractureTree[curr].name == "Internal") /* what about lfrac?*/
	fractureTree[curr].toDelete = true; // this node is a dummy node
    }
    else { // otherwise [and we do not check further, but this _should_ be a node name from the factorization] add "true" vertex as leaf
      TVertex leaf = boost::add_vertex(fractureTree);
      initNode(fractureTree[leaf], *iter, leaf);
      
      fractureTree[leaf].leaf = true; // This should not be needed
      fractureTree[leaf].type = "leaf";
      
      fractureTree[curr].containedNodes.push_back(*iter); //we might need this for module detection
      
      boost::add_edge(curr,leaf,fractureTree);
      //
      if (boost::out_degree(curr,fractureTree) == 1) // we're the first node added
      {
	fractureTree[curr].lfrac = facPos[namelist_g[*iter]]; // neat /* lcutters[namelist_fractureTree[*iter]].pos */
	if (rcutters.find(namelist_g[*iter]) != rcutters.end())
	  fractureTree[curr].rfrac = facPos[namelist_g[*iter]]; // set right fracture. pos? really?
      }
      else 
      {
	fractureTree[curr].rfrac = std::max(fractureTree[curr].rfrac, facPos[rcutters[namelist_g[*iter]]]); // rightmost cutter
      }
      lastChild = namelist_g[*iter]; // this is the last processed child
    }
    
    if (boost::out_degree(root,fractureTree) == 1)
      fractureTree[root].toDelete = true;
    
  }
  
}


/* Detects all modules and deletes dummy nodes which do not represent a (possibly weak) module */
void moduleDetDel() //Detect modules, delete dummies 
{
  // Delete all nodes which flags have been set in the previous step
  tvertex_iter vi, vi_end, next;
  boost::tie(vi,vi_end) = boost::vertices(fractureTree);
  for (next = vi; vi != vi_end; vi = next)
  {
    ++next;
    if (fractureTree[*vi].toDelete && !fractureTree[*vi].root) // we cant do this if we are root (link to parents)
    {     
      edge_iter ei,ei_end, nexte;
      boost::tie(ei,ei_end) = boost::out_edges(*vi, fractureTree);
      for (nexte = ei ; ei != ei_end; ei = nexte) // remove all edges from node
      {
	++nexte;
	boost::add_edge(boost::source(*(boost::in_edges(*vi,fractureTree).first),fractureTree),boost::target(*ei,fractureTree),fractureTree); // link children to parent of vertex to be removed
	boost::remove_edge(*ei,fractureTree);
      }
      boost::remove_edge(boost::source(*(boost::in_edges(*vi,fractureTree).first),fractureTree),*vi,fractureTree); // remove edge to node
      boost::remove_vertex(*vi,fractureTree); // all edges are already gone
    }
    else if (fractureTree[*vi].toDelete && fractureTree[*vi].root) // this is only the case iff root has exactly one child. Otherwise the root is always a module: The one containing all nodes
    {
      fractureTree[boost::target(*(boost::out_edges(*vi,fractureTree).first),fractureTree)].root = true; // we can do this
      boost::remove_edge(*(boost::out_edges(*vi,fractureTree).first),fractureTree); // remove _the_ out edge (can only be one)
      // we dont link child to parent because root does not have parent
      boost::remove_vertex(*vi,fractureTree); // there are no in-edges so we can delete the old root now
    }
  }
}

/* Finally builds the modular decomposition with the reduced fracture tree and the graph structure */
void buildModDecomp(Graph &graph) //we need representative graphs for module merging
{  
  //SubGraph tmp;
  //boost::copy_graph(graph, tmp);
  //boost::property_map<SubGraph, boost::vertex_name_t>::type name = get(boost::vertex_name, tmp);
}

Graph& calcModDecomp(std::vector<std::string> factorization, Graph &graph)
{
  // We should not need this anymore
  /*std::pair<vertex_iter, vertex_iter> vp;
  for (vp = boost::vertices(graph); vp.first != vp.second; ++vp.first)
  {
    namelist_g.insert( std::make_pair(graph[*(vp.first)].name,*vp.first)); // This should be consistent when the graph was generated correctly
  }*/
  
  std::vector<std::string> dyck = parenthesizing(factorization,graph);
  buildFractureTree(dyck);
  moduleDetDel();
  buildModDecomp(graph);
  vertex_iter vi, vi_end;
  for (boost::tie(vi,vi_end) = boost::vertices(fractureTree); vi != vi_end; ++vi)
    std::cout << "Node: " << fractureTree[*vi].name << " of type " << fractureTree[*vi].type << std::endl;
  boost::write_graphviz(std::cout, fractureTree); 
  
  return graph;
  
}

int main()
{
  Graph g;
  boost::property_map<Graph, boost::vertex_name_t>::type name = boost::get(boost::vertex_name, g);
  boost::add_edge(1,2,g);
  boost::add_edge(1,3,g);
  boost::add_edge(1,4,g);
  boost::add_edge(2,4,g);
  boost::add_edge(2,5,g);
  boost::add_edge(2,6,g);
  boost::add_edge(2,7,g);
  boost::add_edge(3,4,g);
  boost::add_edge(3,5,g);
  boost::add_edge(3,6,g);
  boost::add_edge(3,7,g);
  boost::add_edge(4,5,g);
  boost::add_edge(4,6,g);
  boost::add_edge(4,7,g);
  boost::add_edge(5,6,g);
  boost::add_edge(5,7,g);
  boost::add_edge(6,8,g);
  boost::add_edge(6,9,g);
  boost::add_edge(6,10,g);
  boost::add_edge(6,11,g);
  boost::add_edge(7,8,g);
  boost::add_edge(7,9,g);
  boost::add_edge(7,10,g);
  boost::add_edge(7,11,g);
  boost::add_edge(8,9,g);
  boost::add_edge(8,10,g);
  boost::add_edge(8,11,g);
  boost::add_edge(9,10,g);
  boost::add_edge(9,11,g);
  vertex_iter vi, vi_end;
  int ct = 1;
  for (boost::tie(vi,vi_end) = boost::vertices(g),; vi != vi_end; ++vi, ct++)
    boost::put(name, *vi, std::to_string(ct));
  /*Vertex v1 = boost::add_vertex(g);
  setName(g[v1],"1",v1);
  Vertex v2 = boost::add_vertex(g);
  setName(g[v2],"2",v2);
  Vertex v3 = boost::add_vertex(g);
  setName(g[v3],"3",v3);
  Vertex v4 = boost::add_vertex(g);
  setName(g[v4],"4",v4);
  Vertex v5 = boost::add_vertex(g);
  setName(g[v5],"5",v5);
  Vertex v6 = boost::add_vertex(g);
  setName(g[v6],"6",v6);
  Vertex v7 = boost::add_vertex(g);
  setName(g[v7],"7",v7);
  Vertex v8 = boost::add_vertex(g);
  setName(g[v8],"8",v8);
  Vertex v9 = boost::add_vertex(g);
  setName(g[v9],"9",v9);
  Vertex v10 = boost::add_vertex(g);
  setName(g[v10],"10",v10);
  Vertex v11 = boost::add_vertex(g);
  setName(g[v11],"11",v11);
  boost::add_edge(v1,v2,g);
  boost::add_edge(v1,v3,g);
  boost::add_edge(v1,v4,g);
  boost::add_edge(v2,v4,g);
  boost::add_edge(v2,v5,g);
  boost::add_edge(v2,v6,g);
  boost::add_edge(v2,v7,g);
  boost::add_edge(v3,v4,g);
  boost::add_edge(v3,v5,g);
  boost::add_edge(v3,v6,g);
  boost::add_edge(v3,v7,g);
  boost::add_edge(v4,v5,g);
  boost::add_edge(v4,v6,g);
  boost::add_edge(v4,v7,g);
  boost::add_edge(v5,v6,g);
  boost::add_edge(v5,v7,g);
  boost::add_edge(v6,v8,g);
  boost::add_edge(v6,v9,g);
  boost::add_edge(v6,v10,g);
  boost::add_edge(v6,v11,g);
  boost::add_edge(v7,v8,g);
  boost::add_edge(v7,v9,g);
  boost::add_edge(v7,v10,g);
  boost::add_edge(v7,v11,g);
  boost::add_edge(v8,v9,g);
  boost::add_edge(v8,v10,g);
  boost::add_edge(v8,v11,g);
  boost::add_edge(v9,v10,g);
  boost::add_edge(v9,v11,g);*/
  
  std::vector<std::string> par;
  for (unsigned int i = 1; i < 12; i++)
    par.push_back(std::to_string(i));
  
  calcModDecomp(par,g);
  
  return 0;
}
 
