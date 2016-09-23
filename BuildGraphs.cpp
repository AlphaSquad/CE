#include "BuildGraphs.h"

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

BuildGraphs::BuildGraphs(unsigned int seed, bool print)
{
  if (seed == 0)
  {
    std::random_device rd;
    seed = rd();
  }
  if (print)
    std::cout << "Using random seed " << seed << std::endl;
  generator.seed(seed);
}

void BuildGraphs::ReadEdgeList(std::string sFileName, Graph& g, bool& FileExists)
{	
  int edgeID=0;
    std::ifstream fileIn;
    fileIn.open(sFileName.c_str(), std::ifstream::in);

    if(!fileIn.is_open()){
		std::cout << "!error - no file: "<< sFileName << std::endl;
		FileExists = false;
    }
    else{
		std::istream_iterator<int> begin(fileIn), end;
		int      x, y;
		while(begin != end){
			x = *begin;
			y = *(++begin);
			boost::add_edge(x, y, edgeID, g);
			++edgeID;
			begin++;
		}
		FileExists = true;
		fileIn.close();
    }
}

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

std::vector<Graph> BuildGraphs::readGraphs(std::string filename)
{
  std::ifstream ifs;
  bool flag = false;  
  
  ifs.open(filename);
  std::string line;
  std::vector<Graph> graphs;
  while (std::getline(ifs,line))
  {
    if (line == "#OrthologyMatrix" or (line == "#Cograph" and !flag)) // now follow the graphs
    {
      Graph g;
      int i = 0;
      int ct = 0;
      while (std::getline(ifs,line))
      {
	if (line.length() == 0) // This graph is finished
	  break;
	std::vector<std::string> l = split(line, '\t');
	int j = 0;
	for (const auto& w : l)
	{
	  if (w == "1" and j > i)
	    boost::add_edge(i,j,ct++,g);
	  j++;
	}
	i++;
      }
      flag = true;
//       boost::write_graphviz(std::cout, g);
      graphs.push_back(g);
    }
    if (line == "#GeneSpeciesAssociations") // So multiple cographs are not read multiple times
      flag = false;
  }
  ifs.close();
  return graphs;
}

int BuildGraphs::rootedTree(int toAdd, Cotree& t, CVertex& v, bool series, int num)
{
  const auto& name = boost::get(boost::vertex_name, t);
  const auto& root = boost::get(root_t(), t);
  
  bool second = false;
  while (toAdd > 0)
  {
    if (toAdd == 1)
    {
      CVertex leaf = boost::add_vertex(t);
      std::string vname = "v" + std::to_string(++num);
      boost::put(name, leaf, vname);
      boost::put(root, leaf, false);
      boost::add_edge(v,leaf,t);
      return num;
    }
    
    std::uniform_real_distribution<double> distr (2.0,4.0);
    std::normal_distribution<double> distribution(double(toAdd)/distr(generator),2.0); //this is just random numbers. [we want 2.5 children on average, dont know what variance 2.0 means :D ]
    
    int subTree = 0;
    do {
      double child = distribution(generator);
      subTree = rint(child);
    } while ((!second and subTree > toAdd - 1) or subTree > toAdd or subTree < 1); // we must have at least 2 children (1 in toAdd left after first) and every child has to have at least 1 node
    
    if (subTree == 1)
    {
      CVertex leaf = boost::add_vertex(t);
      std::string vname = "v" + std::to_string(++num);
      boost::put(name, leaf, vname);
      boost::put(root, leaf, false);
      boost::add_edge(v,leaf,t);
      toAdd--;
      continue;
    }
    
    CVertex curr = boost::add_vertex(t);
    boost::put(root, curr, false);
    if (series)
      boost::put(name, curr, "Series");
    else
      boost::put(name, curr, "Parallel");
    
    boost::add_edge(v,curr,t);
    num = rootedTree(subTree,t,curr,!series,num); // recurr
    toAdd -= subTree;
    second = true;
  }
  return num;
}

CVertex getLCA (const CVertex& v1, const CVertex& v2, const Cotree& tree)
{
  // this is a trivial implementation but should be fast enough
  const auto& root = boost::get(root_t(), tree);
  std::vector<CVertex> ancestors; // stores ancestors of v1
  CVertex tmp = v1;
  while (!boost::get(root,tmp))
  {
    tmp = *boost::inv_adjacent_vertices(tmp, tree).first;
    ancestors.push_back(tmp);
  }
  tmp = v2;
  while (true)
  {
    const auto& pos = std::find(ancestors.begin(), ancestors.end(), tmp);
    if (pos == ancestors.end())
      tmp = *boost::inv_adjacent_vertices(tmp, tree).first;
    else
      return (*pos);
  }
}

//DEBUG
void printTree(const Cotree& fractureTree)
{
  typedef std::map<CVertex, int> IndexMap;
  IndexMap mapIndex;
  boost::associative_property_map<IndexMap> propmapIndex(mapIndex);
  tvertex_iter vi, vi_end;
  int i = 0; // index
  for (boost::tie(vi,vi_end) = boost::vertices(fractureTree); vi != vi_end; ++vi)
    boost::put(propmapIndex,*vi,i++);
  
  boost::write_graphviz(std::cout, fractureTree, boost::make_label_writer(boost::get(boost::vertex_name,fractureTree)), boost::default_writer(), boost::default_writer(), propmapIndex);
}

Graph BuildGraphs::graphFromTree(const Cotree& cotree)
{
  //printTree(cotree);
  const auto& cname = boost::get(boost::vertex_name, cotree);
  Graph g;
  const auto& name = boost::get(boost::vertex_name, g);
  std::map<Vertex, CVertex> vertices;
  for (const auto& v : boost::vertices(cotree))
    if (boost::get(cname,v) != "Series" and boost::get(cname,v) != "Parallel")
    {
      Vertex v_ = boost::add_vertex(g);
      boost::put(name, v_, boost::get(cname, v));
      vertices.insert(std::make_pair(v_, v));
    }
    
  for (const auto& v1 : vertices)
    for (const auto& v2 : vertices)
      if (v1 < v2)
	if (boost::get(cname, getLCA(v1.second,v2.second,cotree)) == "Series")
	  boost::add_edge(v1.first, v2.first, g);
  
  return g;
}

Cotree BuildGraphs::buildCotree (unsigned int vertices)
{
  Cotree t;
  const auto& name = boost::get(boost::vertex_name, t);
  const auto& root = boost::get(root_t(), t);
  CVertex r = boost::add_vertex(t);
  boost::put(name, r, "Series");
  boost::put(root, r, true);
  rootedTree(vertices, t, r, false, 0); // recursion start
  return t;
}

std::string BuildGraphs::getSubtree(const Tree& t, const TVertex& v)
{
  if (boost::get(&vertex_prop::name,t,v) != "Internal")
  {
    return boost::get(&vertex_prop::name,t,v);
  }
  else
  {
    std::string tmp = "(";
    for (const auto& w : boost::adjacent_vertices(v,t))
    {
      tmp.append(getSubtree(t,w));
      tmp.append(",");
    }
    std::string tmp2 = tmp.substr(0,tmp.length() - 1);
    tmp2.append(")");
    return tmp2;
  }
}

std::string BuildGraphs::getNewick(const Tree& t)
{
  std::string newick;
  for (const auto& v: boost::vertices(t))
  {
    if (boost::get(&vertex_prop::root,t,v))
    {
      newick = getSubtree(t,v);
      newick.append(";");
      return newick;
    }
  }
  std::cerr << "Root not found" << std::endl;
  newick.append(";");
  return newick;
}

void BuildGraphs::writeGraphs(const std::vector<Graph>& graphs, std::string filename, std::vector<std::pair<int,int> > props)
{
  std::ofstream os;
  os.open(filename);
  os << "#NoOfSpecies\n";
  os << boost::num_vertices(graphs[0]);
  os << "\n#SpeciesNames\n";
  for (const auto& v : boost::vertices(graphs[0]))
    os << v << " ";
  os << "\n\n";
  int ct = 1;
  int i = 0;
  bool p = props.size() > 0;
  for (const auto& g : graphs)
  {
    os << "#GeneFamily " << ct++ << "\n";
    os << "#GeneSpeciesAssociations\n";
    for (const auto& v : boost::vertices(g))
      os << v << " ";
    os << std::endl;
    if (p)
    {
    	auto prop = props[i++];
    	os << "#Edits " << prop.first << std::endl;
	os << "#Time " << prop.second << " ms" << std::endl;
    }
    os << "#OrthologyMatrix\n";
    for (const auto& v1 : boost::vertices(g))
    {
      for (const auto& v2 : boost::vertices(g))
      {
	if (boost::edge(v1,v2,g).second)
	  os << 1 << "\t";
	else
	  os << 0 << "\t";
      }
      os << "\n";
    }
    os << "\n";
  }
  os.close();
}

// only for graphs with correct mapping
int BuildGraphs::calculateEditDistance(const Graph& g1, const Graph& g2)
{
  int edits = 0;
  std::stringstream ss;
  for (const auto& v : boost::vertices(g1))
    for (const auto& w : boost::vertices(g1))
    {
      bool e1 = boost::edge(v,w,g1).second;
      bool e2 = boost::edge(v,w,g2).second;
      if (v < w and (e1 xor e2))
        edits++;
    }
  return edits;
}

void BuildGraphs::writeGraph(const Graph& g, std::string filename)
{
  std::cerr << "Not implemented, try writeGraphs" << std::endl;
  std::ofstream ss;
}

Graph BuildGraphs::buildCograph (unsigned int vertices)
{
  Cotree t = buildCotree(vertices);
  Graph g = graphFromTree(t);
  return g;
}

// calculates some graph properties. t should be the corresponding MDT to g and an empty vector of properties
void BuildGraphs::graphProperties(std::vector<float>& prop, const Tree& t, const Graph& g)
{
  float primect = 0;
  float intct = 0;
  float leafct = boost::num_vertices(g);
  float primeleafct = 0;
  float primewayct = 0;
  for (const auto& v : boost::vertices(t))
  {
    if (t[v].name == "Internal")
    {
      intct++;
      if (t[v].type == "Prime")
	primect++;
    }
    else
    {
      if (t[*boost::inv_adjacent_vertices(v,t).first].type == "Prime")
	primeleafct++;
      auto curr = v;
      while (!t[curr].root)
      {
	if (t[*boost::inv_adjacent_vertices(curr,t).first].type == "Prime")
	{
	  primewayct++;
	  break;
	}
	curr = *boost::inv_adjacent_vertices(curr,t).first;
      }
    }
  }
  prop.push_back(primect/intct);
  prop.push_back(primeleafct/leafct);
  prop.push_back(primewayct/leafct);
}

Graph BuildGraphs::removeEdges(Graph& graph, float pc)
{
  std::uniform_real_distribution<float> distr(0.,1.);
  std::vector<Edge> toDelete;
  for (const auto& e : boost::edges(graph))
  {
    float prob = distr(generator);
    if (prob < pc)
      toDelete.push_back(e);
  }
  for (const auto& e : toDelete)
    boost::remove_edge(e,graph);
  return graph;
}

Graph BuildGraphs::editEdges(Graph& graph, float pc)
{ 
  std::uniform_real_distribution<float> distr(0.,1.);
  std::vector<std::pair<Vertex, Vertex> > toAdd;
  std::vector<Edge> toDelete;
  for (const auto& v : boost::vertices(graph))
    for (const auto& w : boost::vertices(graph))
    {
      if (v < w)
      {
	float num = distr(generator);
	if(!boost::edge(v,w,graph).second) // consider each non-edge only once
	{
	  if (num < pc)
	    toAdd.push_back(std::make_pair(v,w));
	}
	else
	{
	  if (num < pc)
	    toDelete.push_back(boost::edge(v,w,graph).first);
	}
      }
    }
  std::cout << (toDelete.size() + toAdd.size()) << std::endl;
  for (const auto& e : toDelete)
  {
    std::cout << "Deleted: " << boost::source(e,graph) << " - " << boost::target(e,graph) << std::endl;
    boost::remove_edge(e,graph);
  }
  for (const auto& e : toAdd)
  {
    std::cout << "Added: " << e.first << " - " << e.second << std::endl;
    boost::add_edge(e.first, e.second, graph);
  }
  std::cout << std::endl;
  return graph;
}

Graph BuildGraphs::addEdges(Graph& graph, float pc)
{
  std::uniform_real_distribution<float> distr(0.,1.);
  std::vector<std::pair<Vertex, Vertex> > toAdd;
  for (const auto& v : boost::vertices(graph))
    for (const auto& w : boost::vertices(graph))
    {
      if (v < w and !boost::edge(v,w,graph).second) // consider each non-edge only once
      {
	float num = distr(generator);
	if (num < pc)
	  toAdd.push_back(std::make_pair(v,w));
      }
    }
  for (const auto& e : toAdd)
    boost::add_edge(e.first, e.second, graph);
  return graph;
}
