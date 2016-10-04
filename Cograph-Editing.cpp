#include "Cograph-Editing.h"
#include "BuildGraphs.h"

#include <boost/lexical_cast.hpp>

// so we can boost iterators as normal iterators (ranged-based for loops!)
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
void printModDecomp(const Tree& fractureTree, std::ostream& strm = std::cout, bool toDelete = false, bool corrV = false)
{
  typedef std::map<TVertex, int> IndexMap;
  IndexMap mapIndex;
  boost::associative_property_map<IndexMap> propmapIndex(mapIndex);
  tvertex_iter vi, vi_end;
  int i = 0; // index
  for (boost::tie(vi,vi_end) = boost::vertices(fractureTree); vi != vi_end; ++vi)
    boost::put(propmapIndex,*vi,i++);
  if (toDelete)
    boost::write_graphviz(strm, fractureTree, boost::make_label_writer(boost::get(&vertex_prop::toDelete,fractureTree)), boost::default_writer(), boost::default_writer(), propmapIndex);
  else if (corrV)
    boost::write_graphviz(strm, fractureTree, boost::make_label_writer(boost::get(&vertex_prop::corrV,fractureTree)), boost::default_writer(), boost::default_writer(), propmapIndex);
  else
    boost::write_graphviz(strm, fractureTree, boost::make_label_writer(boost::get(&vertex_prop::type,fractureTree)), boost::default_writer(), boost::default_writer(), propmapIndex);
}

Vertex getCorr(const TVertex& v, const Tree& t)
{
  return t[v].containedNodes[0];
}

int addEdge(const QVertex& v1, const QVertex& v2, const QGraph& rep, Graph& g)
{
//   std::cout << "Added: ";
  int i = 0;
  for (const auto& v : boost::get(cont_t(),rep,v1))
    for (const auto& w : boost::get(cont_t(), rep, v2))
    {
      i++;
      assert(v != w);
//       std::cout << v << " - " << w << std::endl;
      boost::add_edge(v,w,g);
    }
  
  return i;
}

int removeEdge(const QVertex& v1, const QVertex& v2, const QGraph& rep, Graph& g)
{
//   std::cout << "Removed: " ;
  int i = 0;
  for (const auto& v : boost::get(cont_t(),rep,v1))
    for (const auto& w : boost::get(cont_t(), rep, v2))
    {
      i++;
      assert(v != w);
//       std::cout << v << " - " << w << std::endl;
      boost::remove_edge(v,w,g);
    }
  return i;
}

int editEdge(const QVertex& v1, const QVertex& v2, const QGraph& rep, Graph& g)
{
  assert(v1 != v2);
  if (boost::edge(boost::get(cont_t(),rep,v1)[0],boost::get(cont_t(),rep,v2)[0],g).second)
    return removeEdge(v1,v2,rep,g);
  else
    return addEdge(v1,v2,rep,g);
}

QGraph repGraph (const Tree& t, const TVertex& mod, const Graph& g)
{
  const auto& name = boost::get(boost::vertex_name, g);
  QGraph h;
  const auto& nname = boost::get(boost::vertex_name, h);
  const auto& weight = boost::get(weight_t(), h);
  std::map<QVertex, TVertex> map;
  for (const auto& v : boost::adjacent_vertices(mod,t))
  {
    QVertex w = boost::add_vertex(h);
    std::string n = boost::get(name,getCorr(v,t));
    boost::put(nname, w, n);
    boost::put(weight, w, t[v].containedNodes.size()); // vertex "weight" is number of contained vertices in module
    boost::put(cont_t(), h, w, t[v].containedNodes);
    map.insert(std::make_pair(w,v));
  }
  for (const auto& v : boost::vertices(h))
    for (const auto& w : boost::vertices(h))
    {
      if (v >= w) continue; // do not add edges twice
      TVertex v1 = map[v];
      TVertex v2 = map[w];
      if (boost::edge(getCorr(v1,t),getCorr(v2,t),g).second)
      {
	QEdge e = boost::add_edge(v,w,h).first;
	boost::put(boost::edge_weight,h,e,0);
	boost::put(boost::edge_weight2,h,e,0); // initialize weight of edges with 0
      }
    }
  return h;
}

typedef std::pair<std::pair<QVertex,QVertex>,std::pair<QVertex,QVertex> > p4;

bool compareP4 (p4 p1, p4 p2)
{
  return ((p1.first == p2.first and p1.second == p2.second) or
	  (p1.first.first == p2.second.second and p1.first.second == p2.second.first and
	  p1.second.first == p2.first.second and p1.second.second == p2.first.first));
}

typedef std::pair<std::pair<QVertex,QVertex>,std::pair<QVertex,QVertex> > P4;

// finds all p4s given the representative graph, scores them in the original graph as score of edge
float findP4s (QGraph& rGraph, bool print = false)
{
  std::vector<P4> p4s;
  const auto& weight = boost::get(boost::edge_weight, rGraph);
//   std::map<QEdge,std::set<p4> > p4_map;
  for (const auto& e : boost::edges(rGraph))
  {
    std::vector<QVertex> e_v1;
    std::vector<QVertex> e_v2;
    std::vector<QVertex> n_12;
    QVertex source = boost::source(e,rGraph);
    QVertex target = boost::target(e,rGraph);
    for (const auto& v1 : boost::adjacent_vertices(source,rGraph))
    {
      if (v1 == target) continue;
      if (boost::edge(v1,target,rGraph).second)
	n_12.push_back(v1); // is neighboured to both
      else
        e_v1.push_back(v1);
    }
    for (const auto& v2: boost::adjacent_vertices(target,rGraph))
    {
      if (v2 != source and !boost::edge(v2,source,rGraph).second) // if there is an edge v2-source, then we added that vertex already in the previous loop
	e_v2.push_back(v2);
    }
    for (const auto& v1 : e_v1) 
    {
      auto&& first_edge = boost::edge(v1,source,rGraph);
      for (const auto& v2 : e_v2)
      {
	if (v1 == v2) continue;
	
	auto&& third_edge = boost::edge(target,v2,rGraph);
	
	if (!boost::edge(v1,v2,rGraph).second)  
        {
          // We will find these triples in the complement as well, increase by 0.5 s.t. we count each triple only once
//           std::vector<QVertex> path{v1,source,target,v2};
//           std::sort(path.begin(),path.end());
//           triples[std::make_pair(std::make_pair(path[0],path[1]),path[2])] += 0.5f;
//           triples[std::make_pair(std::make_pair(path[0],path[1]),path[3])] += 0.5f;
//           triples[std::make_pair(std::make_pair(path[0],path[2]),path[3])] += 0.5f;
//           triples[std::make_pair(std::make_pair(path[1],path[2]),path[3])] += 0.5f;
	  boost::put(weight,first_edge.first,boost::get(weight,first_edge.first) + (boost::get(weight_t(), rGraph, v1) * boost::get(weight_t(), rGraph, source) * boost::get(weight_t(), rGraph, target) * boost::get(weight_t(), rGraph, v2)));
	  boost::put(weight,e,boost::get(weight,e) + (boost::get(weight_t(), rGraph, v1) * boost::get(weight_t(), rGraph, source) * boost::get(weight_t(), rGraph, target) * boost::get(weight_t(), rGraph, v2)));
	  boost::put(weight,third_edge.first,boost::get(weight,third_edge.first) + (boost::get(weight_t(), rGraph, v1) * boost::get(weight_t(), rGraph, source) * boost::get(weight_t(), rGraph, target) * boost::get(weight_t(), rGraph, v2)));
          p4s.push_back(std::make_pair(std::make_pair(v1,source),std::make_pair(target,v2)));
          if (print)
          {
            for (const auto& x : boost::get(cont_t(),rGraph,v1))
              std::cout << x << " ";
            std::cout << "- ";
            for (const auto& x : boost::get(cont_t(),rGraph,source))
              std::cout << x << " ";
            std::cout << "- ";
            for (const auto& x : boost::get(cont_t(),rGraph,target))
              std::cout << x << " ";
            std::cout << "- ";
            for (const auto& x : boost::get(cont_t(),rGraph,v2))
              std::cout << x << " ";
            std::cout << std::endl;
          }
	}
      }
    }
  }
  return p4s.size();
}

std::pair<std::map<QVertex,QVertex>,QGraph> getComplement(const QGraph& g)
{
  QGraph h;
  std::map<QVertex,QVertex> m;
  for (const auto& v : boost::vertices(g))
  {
    QVertex w = boost::add_vertex(h);
    boost::put(boost::vertex_name, h, w, boost::get(boost::vertex_name, g, v)); // same name in complement and original
    boost::put(weight_t(),h,w,boost::get(weight_t(),g,v)); // vertex "weight" is number of contained vertices in module
    boost::put(cont_t(),h,w,boost::get(cont_t(),g,v));
    m.insert(std::make_pair(v,w));
    assert(v == w); // otherwise our mapping gets destroyed
  }
  for (const auto& v : boost::vertices(g))
  {
    for (const auto& w : boost::vertices(g))
    {
      if (v >= w) continue; // undirected graph
      if (!boost::edge(v,w,g).second)
      {
	QEdge e = boost::add_edge(m[v],m[w],h).first;
	boost::put(boost::edge_weight,h,e,0);
	boost::put(boost::edge_weight2,h,e,0); // remember initialization
      }
    }
  }
  return std::make_pair(m,h);
}

std::pair<float,float> scoreOperation(const QVertex& v1,const QVertex& v,const QVertex& w, const QGraph& quotient, const QGraph& complement,std::map<QVertex,QVertex>& m,std::vector<std::pair<QVertex, QVertex> >& op)
{
  // score addition w-v1 and deletion v-v1
  // addition w-v1 is score of edge w-v1 in complement
  // deletion v-v1 is score of edge v-v1 in quotient
  
  float d_score = (float)boost::get(boost::edge_weight,quotient,boost::edge(v,v1,quotient).first);
  float a_score = (float)boost::get(boost::edge_weight,complement,boost::edge(m[w],m[v1],complement).first); // i know those exist
  
  // minus the p4s going over all v,w,v1 (triple score)
//   std::vector<QVertex> triplet{v,w,v1};  // This is the base amount of p4s we destroy in any case (all going over v AND w)
//   std::sort(triplet.begin(), triplet.end());
//   d_score -= triples[std::make_pair(std::make_pair(triplet[0],triplet[1]),triplet[2])];
//   a_score -= triples[std::make_pair(std::make_pair(triplet[0],triplet[1]),triplet[2])];

  float div_d = (float)(boost::get(weight_t(),quotient,v) * (float)boost::get(weight_t(),quotient,v1));
  float div_a = (float)(boost::get(weight_t(),complement,m[w]) * (float)boost::get(weight_t(),complement,m[v1]));
  
//   if ((boost::get(cont_t(),quotient,v)[0] == 7 and boost::get(cont_t(),quotient,w)[0] == 8) or (boost::get(cont_t(),quotient,v)[0] == 8 and boost::get(cont_t(),quotient,w)[0] == 7))
//   {
//     std::cout << "delete: ";
//     for (const auto& v_ : boost::get(cont_t(), quotient,v))
//       std::cout << v_ << " "; 
//     std::cout << "- ";
//     for (const auto& v_ : boost::get(cont_t(), quotient,v1))
//       std::cout << v_ << " ";
//     std::cout << ": " << d_score << "/" << div_d << ", add: ";
//     for (const auto& v_ : boost::get(cont_t(), quotient,w))
//       std::cout<< v_ << " ";
//     std::cout << "- ";
//     for (const auto& v_ : boost::get(cont_t(), quotient,v1))
//       std::cout << v_ << " ";
//     std::cout << ": " << a_score << "/" << div_a << std::endl;
//   }
//   else if ((boost::get(cont_t(),quotient,v)[0] == 3 and boost::get(cont_t(),quotient,w)[0] == 8) or (boost::get(cont_t(),quotient,v)[0] == 8 and boost::get(cont_t(),quotient,w)[0] == 3))
//   {
//     std::cout << "3|8 delete: ";
//     for (const auto& v_ : boost::get(cont_t(), quotient,v))
//       std::cout << v_ << " "; 
//     std::cout << "- ";
//     for (const auto& v_ : boost::get(cont_t(), quotient,v1))
//       std::cout << v_ << " ";
//     std::cout << ": " << d_score << "/" << div_d << ", add: ";
//     for (const auto& v_ : boost::get(cont_t(), quotient,w))
//       std::cout<< v_ << " ";
//     std::cout << "- ";
//     for (const auto& v_ : boost::get(cont_t(), quotient,v1))
//       std::cout << v_ << " ";
//     std::cout << ": " << a_score << "/" << div_a << std::endl;
//   }
  
  
  if (d_score/div_d > a_score/div_a) // maximize number of p4s per edge edition
  {
//     std::cout << "Remove " << boost::get(cont_t(), quotient,v)[0] << "/" << boost::get(cont_t(), quotient,v1)[0] << " (" << d_score << "/" << div_d << ">" << a_score << "/" << div_a << ")" << std::endl;
    op.push_back(std::make_pair(v,v1));
    return std::make_pair(d_score,div_d);
  }
  else
  {
//     std::cout << "Add " << boost::get(cont_t(), quotient,w)[0] << "/" << boost::get(cont_t(), quotient,v1)[0] << " (" << a_score << "/" << div_a << ">" << d_score<< "/" << div_d << ")" << std::endl;
    op.push_back(std::make_pair(w,v1)); // either add w-v1 or delete v-v1
    return std::make_pair(a_score,div_a);
  }
}

// sets the symmetric difference in neighborhoods of v and w into diff1 and diff2
void getSymmetricDifference(const QVertex& v,const QVertex& w, const QGraph& quotient, std::vector<QVertex>& diff1, std::vector<QVertex>& diff2)
{
  std::vector<QVertex> n_v{boost::adjacent_vertices(v,quotient).first, boost::adjacent_vertices(v,quotient).second};
  std::vector<QVertex> n_w{boost::adjacent_vertices(w,quotient).first, boost::adjacent_vertices(w,quotient).second};
  std::sort(n_v.begin(), n_v.end());
  std::sort(n_w.begin(), n_w.end());
  // get the difference in the neighborhoods
  std::set_difference(n_v.begin(), n_v.end(), n_w.begin(), n_w.end(), std::back_inserter(diff1));
  std::set_difference(n_w.begin(), n_w.end(), n_v.begin(), n_v.end(), std::back_inserter(diff2));
}

float calcCase2(const QVertex& v, const QVertex& v1, const QVertex& v2, const QGraph& quotient)
{
  float double_p4s = 0.;
  std::vector<QVertex> tmp_1;
  std::vector<QVertex> tmp_2;
  getSymmetricDifference(v1, v2, quotient, tmp_1, tmp_2);
  // two cases: double remove if neighbor to p->first and one of the neighbors
  for (const auto& vi : tmp_2)
    if (boost::edge(v,vi,quotient).second) // v == p->first
      double_p4s += boost::get(weight_t(),quotient,v) * boost::get(weight_t(), quotient, v1) * boost::get(weight_t(), quotient, v2) * boost::get(weight_t(), quotient, vi); 
    // or no remove at all (-2) if neighbor to both neighbors and not p->first
  for (const auto& vi : boost::adjacent_vertices(v1,quotient))
    if (!boost::edge(v,vi,quotient).second and boost::edge(v2,vi,quotient).second)
      double_p4s += 2 * boost::get(weight_t(),quotient,v) * boost::get(weight_t(), quotient, v1) * boost::get(weight_t(), quotient, v2) * boost::get(weight_t(), quotient, vi); // number of mutual neighbors of v1, v2
  return double_p4s;
}

float calcCase1(const QVertex& v, const QVertex& v1, const QVertex& v2, const QGraph& quotient, bool edge)
{
  float double_p4s = 0.;
  std::vector<QVertex> tmp_1;
  std::vector<QVertex> tmp_2;
  getSymmetricDifference(v1, v2, quotient, tmp_1, tmp_2);
  for (const auto& vi : tmp_1)
    if (!boost::edge(v,vi,quotient).second xor edge)
      double_p4s += boost::get(weight_t(),quotient,v) * boost::get(weight_t(), quotient, v1) * boost::get(weight_t(), quotient, v2) * boost::get(weight_t(), quotient, vi);
  for (const auto& vi : tmp_2)
    if (!boost::edge(v,vi,quotient).second xor edge)
      double_p4s += boost::get(weight_t(),quotient,v) * boost::get(weight_t(), quotient, v1) * boost::get(weight_t(), quotient, v2) * boost::get(weight_t(), quotient, vi);
    return double_p4s;
}

float checkDoubleP4s(const QVertex& v, const QVertex& v1, const QVertex& v2, const QGraph& quotient)
{
  const auto& e1 = boost::edge(v,v1, quotient);
  const auto& e2 = boost::edge(v,v2, quotient);
  const auto& e3 = boost::edge(v1, v2, quotient);
  std::vector<QVertex> tmp_1;
  std::vector<QVertex> tmp_2;
  float double_p4s = 0.f;
  if ((e1.second and e2.second and !e3.second) or (e1.second xor e2.second and e3.second))
  {
    double_p4s += calcCase1(v,v1,v2,quotient, false);
  }
  else if (e1.second and !e2.second and !e3.second) // we add one and remove one, p4s go vi-x-v-vj with x neighbour of v,vi, not of vj
  {
    double_p4s += calcCase2(v,v1,v2, quotient); 
  }
  else if (!e1.second and e2.second and !e3.second) // we add one and remove one, p4s go vi-x-v-vj with x neighbour of v,vi, not of vj (second way)
  {
    double_p4s += calcCase2(v,v2,v1, quotient);
  }
  else if (!e1.second and !e2.second and e3.second)
  {
    double_p4s += calcCase1(v,v1,v2,quotient, true);
  }
  else
  {
    //!e1.second and !e2.second and !e3.second cannot be a p4 on 4 vertices (triangular edges missing already)
    //e1.second and e2.second and e3.second is same case
    double_p4s += 0;
  }
  return double_p4s;
}

float checkDoubleP4s_case2(const QVertex& v, const QVertex& w, const QVertex& v1, const QVertex& v2, const QGraph& quotient)
{
  const auto& e1 = boost::edge(v,v1, quotient);
  const auto& e2 = boost::edge(v,v2, quotient);
  const auto& e3 = boost::edge(w,v1, quotient);
  const auto& e4 = boost::edge(w,v2, quotient);
  const auto& e5 = boost::edge(v1, v2, quotient);
  const auto& e6 = boost::edge(v,w,quotient);
  if (e1.second and e4.second)
    return (float)(((e5.second xor e6.second) and !(e2.second or e3.second)) or ((e2.second xor e3.second) and !(e5.second or e6.second)));
  else if (e1.second xor e4.second)
  {
    if (e5.second and e6.second and !(e2.second or e3.second)) // p4 v-w-v2-v1 is replace by p4 w-v-v1-v2
      return 2.f * boost::get(weight_t(),quotient,v) * boost::get(weight_t(), quotient, v1) * boost::get(weight_t(), quotient, v2) * boost::get(weight_t(), quotient, w);
    else if (!(e5.second or e6.second) and e2.second and e3.second)
      return 2.f * boost::get(weight_t(),quotient,v) * boost::get(weight_t(), quotient, v1) * boost::get(weight_t(), quotient, v2) * boost::get(weight_t(), quotient, w); // symmetric version
  }
  else if (!e1.second and !e4.second)
  {
    return (float)
    ((e2.second and e3.second and e5.second and !e6.second) 
    or (e2.second and e3.second and !e5.second and e6.second) 
    or (e2.second and !e3.second and e5.second and e6.second) 
    or (!e2.second and e3.second and e5.second and e6.second)) * boost::get(weight_t(),quotient,v) * boost::get(weight_t(), quotient, v1) * boost::get(weight_t(), quotient, v2) * boost::get(weight_t(), quotient, w);
  }
  return 0.f;
}

// Given a prime module mod, merge two childrens which are scored best
int mergeModules_adjusted(const Tree& tree, const TVertex& mod, Graph& g, unsigned int seed, bool debug = false)
{
  QGraph quotient = repGraph(tree,mod,g); // create representative graph/quotient graphs
  findP4s(quotient); // scores edges in the quotient graph
  std::pair<std::map<QVertex,QVertex>,QGraph> c = getComplement(quotient);
  QGraph complement = c.second;
  std::map<QVertex,QVertex> m = c.first;
  findP4s(complement); // complement gets scored as well
  
  float score = -(float)boost::num_vertices(g) * boost::num_vertices(g) * boost::num_vertices(g) * boost::num_vertices(g);
  float final_size = 0.f;
  std::vector<std::pair<QVertex,QVertex> > modules;
	std::pair<QVertex, QVertex> modulesf;
	std::vector<std::pair<QVertex,QVertex> > min_opf;
  std::vector<std::vector<std::pair<QVertex,QVertex> > > min_op;
  std::vector<float> score_list;
  
  for (const auto& v : boost::vertices(quotient))
    for (const auto& w : boost::vertices(quotient))
    {
      if (v >= w) continue; // do not merge a module with itself, nor consider pairs twice
//       std::cout << "Scoring " << boost::get(cont_t(), quotient,v)[0] << " " << boost::get(cont_t(), quotient,w)[0] << std::endl;
      
      std::vector<QVertex> diff1;
      std::vector<QVertex> diff2;
      getSymmetricDifference(v,w,quotient, diff1, diff2);
      
      // we dont want the other vertex to be part of this set
      if (diff1.size() == 0 and diff2.size() == 0)
      {
				std::cout << boost::get(cont_t(),quotient,v)[0] << " " << boost::get(cont_t(),quotient,w)[0] << std::endl;
				boost::write_graphviz(std::cout, quotient);
      }
      assert(diff1.size() > 0 or diff2.size() > 0); // the symmetric difference has to be > 0 else they should be the same module
      float tmp_score = 0.f;
      std::vector<std::pair<QVertex, QVertex> > tmp_op;
      std::vector<float> tmp_score_list;
      
      float sizes = 0.f;
      for (const auto& v1 : diff1)
      {
	if (v1 == w) continue; // we will be a series module
	auto pscore = scoreOperation(v1,v,w,quotient,complement,m,tmp_op);
        float p4s = pscore.first;
        for (const auto& p : tmp_op)
        {
          if (p == tmp_op.back()) break; // this is the one we just added
          if (p.first == tmp_op.back().first)
          {
            if (p.first == v) // both edges are neighboured to v
              p4s -= checkDoubleP4s(v,p.second, tmp_op.back().second, quotient);
            else if (p.first == w) // both edges are neighboured to w
              p4s -= checkDoubleP4s(w,p.second, tmp_op.back().second, quotient);
          }
          else
          { // one edge is neighboured to v the other to w
            p4s -= checkDoubleP4s_case2(v,w,p.second, tmp_op.back().second, quotient);
          }
        }
				tmp_score += p4s;
        tmp_score_list.push_back(p4s);
        sizes += pscore.second;
      }
      
      // same for the other way round ;)
      for (const auto& v1 : diff2)
      {
	if (v1 == v) continue; // we will be a series module
	auto pscore = scoreOperation(v1,w,v,quotient,complement,m,tmp_op);
        float p4s = pscore.first;
        for (const auto& p : tmp_op)
        {
          if (p == tmp_op.back()) break; // this is the one we just added
            if (p.first == tmp_op.back().first)
            {
              if (p.first == v)
                p4s -= checkDoubleP4s(v,p.second, tmp_op.back().second, quotient);
              else if (p.first == w)
                p4s -= checkDoubleP4s(w,p.second, tmp_op.back().second, quotient);
            }
            else
            {
              p4s -= checkDoubleP4s_case2(v,w,p.second, tmp_op.back().second, quotient);
            }
        }
        tmp_score += p4s; // real number of p4s removed (takes quotient graph into account)
        tmp_score_list.push_back(p4s);
        sizes += pscore.second;
      }
      
      tmp_score /= sizes;
      //std::cout << "Modules " << boost::get(cont_t(), quotient,v)[0] << " / " << boost::get(cont_t(), quotient,w)[0] << " with score: " << tmp_score << " and size: " << sizes << " (tmp score/size is: " << score << "/" << final_size << ")" <<  std::endl;
			if (tmp_score > score or (tmp_score == score and sizes <= final_size))
      {
//         std::cout << "Modules " << boost::get(cont_t(), quotient,v)[0] << " / " << boost::get(cont_t(), quotient,w)[0] << " with score: " << tmp_score << " and size: " << sizes <<  std::endl;
				if (tmp_score == score and sizes == final_size)
				{
					min_op.push_back(tmp_op);
					modules.push_back(std::make_pair(v,w));
				} 
				else
				{
					min_op = {tmp_op};
					modules = {std::make_pair(v,w)};
				}
				final_size = sizes;
				score = tmp_score;
        score_list = tmp_score_list; // this is currently not correctly being updated TODO
			}
    }
	if (modules.size() > 1)
	{
		std::mt19937 generator;
		generator.seed(seed);
		std::uniform_int_distribution<unsigned int> uni(0,modules.size() - 1);
		unsigned int pos = uni(generator);
		modulesf = modules[pos];
		min_opf = min_op[pos];
	}
	else
	{
		modulesf = modules[0];
		min_opf = min_op[0];
	}
//   std::cout << "Score: " << score*min_op.size() << std::endl;
  int edits = 0;
  // DEBUG
  if (false)
  {
    for (const auto& v : boost::get(cont_t(), quotient, modulesf.first))
      std::cout << v << " ";
    std::cout << "merged with ";
    for (const auto& v : boost::get(cont_t(), quotient, modulesf.second))
      std::cout << v << " ";
    std::cout << std::endl;
    std::cout << "(scored " << score << ", individual scores ";
    for (const auto& s : score_list)
      std::cout << s << " ";
    std::cout << " with edits: ";
    for (const auto& e : min_opf)
      std::cout << ((boost::edge(e.first,e.second,quotient).second) ? "Removed: " : "Added: ") << boost::get(cont_t(),quotient, e.first)[0] << " - " << boost::get(cont_t(), quotient, e.second)[0] << " ";
//      std::cout << " (score reduced by: " << reduced << ")" << std::endl;
    std::cout << ")" << std::endl;
  }
  assert(modulesf.first != modulesf.second);
  
  for (const auto& e : min_opf)
    edits += editEdge(e.first, e.second, quotient, g);
  
//   std::cout << std::endl;
  return edits;
}

bool isCograph(const Tree& decomp)
{
  const auto& modType = boost::get(&vertex_prop::type, decomp);
  
  tvertex_iter vi, vi_end;
  for (boost::tie(vi,vi_end) = boost::vertices(decomp); vi != vi_end; ++vi)
  {
    if (boost::get(modType, *vi) == "Prime" or boost::get(modType,*vi) == "Spider")
    {
      return false;
    }
  }
  return true;
}

int editSpider (const Tree& tree, const TVertex& mod, Graph& g)
{
  int min_degree = boost::num_vertices(g);
  int max_degree = 0;
  taiter tai, tai_end, tai_;
  for (const auto& v : boost::adjacent_vertices(mod,tree))
  {
    if (tree[v].name != "Internal")
    {
      int deg = boost::degree(getCorr(v,tree),g);
      if (deg < min_degree) // Find out legs [thin spider]
				min_degree = deg;
      if (deg > max_degree)
				max_degree = deg;
    }
  }
  std::vector<TVertex> legs;
  std::vector<TVertex> body;
  for (const auto& v : boost::adjacent_vertices(mod,tree))
  {
    if (tree[v].name != "Internal" and boost::degree(getCorr(v,tree),g) == min_degree) // we are a leg
      legs.push_back(v);
    else if (tree[v].name != "Internal" and boost::degree(getCorr(v,tree),g) == max_degree)
      body.push_back(v);
  }
  int ct = legs.size() - 1;
  for (const auto& v : legs)
  {
    if (ct == 0) // remove all but ONE leg-body edge
      break;
    for (const auto& v_ : boost::adjacent_vertices(mod,tree))
    {
      if (boost::edge(getCorr(v,tree),getCorr(v_,tree),g).second)
      {
// 	std::cout << "Removed: " << getCorr(v,tree) << " - " << getCorr(v_,tree) << std::endl;
				boost::remove_edge(getCorr(v,tree), getCorr(v_,tree),g);
				ct--;
				break;
      }
    }
  }
  return (legs.size() - 1);
}

Graph cograph_editing_adjusted(Graph g, unsigned int seed)
{  
  ModDecomp mod;
  clock_t w = clock();
  Tree decomp = mod.decompose(g);
//   printModDecomp(decomp);
  w = clock() - w;
//   std::cout << "MD took " << w/1000 << "ms" << std::endl;
  
  tvertex_iter vi, vi_end;
  const auto& modType = boost::get(&vertex_prop::type, decomp);
  
  int edits = 0;
  int iters = boost::num_vertices(g); // this should be sufficient (for graphs <1000 vertices)
  while (!isCograph(decomp) and iters > 0)
  {    
//     std::cout << "checked cograph" << std::endl;
//     printModDecomp(decomp);
    for (boost::tie(vi,vi_end) = boost::vertices(decomp); vi != vi_end; ++vi)
    {
      if (boost::get(modType, *vi) == "Spider")
      {
				// EDP4
				TVertex module = *vi;
				edits += editSpider(decomp, module, g);
      }
      if (boost::get(modType, *vi) == "Prime") // edit here!
      {
				TVertex module = *vi;
				w = clock();
				int n_edits = mergeModules_adjusted(decomp, module, g, seed);
// 	int n_edits = edit(decomp, module, g);
				edits += n_edits;
				if (n_edits == 0)
				{
	  			std::cerr << "Unable to edit graph: " << std::endl;
	  			printModDecomp(decomp);
	  			return g;
				}
			w = clock() - w;
// 	std::cout << "Merging took " << w/1000 << "ms" << std::endl;
      }
    }
//     boost::write_graphviz(std::cout, g);
    decomp = mod.decompose(g);
//     std::ofstream ss(std::string("dump/dump_") + boost::lexical_cast<std::string>(500-iters) + ".dot");
//     printModDecomp(decomp,ss);
//     std::ofstream s2(std::string("dump/dump_graph_")  + boost::lexical_cast<std::string>(500-iters) + ".dot");
//     boost::write_graphviz(s2,g);
    iters--;
  }
  if (iters == 0)
    std::cerr << "Iteration limit reached" << std::endl; // this means something went wrong (orly)
//   printModDecomp(decomp);
//   boost::write_graphviz(std::cout,g);
//   std::cout << edits << std::endl;
  //Tree after = mod.decompose(g);
  return g;
}

int main(int, char* argv[])
{
  
  /* build and edit */
//   BuildGraphs builder(atoi(argv[3]));
//   ModDecomp dec;
//   Graph g = builder.buildCograph(atoi(argv[1]));
//   std::ofstream ss(std::string("dump/Original.dot"));
//   printModDecomp(dec.decompose(g),ss);
//   g = builder.editEdges(g,atof(argv[2]));
//   ss.close();
//   ss.open(std::string("dump/Edited.dot"));
//   printModDecomp(dec.decompose(g),ss);
//   builder.writeGraphs(std::vector<Graph>{g}, "dump/graph.txt");
//   ss.close();
//   clock_t c = clock();
//   Graph h = cograph_editing_adjusted(g);
//   c = clock() - c;
//   std::cout << c/1000 << " millisec" << std::endl;
//   Tree t = dec.decompose(h);
//   ss.open(std::string("dump/Final.dot"));
//   printModDecomp(t,ss);
//   ss.close();
   
  /* run */
  // random is just a quick&dirty solution
  unsigned int seed = atoi(argv[3]);
  if (!seed)
  {
  	std::random_device rd;
  	seed = rd();
  }

  BuildGraphs builder(seed,false); // no random required
  std::vector<Graph> graphs = builder.readGraphs(argv[1]);
  clock_t t = clock();
  int i = 0;
  int tot = 0;
	std::vector<Graph> out;
  std::vector<std::pair<int,int> > props;
  for (auto& g : graphs)
  {
    //std::cerr << "Component: " << (++i) << std::endl;
    t = clock();
    Graph h = cograph_editing_adjusted(g, seed);
    t = (clock() - t)/1000;
    int edits = builder.calculateEditDistance(g,h);
    tot += edits;
		props.push_back(std::make_pair(edits,t));
    out.push_back(h);
  } 
	std::cout << tot << std::endl;
  builder.writeGraphs(out, argv[2], props);
  
  /* run MD */
//   BuildGraphs builder(0, false); // no random required
//   std::vector<Graph> graphs = builder.readGraphs(argv[1]);
//   ModDecomp dec;
//   for (auto& g : graphs)
//   {
//     Tree t = dec.decompose(g);
//     std::cout << builder.getNewick(t) << std::endl;
//   }


  /* run single */
//  ModDecomp mod;
//  BuildGraphs builder(0); // no random required
//  std::vector<Graph> graphs = builder.readGraphs(argv[1]);
//   std::vector<Graph> origs = builder.readGraphs(argv[2]);
//  clock_t t = clock();
//  Graph g = graphs[atoi(argv[2])];
//   Graph h = origs[atoi(argv[3])];
//  Tree orig = mod.decompose(g);
//  printModDecomp(orig);
//  Graph edited = cograph_editing_adjusted(g);
//  int edits = builder.calculateEditDistance(g,edited);
//  Tree two = mod.decompose(edited);
//  printModDecomp(two);
//  t = (clock() - t)/1000;
//  builder.writeGraphs(std::vector<Graph>{g},argv[3],edits,t);

  /* create testset */
//   BuildGraphs builder(atoi(argv[3]));
//   std::vector<Graph> graphs;
//   std::vector<Graph> cographs;
//   for (unsigned int i = 0; i < 100; i++)
//   {
//     Graph g = builder.buildCograph(atoi(argv[1]));
//     cographs.push_back(g);
//     Graph h = builder.editEdges(g, atof(argv[2]));
//     graphs.push_back(h);
//   }
//   float per = atof(argv[2]) * 100;
//   builder.writeGraphs(cographs, std::string("../Test/Cographs/Cographs_") + argv[1] + "_" + boost::lexical_cast<std::string>(per) + ".txt");
//   builder.writeGraphs(graphs, std::string("../Test/Graphs_") + argv[1] + "_" + boost::lexical_cast<std::string>(per) + ".txt");
  
  
  return 0;
}
