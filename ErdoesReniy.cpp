
//  compile with: g++ -Wall -ansi -pedantic main.cpp -o bla


#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/erdos_renyi_generator.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/random.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/edge_connectivity.hpp>


#include <boost/random.hpp>
#include <sys/time.h>


using namespace boost;
using namespace std;

//typedef boost::adjacency_list<> Graph;
typedef adjacency_list<listS,
                       vecS,
                       undirectedS,
                       no_property > Graph;
typedef boost::erdos_renyi_iterator<boost::minstd_rand, Graph> ERGen;
typedef graph_traits<Graph>::vertex_iterator vIt;

void Deg(const Graph& G){


	int maxD = 0;
	vector<int> DegSeq(num_vertices(G));
	vIt v = vertices(G).first; 
  for (; v != vertices(G).second; ++v)
		DegSeq[degree(*v,G)] = 0;	


	v = vertices(G).first; 
  for (; v != vertices(G).second; ++v){
		if(maxD < degree(*v,G)) maxD = degree(*v,G);
		DegSeq[degree(*v,G)]++;
	}


	cout << endl <<"Degree Seq:"<< endl;
	for(unsigned int i=0; i<=maxD; i++){
    if(i<10)	cout << " " << i << "|";
		else cout << i << "|";
    if(i%20==0 && i!=0) cout << endl;
	}

	cout << endl;
	for(unsigned int i=0; i<=maxD; i++){
		if(DegSeq[i]<10) cout << " " << DegSeq[i] << "|";
		else cout << DegSeq[i] << "|";
    if(i%20==0 && i!=0) cout << endl;
	}

	cout << endl;
	cout << endl << num_vertices(G) << " = |V|  " << num_edges(G) << " = |E|";
	cout << endl << "*********Thin Spider Test *********";
	cout << endl << DegSeq[1] <<  "= num_deg(1) = s"; 
	int v_s = num_vertices(G)-DegSeq[1];
	cout << endl << DegSeq[v_s] <<" = num_deg(v-s) = s?"; 
	
	cout << endl << "*********Thick Spider Test *********";
	cout << endl << DegSeq[num_vertices(G)-2] <<  "= num_deg(v-2) = s"; 
	int s = DegSeq[num_vertices(G)-2];
  if(s>0)	cout << endl << DegSeq[s-1] <<" = num_deg(s-1) = s?" ;
	else cout << endl << s-1;

}

void ConstructSpider(Graph& spider){

	int K = 10;	
	int S = K;
	int R=;

	for(int i=0; i<K; i++)
		for(int j=i+1; j<K; j++)
			add_edge(i, j, spider);

	int k = K;
	for(int i=0; i<K; i++){
		add_edge(i, k, spider);
		k++;
	}
	
	for(int i=0; i<K; i++)
		for(int j=k; j<k+R; j++)
	    		add_edge(i, j, spider);
   

}


int main()
{
  
	unsigned int seed(time(0));
  struct timeval tv;
  gettimeofday(&tv, NULL);
  srand(tv.tv_usec);
	//srand (time(NULL));

 //  mt19937 gen(time(0));
 // Graph g2; generate_random_graph(g2, 10, 10, gen, false);

  // Create graph with 100 nodes and edges with probability 0.05
  boost::minstd_rand gen(time(0));
  Graph g(ERGen(gen, 100, 0.05), ERGen(), 100);
  //boost::print_graph(g);

//	  mt19937 gen(time(NULL));
	  Graph g2;
	  generate_random_graph(g2, 100, 1500, gen, false);

	Graph spider;		
	ConstructSpider(spider);


	cout << endl << "++++++++++++++++++++++Erdoes-Renyi:++++++++++++++++++++++";
	Deg(g);


	cout << endl << endl << "++++++++++++++++++++++other++++++++++++++++++++++";
	Deg(g2);


	cout << endl << endl << "++++++++++++++++++++++SPIDER++++++++++++++++++++++";
	//cout << endl;	print_graph(spider); cout << endl;
	Deg(spider);


  return 0;
}
