#include <argtable2.h>
#include <cstdio>
#include <cmath>
#include <cstring> 
#include <fstream>
#include <iostream>
#include <iomanip>
#include <iostream>
#include <list>
#include <regex.h>
#include <sstream>
#include <string>

#include "./parse_parameters.h"
#include "data_structure/graph_access.h"
#include "quality_metrics.h"
#include "../lib/karma_config.h"
#include "configuration.h"
#include "graph_io.h"
#include "macros_assertions.h"
#include "random_functions.h"
#include "timer.h"

#define SQR_NR_BINS 10

/********************************************************************/
/*                                dfs                               */
/********************************************************************/
void dfs(NodeID n, graph_access & G, std::vector<bool> & marked){
  marked[n] =true;
  forall_out_edges(G, e, n){
    NodeID nn = G.getEdgeTarget(e);
    if(marked[nn] == false){
      dfs(nn, G, marked);
    }
  } endfor
}

/********************************************************************/
/*                 getNumberConnectedComponents                     */
/********************************************************************/
NodeID getNumberConnectedComponents(graph_access & G){
  std::vector<bool> marked(G.number_of_nodes(), false);
  NodeID ret = 0;
  forall_nodes(G,n){
    if(marked[n] == false){
      ret++;
      dfs(n, G, marked);
    }
  } endfor
  return(ret);
}

/********************************************************************/
/*                     getClusteringCoefficient                     */
/********************************************************************/
double getClusteringCoefficient(graph_access & G){
  long int nrTriplets = 0;
  long int nrClosedTriplets = 0;
  forall_nodes(G,n) {
    NodeID deg = G.getNodeDegree(n);
    if(deg > 1) {
      nrTriplets+= (deg * (deg - 1)) /2;
      std::vector<NodeID> NeighbsOfn(deg, UNDEFINED_NODE);
      NodeID i = 0;
      forall_out_edges(G, e, n) {
	NodeID nn = G.getEdgeTarget(e);
	NeighbsOfn[i] = nn;
	i++;
      } endfor
      for(NodeID i = 0; i < deg; i++) {
	for(NodeID j = i; j < deg; j++) {
	  //dealing with triplet n, NeighbsOfn[i], NeighbsOfn[j]
	  forall_out_edges(G, ee, NeighbsOfn[i]){
	    if(G.getEdgeTarget(ee) == NeighbsOfn[j]) {
	      nrClosedTriplets++;
	    }
	  } endfor
        }
      }
    }
  } endfor
  if(nrClosedTriplets > 0) {
    return(((double)nrClosedTriplets) / ((double)nrTriplets));
  } else {
    return(0.0);
  }
}

/********************************************************************/
/*              getAverageLocalClusteringCoefficient                */
/********************************************************************/
double getAverageLocalClusteringCoefficient(graph_access & G){
  std::vector<bool> isNeighbOfn(G.number_of_nodes(), false);
  double ret = 0.0;
  forall_nodes(G,n) {
    EdgeWeight deg = G.getNodeDegree(n);
    if(deg > 1) {
      forall_out_edges(G, e, n) {
	NodeID nn = G.getEdgeTarget(e);
	isNeighbOfn[nn] = true;
      } endfor
      unsigned long int nrLinks = 0;
      forall_out_edges(G, e, n) {
	NodeID nn = G.getEdgeTarget(e);
	forall_out_edges(G, ee, nn){
	  NodeID nnn = G.getEdgeTarget(ee);
	  if(isNeighbOfn[nnn] == true) {
	    nrLinks++;
	  }
	} endfor
      } endfor
      ret+= ((double) nrLinks) / ((double)(deg * (deg - 1)));

      //reset isNeighbOfn
      forall_out_edges(G, e, n) {
	NodeID nn = G.getEdgeTarget(e);
	isNeighbOfn[nn] = false;
      } endfor
    }
  } endfor 
  return(ret / ((double) G.number_of_nodes()));
}

/********************************************************************/
/*                        getMinMaxAvgDegrees                          */
/********************************************************************/
void getMinMaxAvgDegrees(graph_access & G, EdgeWeight & minDegree, EdgeWeight & maxDegree, EdgeWeight & avgDegree) {
  minDegree = G.number_of_nodes();
  maxDegree = 0;
  avgDegree = 0;
  forall_nodes(G,n) {
    //printf("Node %d has degree=%d\n",n,G.getNodeDegree(n));  
    avgDegree = avgDegree + G.getNodeDegree(n);
    if(G.getNodeDegree(n) < minDegree) {
      minDegree = G.getNodeDegree(n);
    }
    if(G.getNodeDegree(n) > maxDegree) {
      maxDegree = G.getNodeDegree(n);
    }
  } endfor
  avgDegree = avgDegree / G.number_of_nodes();
}

/********************************************************************/
/*                   writeUnweightedDegreeHisto                     */
/********************************************************************/
void writeUnweightedDegreeHisto(graph_access & G, bool cumul, EdgeWeight minDeg, EdgeWeight maxDeg){
  double spanDeg = (double) (maxDeg - minDeg);
  int nrBins = SQR_NR_BINS * SQR_NR_BINS;
  int nrBinsDec = nrBins - 1;
  double stretch = (double) nrBins;
  std::vector<long unsigned int> histo(nrBins, 0);
  //EdgeWeight notTooHigh = 0;
  forall_nodes(G,n){
    EdgeWeight deg = G.getNodeDegree(n);
    int bin = (int) (stretch * ((double)(deg - minDeg)) / spanDeg);
    if(bin < 0) bin = 0;
    if(bin > nrBinsDec) bin = nrBinsDec;
    (histo[bin])++;
  } endfor
  if(cumul == true){
    for(int i = 1; i < nrBins; i++){
      histo[i] += histo[i - 1];
    }
  }
  int bin = 0;
  for(int i = 0; i < SQR_NR_BINS; i++){
    for(int j = 0; j < SQR_NR_BINS; j++){
      printf("%f  ", ((double) histo[bin]) / ((double) G.number_of_nodes()));
      bin++;
    }
    printf("\n");
  }
//   printf("Fraction of nodes with degree <= MAX_DEG_FOR_KEEP_SEARCHING is %f\n",
// 	 ((double) notTooHigh) / ((double)G.number_of_nodes()));
}

/********************************************************************/
/*                          characterizeG                           */
/********************************************************************/
void characterizeG(graph_access & G, bool cumul){
  NodeID ncc = getNumberConnectedComponents(G);
  if(ncc == 1){
    printf("Graph is CONNECTED.\n");
  }
  else{
    printf("Graph is disconnected and has %u connected components.\n", ncc);
  }

  printf("Number of nodes: %u\n", G.number_of_nodes());
  printf("Number of edges: %u\n", G.number_of_edges() / 2);
  printf("Density is %f\n", ((double) (G.number_of_edges()) / 2) / ((double) G.number_of_nodes()));
  EdgeWeight minDeg, maxDeg, avgDeg;
  printf("Clustering coefficient is %f\n", getClusteringCoefficient(G));
  printf("Average local clustering coefficient is %f\n", getAverageLocalClusteringCoefficient(G));
  getMinMaxAvgDegrees(G, minDeg, maxDeg, avgDeg);
  printf("Minimum node degree: %d\n", minDeg);
  printf("Maximum node degree: %d\n", maxDeg);
  printf("Average node degree: %d\n", avgDeg);
  if(cumul == true){
    printf("Cummulative distribution of unweighted node degrees:\n");
  }
  else{
    printf("Frequency distribution of unweighted node degrees:\n");
  }
  writeUnweightedDegreeHisto(G, true, minDeg, maxDeg);
}

/********************************************************************/
/******************************* main *******************************/
/********************************************************************/
int main(int argn, char** argv) {
  PartitionConfig partition_config, repartition_config;
  std::string graph_filename;
  graph_access G;

  //parsing the parameters
  bool is_graph_weighted = false;
  bool suppress_output   = false;
  bool recursive         = false;
  
  int ret_code = parse_parameters(argn, argv, 
				  partition_config, 
				  graph_filename, 
				  is_graph_weighted, 
				  suppress_output, recursive);
  if(ret_code) {
    return(0);
  }
  
  //building G
  graph_io::readGraphWeighted(G, graph_filename);                                                         
  printf("%s\n",graph_filename.c_str());
  characterizeG(G, true);
  printf("\n");
  printf("\n");


  return 0;
}
