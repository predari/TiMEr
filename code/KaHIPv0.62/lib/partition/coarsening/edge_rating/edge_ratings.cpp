/******************************************************************************
 * edge_ratings.cpp 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 ******************************************************************************
 * Copyright (C) 2013 Christian Schulz <christian.schulz@kit.edu>
 * From "djoko" on from Roland Glantz 
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

#include <math.h>
#include "edge_ratings.h"
#include "partition_config.h"       
#include "random_functions.h"
#include "timer.h"
#include "data_structure/union_find.h"
  
double totalWeight = 0.0;
double halfTotalWeight = 0.0;
std::vector<NodeID> BFSparent, MSTparent;
std::vector<EdgeID> MSTparentEdge;
std::vector<double> cutWeight;
edge_ratings::edge_ratings(const PartitionConfig & _partition_config) : partition_config(_partition_config){

}

edge_ratings::~edge_ratings() {

}

void edge_ratings::rate(graph_access & G, unsigned level) {
  //rate the edges
  if(level == 0 && partition_config.first_level_random_matching) {
    return;
  } else if(partition_config.matching_type == MATCHING_RANDOM_GPA && level < partition_config.aggressive_random_levels) {
    return;
  } 
  if(level == 0 && partition_config.rate_first_level_inner_outer && 
     partition_config.edge_rating != EXPANSIONSTAR2ALGDIST ) {
    
    rate_inner_outer(G);
    
  } else if(partition_config.matching_type != MATCHING_RANDOM) {
    switch(partition_config.edge_rating) {
    case WEIGHT:
      rate_weight(G);
      break;
    case EXPANSIONSTAR:
      rate_expansion_star(G);
      break;
    case PSEUDOGEOM:
      rate_pseudogeom(G);
      break;
    case EXPANSIONSTAR2:
      rate_expansion_star_2(G);
      break;
    case EXPANSIONSTAR2ALGDIST:
      rate_expansion_star_2_algdist(G);
      break;
    case RESISTANCE:
      rate_resistance(G, level);
      break;
    case ZERO:
      rate_zero(G);
      break;
    }
  }
}

//simd implementation is possible
void edge_ratings::compute_algdist(graph_access & G, std::vector<float> & dist) {
        for( unsigned R = 0; R < 3; R++) {
                std::vector<float> prev(G.number_of_nodes(), 0);
                forall_nodes(G, node) {
                        prev[node] = random_functions::nextDouble(-0.5,0.5); 
                } endfor

                std::vector<float> next(G.number_of_nodes(), 0);
                float w = 0.5;

                for( unsigned k = 0; k < 7; k++) {
                        forall_nodes(G, node) {
                                next[node] = 0;

                                forall_out_edges(G, e, node) {
                                        NodeID target = G.getEdgeTarget(e);
                                        next[node] += prev[target] * G.getEdgeWeight(e);
                                } endfor

                                float wdegree = G.getWeightedNodeDegree(node);
                                if(wdegree > 0) {
                                        next[node] /= (float)wdegree;

                                }
                        } endfor

                        forall_nodes(G, node) {
                                prev[node] = (1-w)*prev[node] + w*next[node];
                        } endfor

                }

                forall_nodes(G, node) {
                        forall_out_edges(G, e, node) {
                                NodeID target = G.getEdgeTarget(e);
                                //dist[e] = max(dist[e],fabs(prev[node] - prev[target]));
                                dist[e] += fabs(prev[node] - prev[target]) / 7.0;
                        } endfor
                } endfor
        }

        forall_nodes(G, node) {
                forall_out_edges(G, e, node) {
                        dist[e] += 0.0001;
                } endfor
        } endfor

}


void edge_ratings::rate_weight(graph_access & G) {
  forall_nodes(G,n) {
    forall_out_edges(G, e, n) {
      EdgeRatingType rating = ((EdgeRatingType) G.getEdgeWeight(e));
      G.setEdgeRating(e, rating);
    } endfor
  } endfor
}

void edge_ratings::rate_expansion_star_2_algdist(graph_access & G) {

        std::vector<float> dist(G.number_of_edges(), 0);
        compute_algdist(G, dist);

        forall_nodes(G,n) {
                NodeWeight sourceWeight = G.getNodeWeight(n);
                forall_out_edges(G, e, n) {
                        NodeID targetNode = G.getEdgeTarget(e);
                        NodeWeight targetWeight = G.getNodeWeight(targetNode);
                        EdgeWeight edgeWeight = G.getEdgeWeight(e);

                        EdgeRatingType rating = 1.0*edgeWeight*edgeWeight / (targetWeight*sourceWeight*dist[e]);
                        G.setEdgeRating(e, rating);
                } endfor
        } endfor
}


void edge_ratings::rate_expansion_star_2(graph_access & G) {
        forall_nodes(G,n) {
                NodeWeight sourceWeight = G.getNodeWeight(n);
                forall_out_edges(G, e, n) {
                        NodeID targetNode = G.getEdgeTarget(e);
                        NodeWeight targetWeight = G.getNodeWeight(targetNode);
                        EdgeWeight edgeWeight = G.getEdgeWeight(e);

                        EdgeRatingType rating = 1.0*edgeWeight*edgeWeight / (targetWeight*sourceWeight);
                        G.setEdgeRating(e, rating);
                } endfor
        } endfor
}

void edge_ratings::rate_inner_outer(graph_access & G) {
        forall_nodes(G,n) {
#ifndef WALSHAWMH
                EdgeWeight sourceDegree = G.getWeightedNodeDegree(n);
#else
                EdgeWeight sourceDegree = G.getNodeDegree(n);
#endif
                if(sourceDegree == 0) continue;

                forall_out_edges(G, e, n) {
                        NodeID targetNode = G.getEdgeTarget(e);
#ifndef WALSHAWMH
                        EdgeWeight targetDegree = G.getWeightedNodeDegree(targetNode);
#else
                        EdgeWeight targetDegree = G.getNodeDegree(targetNode);
#endif
                        EdgeWeight edgeWeight = G.getEdgeWeight(e);
                        EdgeRatingType rating = 1.0*edgeWeight/(sourceDegree+targetDegree - edgeWeight);
                        G.setEdgeRating(e, rating);
                } endfor
        } endfor
}

void edge_ratings::rate_expansion_star(graph_access & G) {
        forall_nodes(G,n) {
                NodeWeight sourceWeight = G.getNodeWeight(n);
                forall_out_edges(G, e, n) {
                        NodeID targetNode       = G.getEdgeTarget(e);
                        NodeWeight targetWeight = G.getNodeWeight(targetNode);
                        EdgeWeight edgeWeight   = G.getEdgeWeight(e);

                        EdgeRatingType rating = 1.0 * edgeWeight / (targetWeight*sourceWeight);
                        G.setEdgeRating(e, rating);
                } endfor
        } endfor
}

void edge_ratings::rate_pseudogeom(graph_access & G) {
        forall_nodes(G,n) {
                NodeWeight sourceWeight = G.getNodeWeight(n);
                forall_out_edges(G, e, n) {
                        NodeID targetNode       = G.getEdgeTarget(e);
                        NodeWeight targetWeight = G.getNodeWeight(targetNode);
                        EdgeWeight edgeWeight   = G.getEdgeWeight(e);
                        double random_term      = random_functions::nextDouble(0.6,1.0);
                        EdgeRatingType rating   = random_term * edgeWeight * (1.0/(double)sqrt((double)targetWeight) + 1.0/(double)sqrt((double)sourceWeight));
                        G.setEdgeRating(e, rating);
                } endfor
        } endfor
}

/********************************* djoko *********************************/

double edge_ratings::rand01(void){
  return(((double) rand()) / ((double)(RAND_MAX + 1.0)));
}

EdgeID edge_ratings::getEdge(NodeID s, NodeID t, graph_access & G){
  EdgeID retE = UNDEFINED_EDGE;
  forall_out_edges(G, e, s){
    if(G.getEdgeTarget(e) == t){
      retE = e;
      return(retE);
    }
  } endfor
  return(retE);
}

void edge_ratings::getReverseEdges(graph_access & G, std::vector<EdgeID> & revEdge){
  forall_nodes(G, s){
    forall_out_edges(G, e, s){
      NodeID t = G.getEdgeTarget(e);
      if(s < t){
	EdgeID revE = getEdge(t, s, G);
	revEdge[e] = revE;
	revEdge[revE] = e;
      }
    } endfor
  } endfor
}

void edge_ratings::globalBFS(NodeID rootBFS, graph_access & G, std::vector<bool> & inBFSTree){
  std::vector<bool> marked(G.number_of_nodes(), false);
  std::queue<NodeID> q;
  marked[rootBFS] = true;
  q.push(rootBFS);
  while(!q.empty()){
    NodeID v = q.front();
    q.pop();
    forall_out_edges_in_random_order(G, e, v){
      NodeID w = G.getEdgeTarget(e);
      if(marked[w] == false){
	marked[w] = true;
	BFSparent[w] = v;
	inBFSTree[e] = true;
	q.push(w);
      }
    } endfor
  }
}

DirEdge* edge_ratings::buildSortedEdgeList(graph_access & G,
					   std::vector<unsigned int> & bottleneckIndex1,
					   unsigned int maxBI1){
  DirEdge* deliArr[maxBI1 + 1];
  for(unsigned int i = 0; i <= maxBI1; i++){
    deliArr[i] = NULL;
  }
  
  forall_nodes(G, source) {
    forall_out_edges(G, e, source){
      EdgeID target = G.getEdgeTarget(e);
      if(source < target){
	DirEdge* newDirEdge_ptr = (DirEdge*)malloc(sizeof(DirEdge));
	newDirEdge_ptr->source = source; 
	newDirEdge_ptr->target = target;
	newDirEdge_ptr->edID = e;
	newDirEdge_ptr->next = deliArr[bottleneckIndex1[e]];
	deliArr[bottleneckIndex1[e]] = newDirEdge_ptr;
      }    
    } endfor
  } endfor

  unsigned int i = 0;
  while((i <= maxBI1) && (deliArr[i] == NULL)){
    i++;
  }
  DirEdge* ret = NULL;
  if(i <= maxBI1){
    ret = deliArr[i];
    DirEdge* lastNotNULL = deliArr[i];
    while(lastNotNULL->next != NULL){
      lastNotNULL = lastNotNULL->next;
    }
    for(unsigned int j = i + 1; j <= maxBI1; j++){
      if(deliArr[j] != NULL){
  	lastNotNULL->next = deliArr[j];
  	lastNotNULL = deliArr[j];
  	while(lastNotNULL->next != NULL){
  	  lastNotNULL = lastNotNULL->next;
  	}      
      }
    }
  }
  return(ret);
}

void edge_ratings::parentMaker(NodeID w, graph_access & G, std::vector<bool> & inMST){
  forall_out_edges(G, e, w){
    if(inMST[e] == true){
      NodeID t = G.getEdgeTarget(e);
      if(MSTparent[t] == UNDEFINED_NODE){
	MSTparent[t] = w;
	MSTparentEdge[t] = e;
	parentMaker(t, G, inMST);
      }
    }
  } endfor
}

void edge_ratings::globalMST(NodeID rootMST, graph_access & G, 
			     std::vector<EdgeID> & revEdge,
			     std::vector<unsigned int> & bottleneckIndex1,
			     unsigned int maxBI1,
			     std::vector<bool> & inMST){
  DirEdge *deli = buildSortedEdgeList(G, bottleneckIndex1, maxBI1);
  // EdgeID count = 0;
  // DirEdge *runDeli = deli;
  // while(runDeli != NULL){
  //   count++;
  //   runDeli = runDeli->next;
  // }
  // printf("Number of sorted edges is %d\n", count);
 
  //construct MST
  union_find uf(G.number_of_nodes());
  DirEdge *runDeli = deli;
  while(runDeli != NULL){
    NodeID source = runDeli->source;
    NodeID target = runDeli->target;
    int sourceComp = uf.Find(source);
    int targetComp = uf.Find(target);
    if(sourceComp != targetComp){
      uf.Union(source, target);
      inMST[runDeli->edID] = true;
    }
    runDeli = runDeli->next;
  }

  //insert reverse edges into MST
  for(unsigned int e = 0; e < G.number_of_edges(); e++){
    if(inMST[e] == true){
      inMST[revEdge[e]] = true;
    }
  }    
  //make MST a rooted tree
  MSTparent[rootMST] = 0;
  MSTparentEdge[rootMST] = 0;
  parentMaker(rootMST, G, inMST);
  MSTparent[rootMST] = UNDEFINED_NODE;
  MSTparentEdge[rootMST] = UNDEFINED_EDGE;
  // forall_nodes(G, u){
  //   forall_out_edges(G, e, u){
  //     if(inMST[e] == true){
  // 	NodeID t = G.getEdgeTarget(e);
  // 	printf("u, t and MSTparent[t] are %d, %d and %d\n", u, t, MSTparent[t]);
  //    }
  //   } endfor
  // } endfor
}

void edge_ratings::resetInTree(graph_access & G, std::vector<bool> & inTree){
  for(EdgeID e = 0; e < G.number_of_edges(); e++){
    inTree[e] = false;
  }
}

void edge_ratings::getLabels(NodeID w, graph_access & G,
			     std::vector<bool> & inMST, NodeID* currLabel_ptr,
			     std::vector<NodeID> & label, std::vector<NodeID> & maxLabelDescendants){
  forall_out_edges(G, e, w){
    if(inMST[e] == true){
      NodeID z = G.getEdgeTarget(e);
      if(label[z] == UNDEFINED_NODE){
	(*currLabel_ptr)++;
	label[z] = *currLabel_ptr;
	getLabels(z, G, inMST, currLabel_ptr, label, maxLabelDescendants);
	maxLabelDescendants[z] = *currLabel_ptr;
      }
    }
  } endfor
}

NodeID edge_ratings::commonAncestor(NodeID v1, NodeID v2,
				    std::vector<NodeID> & label, std::vector<NodeID> & maxLabelDescendants){
  NodeID l = label[v2];
  NodeID m = maxLabelDescendants[v2];
  NodeID a = v1;
  while((maxLabelDescendants[a] < m) || (label[a] > l)){
    a = MSTparent[a];
  }
  return(a);
}

void edge_ratings::leaf(NodeID v, graph_access & G, std::vector<bool> & inMST,
			std::vector<NodeID> & label, std::vector<NodeID> & maxLabelDescendants,
			std::vector<double> & intraWeight, std::vector<double> & interWeight,
			std::vector<double> & subtreeWeight){
  EdgeID parentEdge = UNDEFINED_EDGE;
  forall_out_edges(G, f, v){
    NodeID t = G.getEdgeTarget(f);
    if(inMST[f] == true){
      parentEdge = f;
    }
    else{
      double edgeWeight = (double) G.getEdgeWeight(f);
      interWeight[v] += edgeWeight;
      intraWeight[commonAncestor(v, t, label, maxLabelDescendants)] += edgeWeight;
    }
  } endfor
  subtreeWeight[v] = (double) G.getWeightedNodeDegree(v);
  if(parentEdge != UNDEFINED_EDGE){
    cutWeight[v] = interWeight[v] + G.getEdgeWeight(parentEdge);
  }
}

void edge_ratings::update(NodeID u, graph_access & G, std::vector<bool> & inMST, 
			  std::vector<NodeID> & label, std::vector<NodeID> & maxLabelDescendants,
			  std::vector<double> & intraWeight, std::vector<double> & interWeight,
			  std::vector<double> & subtreeWeight){
  EdgeID parentEdge = UNDEFINED_EDGE;
  forall_out_edges(G, f, u){
    NodeID t = G.getEdgeTarget(f);
    if(inMST[f] == true){
      if(label[t] < label[u]){
	parentEdge = f;
      }
      else{
	interWeight[u] += interWeight[t];
	subtreeWeight[u] += subtreeWeight[t];
     }
    }
    else{
      if((label[t] < label[u]) || (label[t] > maxLabelDescendants[u])){
	double edgeWeight = (double) G.getEdgeWeight(f);
	intraWeight[commonAncestor(u, t, label, maxLabelDescendants)] += edgeWeight;
	interWeight[u] += edgeWeight;
      }
    }
  } endfor  
  
  interWeight[u] -= intraWeight[u];
  subtreeWeight[u] += (double) G.getWeightedNodeDegree(u);
  if(parentEdge != UNDEFINED_EDGE){
    cutWeight[u] = interWeight[u] + G.getEdgeWeight(parentEdge);
  }
}

void edge_ratings::mainDFT(NodeID v, graph_access & G, std::vector<bool> & inMST,
			   std::vector<NodeID> & label, std::vector<NodeID> & maxLabelDescendants,
			   std::vector<double> & intraWeight, std::vector<double> & interWeight,
			   std::vector<double> & subtreeWeight){
  unsigned int numberChildren = 0;
  forall_out_edges(G, e, v){
    if(inMST[e] == true){
      NodeID w = G.getEdgeTarget(e);
      if(label[w] > label[v]){
	numberChildren++;
	mainDFT(w, G, inMST, label, maxLabelDescendants, intraWeight, interWeight, subtreeWeight);
      }
    }
  } endfor 
  if(numberChildren == 0){
    leaf(v, G, inMST, label, maxLabelDescendants, intraWeight, interWeight, subtreeWeight);
  }
  else{
    update(v, G, inMST, label, maxLabelDescendants, intraWeight, interWeight, subtreeWeight);
  }
}

unsigned int edge_ratings::computeBottleneckIndices_one(graph_access & G, 
							std::vector<EdgeID> & revEdge,
							std::vector<unsigned int> & bottleneckIndex1){
  std::vector<bool> inBFSTree(G.number_of_edges(), false);
  BFSparent.resize(G.number_of_nodes(), UNDEFINED_NODE);
  for(unsigned int i = 0; i < MAX_NUMBER_TREES; i++){
    NodeID randomVertex = (unsigned int) (0.5 + (rand01() * ((double) G.number_of_nodes() - 1)));
    BFSparent[randomVertex] = UNDEFINED_NODE; 
    globalBFS(randomVertex, G, inBFSTree);
    for(unsigned int e = 0; e < G.number_of_edges(); e++){
      if(inBFSTree[e] == true){
	inBFSTree[revEdge[e]] = true;
      }
    }    
    std::vector<NodeID> label(G.number_of_nodes(), UNDEFINED_NODE);
    std::vector<NodeID> maxLabelDescendants(G.number_of_nodes(), UNDEFINED_NODE);
    label[randomVertex] = 0;
    maxLabelDescendants[randomVertex] = G.number_of_nodes() - 1;
    NodeID currLabel = 0;
    getLabels(randomVertex, G, inBFSTree, &currLabel, label, maxLabelDescendants);
    forall_nodes(G, u){
      forall_out_edges(G, e, u){
	if(inBFSTree[e] == true){
	  NodeID t = G.getEdgeTarget(e);
	  if(label[u] < label[t]){//u is parent of t
	    (bottleneckIndex1[e]) ++;
	  }
	}
      } endfor
    } endfor
  resetInTree(G, inBFSTree);
  }  
  forall_nodes(G, u){
    forall_out_edges(G, e, u){
      NodeID t = G.getEdgeTarget(e);
      if(u < t){
	if(bottleneckIndex1[e] > bottleneckIndex1[revEdge[e]]){
	  bottleneckIndex1[e] = bottleneckIndex1[revEdge[e]];
	}
	else{
	  bottleneckIndex1[revEdge[e]] = bottleneckIndex1[e];
	}
      }
    } endfor
  } endfor

  unsigned int maxBI1 = 0;
  for(unsigned int e = 0; e < G.number_of_edges(); e++){
    if(bottleneckIndex1[e] > maxBI1){
      maxBI1 = bottleneckIndex1[e];
    }
  }
  return(maxBI1);
}

EdgeID edge_ratings::findProxyEdge1(graph_access & G,
				    std::vector<NodeID> & label, std::vector<NodeID> & maxLabelDescendants,
				    std::vector<bool> & inMST,
				    std::vector<unsigned int> & bottleneckIndex1,
				    NodeID source,
				    NodeID target){
  NodeID caNode = commonAncestor(source, target, label, maxLabelDescendants);
  unsigned int maxBottleneckIndex1 = 0;
  EdgeID proxyEdge = UNDEFINED_EDGE;
  while(source != caNode){
    EdgeID parentEdge = MSTparentEdge[source];
    if(bottleneckIndex1[parentEdge] >= maxBottleneckIndex1){
      maxBottleneckIndex1 = bottleneckIndex1[parentEdge];
      proxyEdge = parentEdge;     
    }
    source = MSTparent[source];
  }
  while(target != caNode){
    EdgeID parentEdge = MSTparentEdge[target];
    if(bottleneckIndex1[parentEdge] >= maxBottleneckIndex1){
      maxBottleneckIndex1 = bottleneckIndex1[parentEdge];
      proxyEdge = parentEdge;     
    }
    target = MSTparent[target];
  }
  return(proxyEdge);
}

EdgeID edge_ratings::findProxyEdgeR(graph_access & G,
				    std::vector<NodeID> & label, std::vector<NodeID> & maxLabelDescendants,
				    std::vector<bool> & inMST, std::vector<EdgeRatingType> & bottleneckIndexR,
				    NodeID source,
				    NodeID target){
  NodeID caNode = commonAncestor(source, target, label, maxLabelDescendants);
  unsigned int minBottleneckIndexR = 1.0;
  EdgeID proxyEdge = UNDEFINED_EDGE;
  while(source != caNode){
    EdgeID parentEdge = MSTparentEdge[source];
    if(bottleneckIndexR[parentEdge] <= minBottleneckIndexR){
      minBottleneckIndexR = bottleneckIndexR[parentEdge];
      proxyEdge = parentEdge; 
      //printf("source:::new proxyEdge\n");
    }
    source = MSTparent[source];
  }
  while(target != caNode){
    EdgeID parentEdge = MSTparentEdge[target];
    //printf("target:::bottleneckIndexR[parentEdge] is %f\n", bottleneckIndexR[parentEdge]);
    if(bottleneckIndexR[parentEdge] <= minBottleneckIndexR){
      minBottleneckIndexR = bottleneckIndexR[parentEdge];
      proxyEdge = parentEdge;     
      //printf("target:::new proxyEdge\n");
    }
    target = MSTparent[target];
  }
  return(proxyEdge);
}

void edge_ratings::computeBottleneckIndices_resistance(graph_access & G,
						       std::vector<EdgeID> & revEdge,
						       std::vector<unsigned int> & bottleneckIndex1,
						       unsigned int maxBI1,
						       std::vector<EdgeRatingType> & bottleneckIndexR){
  //compute globalMST
  std::vector<bool> inMST(G.number_of_edges(), false);
  MSTparent.resize(G.number_of_nodes(), UNDEFINED_NODE);
  for(NodeID i = 0; i < G.number_of_nodes(); i++){
    MSTparent[i] = UNDEFINED_NODE;
  }
  MSTparentEdge.resize(G.number_of_nodes(), UNDEFINED_EDGE);
  for(NodeID i = 0; i < G.number_of_nodes(); i++){
    MSTparentEdge[i] = UNDEFINED_EDGE;
  }
  NodeID randomVertex = (unsigned int) (0.5 + (rand01() * ((double) G.number_of_nodes() - 1)));//root of MST
  globalMST(randomVertex, G, revEdge, bottleneckIndex1, maxBI1, inMST);

  //compute bottleneckIndexR
  cutWeight.resize(G.number_of_nodes(), 0.0);
  cutWeight[randomVertex] = 0.0;    
  std::vector<NodeID> label(G.number_of_nodes(), UNDEFINED_NODE);
  std::vector<NodeID> maxLabelDescendants(G.number_of_nodes(), UNDEFINED_NODE);
  label[randomVertex] = 0;
  maxLabelDescendants[randomVertex] = G.number_of_nodes() - 1;

  NodeID currLabel = 0;
  getLabels(randomVertex, G, inMST, &currLabel, label, maxLabelDescendants);


  std::vector<double> intraWeight(G.number_of_nodes(), 0.0);
  std::vector<double> interWeight(G.number_of_nodes(), 0.0);
  std::vector<double> subtreeWeight(G.number_of_nodes(), 0.0);
  mainDFT(randomVertex, G, inMST, label, maxLabelDescendants, intraWeight, interWeight, subtreeWeight);

  forall_nodes(G, sourceE){
    forall_out_edges(G, e, sourceE){
      NodeID targetE = G.getEdgeTarget(e);
      if(label[sourceE] < label[targetE]){//sourceE is parent of targetE
	if(inMST[e] == true){
	  double minSWT = subtreeWeight[targetE];
	  if(minSWT > halfTotalWeight){
	    minSWT = totalWeight - minSWT;
	  }
	  bottleneckIndexR[e] = (cutWeight[targetE]) / minSWT;
	  bottleneckIndexR[revEdge[e]] = bottleneckIndexR[e];
	}
      }
    } endfor
  } endfor

  forall_nodes(G, sourceE){
    forall_out_edges(G, e, sourceE){
      NodeID targetE = G.getEdgeTarget(e);
      if(label[sourceE] < label[targetE]){//sourceE is parent of targetE
	if(inMST[e] == false){
	    EdgeID proxyEdge = findProxyEdgeR(G, label, maxLabelDescendants, inMST, bottleneckIndexR,
					      sourceE, targetE);
	    //<<(directed) proxyEdge must be such that source is parent of target!!!
	    NodeID targetProxyE = G.getEdgeTarget(proxyEdge);	
	    double minSWT = subtreeWeight[targetProxyE];
	    if(minSWT > halfTotalWeight){
	      minSWT = totalWeight - minSWT;
	    }
	    bottleneckIndexR[e] = (cutWeight[targetProxyE]) / minSWT;
	    //bottleneckIndexR[e] = (cutWeight[targetProxyE] + G.getEdgeWeight(e) - G.getEdgeWeight(proxyEdge)) / minSWT;
	    bottleneckIndexR[revEdge[e]] = bottleneckIndexR[e];
	}
      }
    } endfor
  } endfor

  // //make bottleneckIndexR symmetric 
  // for(EdgeID e = 0; e < G.number_of_edges(); e++){
  //   if(bottleneckIndexR[e] == 0.0){
  //     bottleneckIndexR[e] = bottleneckIndexR[revEdge[e]];
  //   }
  // }
  
  // //multiply bottleneckIndexR with (1 + bottleneckIndex1)
  // for(EdgeID e = 0; e < G.number_of_edges(); e++){
  //   bottleneckIndexR[e] *= (1.0 + ((double)(bottleneckIndex1[e])));
  // }

  // // find maxBottleneckIndexR
  // EdgeRatingType maxBottleneckIndexR = 0.0;
  // for(EdgeID e = 0; e < G.number_of_edges(); e++){
  //   if(maxBottleneckIndexR < bottleneckIndexR[e]){
  //     maxBottleneckIndexR = bottleneckIndexR[e];
  //   }
  // }

  // //flip bottleneckIndexR 
  // maxBottleneckIndexR += 1.0;
  // for(EdgeID e = 0; e < G.number_of_edges(); e++){
  //   bottleneckIndexR[e] = maxBottleneckIndexR - bottleneckIndexR[e];
  // }
}

void edge_ratings::rate_resistance(graph_access & G, unsigned level){
  if(level == 0) {
    totalWeight = 0.0;
    forall_nodes(G, u)
      {
	totalWeight += (double) G.getWeightedNodeDegree(u);
      } endfor
    halfTotalWeight = totalWeight / 2.0;
    
    std::vector<EdgeID> revEdge(G.number_of_edges(), UNDEFINED_EDGE);
    getReverseEdges(G, revEdge);
    
    std::vector<unsigned int> bottleneckIndex1(G.number_of_edges(), 0);
    unsigned int maxBI1 = computeBottleneckIndices_one(G, revEdge, bottleneckIndex1);
    
    std::vector<EdgeRatingType> bottleneckIndexR(G.number_of_edges(), 0.0);
    computeBottleneckIndices_resistance(G, revEdge, bottleneckIndex1, maxBI1, bottleneckIndexR);
    //writeBottleneckIndexHistoR(G, bottleneckIndexR, 1);
    
    forall_nodes(G,n) {
      NodeWeight sourceWeight = G.getNodeWeight(n);
      forall_out_edges(G, e, n) {
	NodeID targetNode = G.getEdgeTarget(e);
	NodeWeight targetWeight = G.getNodeWeight(targetNode);
	EdgeWeight edgeWeight = G.getEdgeWeight(e);
	//EdgeRatingType rating = (bottleneckIndexR[e] * edgeWeight) / ((double)(targetWeight*sourceWeight));
	EdgeRatingType rating = bottleneckIndexR[e];
	G.setEdgeRating(e, rating);
      } endfor
    } endfor
  } else {
    rate_zero(G);
  }
}

void edge_ratings::rate_zero(graph_access & G) {
  forall_nodes(G,n) {
    forall_out_edges(G, e, n) {
      EdgeRatingType rating = 0.0;
      G.setEdgeRating(e, rating);
    } endfor
  } endfor
}

