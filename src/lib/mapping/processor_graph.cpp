/*
 * graph.cpp
 *
 *  Created on: Feb 16, 2017
 *      Author: Roland */

#include "./processor_graph_access.h"

#include <stdlib.h>
#include <math.h>
#include <map>
#include <bitset>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <queue>

using namespace std;

/*****************************************/
/*****************************************/
void processor_graph_access::getDist2v(processor_graph_access & P, vector<int64_t> &dist2v, int64_t v) {
  forall_nodes(P, v) {
    dist2v[v] = UNDEFINED_NODE;
  } endfor
    
      queue<int64_t> q;
  q.push(v);
  dist2v[v] = 0;
  while(q.empty() == false) {
    int64_t f = q.front();
    int64_t newDist = 1 + dist2v[f];
    q.pop();
    forall_out_edges(P, e, f) {
      NodeID t = P.getEdgeTarget(e);
      if(dist2v[t] == UNDEFINED_NODE) {
	q.push(t);
	dist2v[t] = newDist;
      }
    } endfor
  }
}

/*****************************************/
/*****************************************/
bool processor_graph_access::setHammingLabels(processor_graph_access & P) {
  bool ret = true;
  int64_t edgeID_outer = 0;
  int64_t digit = 0;
  vector<bool> edgeMarked(P.number_of_edges(), false);
  vector<int64_t> dist2v(P.number_of_nodes());
  vector<int64_t> dist2w(P.number_of_nodes());
  P.resize_HammingLabels(P.number_of_nodes());
  forall_nodes(P, v) {
    forall_out_edges(P, e, v) {
      NodeID w = P.getEdgeTarget(e);
      if((v < w) && (edgeMarked[edgeID_outer] == false)) {
	//cout << "v and w are " << v << " and " << w << endl;
	P.getDist2v(P, dist2v, v);
	P.getDist2v(P, dist2w, w);
	EdgeID edgeID_inner = 0;
	forall_nodes(P, vv) {
	  forall_out_edges(P, ee, vv) {
	    NodeID ww = P.getEdgeTarget(ee);
	    if(vv < ww) {
	      //cout << "   vv and ww are " << vv << " and " << ww << endl;
	      if(dist2v[vv] == dist2v[ww]) {
		cout << "vertices " << vv << " and " << ww << " both have distance " << dist2v[vv] << " to " << v << endl;
		cout << "NOT A PARTIAL CUBE" << endl;
		return(false);
	      }
	      if(dist2w[vv] == dist2w[ww]) {
		cout << "vertices " << vv << " and " << ww << " both have distance " << dist2w[vv] << " to " << w << endl;
		cout << "NOT A PARTIAL CUBE" << endl;
		return(false);
	      }
	      if((dist2v[vv] + dist2w[ww]) != (dist2w[vv] + dist2v[ww])) {
		if(edgeMarked[edgeID_inner] == true) {
		  cout << "Edge " << edgeID_inner << " is in more than one cutset" << endl;
		  cout << "NOT A PARTIAL CUBE" << endl;
		  return(false);
		}
		//cout << "      Marking the edge between " << vv << " and " << ww << endl;
		edgeMarked[edgeID_inner] = true;
	      }
	    }
	    edgeID_inner++;
	  } endfor
	} endfor
	forall_nodes(P, u) {
	  if(dist2w[u] < dist2v[u]) {
	    P.setHammingLabel(u, P.getHammingLabel(u) + (1 << digit));
	  }
	} endfor
	digit++;
	//cout << "digit ist jetzt " << digit << endl;
      }
      edgeID_outer++;
    } endfor
  } endfor
  return(ret);
}

/*****************************************/
/*****************************************/
processor_graph_access::processor_graph_access(void) {
}

/*****************************************/
/*****************************************/
processor_graph_access::~processor_graph_access(void) {
}
