/*
 * permutation.cpp
 *
 *  Created on: Jan 02, 2017
 *      Author: Roland
*/

#include "hierarchy.h"
#include "../../alma_graph_access.h"
#include "../utils/arith.h"

using namespace std;

/*****************************************/
/*****************************************/
Hierarchy::Hierarchy(int64_t nl) {
    numLevels = nl;
    graphsVec.resize(nl);
    parentsVec.resize(nl);
}

/*****************************************/
/*****************************************/
Hierarchy::~Hierarchy(void) {
  // for(int64_t l = 0; l < numLevels; l++) {
  //   (parentsVec[l]).resize(0);
  // }
  // parentsVec.resize(0);

  for(int64_t l = 1; l < numLevels; l++) {
    delete graphsVec[l];
  }
  // graphsVec.resize(0);
}

/***************************************************************/
/* All graphs, except the one on the finest level, are deleted */
/***************************************************************/
// Hierarchy::~Hierarchy(void) {
//     for(int64_t l = 1; l < numLevels; l++) {
//       delete(graphsVec[l]);
//       (parentsVec[l]).resize(0);
//     }
// }

/*****************************************/
/*****************************************/
int64_t Hierarchy::getNumberOfLevels(void) {
    return(numLevels);
}

/*****************************************/
/*****************************************/
NodeID Hierarchy::getNumberVerticesOnLevel(int64_t l) {
    return((graphsVec[l])->number_of_nodes());
}

/*****************************************/
/*****************************************/
void Hierarchy::displayParentsOnLevel(int64_t l) {
  for(int64_t i = 0; i < (int64_t)((parentsVec[l]).size()); i++) {
        cout << (parentsVec[l])[i] << " ";
    }
    cout << endl;
}

/*****************************************/
/*****************************************/
void Hierarchy::insertGraph(alma_graph_access* newG_ptr, int64_t pos) {
    graphsVec[pos] = newG_ptr;
}

/*****************************************/
/*****************************************/
alma_graph_access Hierarchy::retFinestGraph(void) {
  return(*(graphsVec[0]));
}

/*****************************************/
/*****************************************/
void Hierarchy::insertParents(vector<NodeID>& parents, int64_t pos) {
    parentsVec[pos] = parents;
}

/*****************************************/
/*****************************************/
void Hierarchy::updateLabelsOnFinestLevel(void) {
  alma_graph_access & f = *(graphsVec[0]);
  int64_t incFromMostSigDigit = 1LL << numLevels;
  forall_nodes(f, v) {
    int64_t newLabel_v = (f.getLabel(v)) % 2;//least significant digit
    int64_t oldParentVertex = v;
    int64_t exp = 0;
    while(exp < numLevels - 1) {
      int64_t parentVertex = (parentsVec[exp])[oldParentVertex];
      exp++;
      int64_t parentLabel = (graphsVec[exp])->getLabel(parentVertex);
      int64_t prefCompleteLabel = newLabel_v + (parentLabel << exp);
      if(f.getVertexWithLabel(prefCompleteLabel) != UNDEFINED_NODE) {//vertex with such a label exists
	newLabel_v += (parentLabel % 2) << exp;	      
      } else {
	//cout << "second choice" << endl;
	int64_t otherCompleteLabel = flipPosition(exp, prefCompleteLabel);
	if(f.getVertexWithLabel(otherCompleteLabel) == UNDEFINED_NODE) {
	  cout << "stuck!" << endl;
	}
	newLabel_v += (1 - (parentLabel % 2)) << exp;
      }
      oldParentVertex = parentVertex;
    }
    if(f.getLabel(v) >= incFromMostSigDigit) {
      newLabel_v += incFromMostSigDigit;//most significant digit
    }
    f.setLabel(v, newLabel_v);
    if(v == 0) {
      //cout << "New Label of v is " << newLabel_v << endl;
    }
  } endfor
  f.makeLabel2Vertex(f);
}
