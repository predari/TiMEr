#include "./build_processor_graphs.h"
#include "./alma/utils/arith.h"

/******************************************************/
/*                       random01                     */
/******************************************************/
double random01(void){
  return(((double) rand()) / ((double)(RAND_MAX + 1.0)));
}


/******************************************************/
/*                        getPower                    */
/******************************************************/
int getPower(int size) {
  int twopowers = 2;
  int power = 1;

  while (twopowers < size) {
    twopowers *= 2;
    power++;
  }

  if (twopowers > size)  error("Invalid size.", 11);

  return power;

}

/******************************************************/
/*                     buildFatTree                   */
/******************************************************/
void buildFatTree(processor_graph_access & proc, vector<int> & extensions, vector<int> & bandwidths) {
  vector<int> nodes;
  vector<int> nodeWeights;
  vector<int> edges;
  vector<int> edgeWeights;

  int sizeOfThisLevel = 1;
  int startOfLastLevel = 0;
  int startOfThisLevel = 0;
  int startOfNextLevel = 1;
  int numEdges = 0;
  int numNodes = 1;
  
  nodes.push_back(0);
  nodeWeights.push_back(0);
  

  for (unsigned int i = 0; i < extensions.size(); i++) {    
    //add edges for all nodes in level i
    for (int j = 0; j < sizeOfThisLevel; j++) {
      //edge to parent
      if (i > 0) {
	edges.push_back(startOfLastLevel + (j / extensions[i-1]));
	edgeWeights.push_back(bandwidths[i-1]);
      }

      for (int k = 0; k < extensions[i]; k++) {
	edges.push_back(startOfNextLevel + k + j * extensions[i]);
	edgeWeights.push_back(bandwidths[i]);
	numEdges++;
      }
    }
    //add new nodes for level i+1
    for (int j = 0; j < (sizeOfThisLevel * extensions[i]) ; j++) {
      if (i + 1 == extensions.size()) {
	nodes.push_back(numEdges);
	nodeWeights.push_back(1);
      } else {
	nodes.push_back(numEdges + j * extensions[i + 1]);
	nodeWeights.push_back(0);
      }
      numNodes++;
      numEdges++;
    
    }
    startOfLastLevel = startOfThisLevel;
    startOfThisLevel += sizeOfThisLevel;
    sizeOfThisLevel *= extensions[i];
    startOfNextLevel += sizeOfThisLevel;
  }
  vector<int> vec;
  proc.treeSort = vec;

  for (int j = 0; j < sizeOfThisLevel; j++) {
    edges.push_back(startOfLastLevel + (j / extensions[extensions.size() - 1]));
    edgeWeights.push_back(bandwidths[bandwidths.size() - 1]);
    proc.treeSort.push_back(j);
  }
  nodes.push_back(numEdges);
  proc.splitDeg = extensions;

  proc.build_from_metis_weighted(numNodes, &nodes[0], &edges[0], &nodeWeights[0], &edgeWeights[0]);
  int nodesNo = proc.number_of_nodes();

  for (int i = 0; i < nodesNo; i++) {
    vector<int> row(nodesNo);
    proc.distanceMatrix.push_back(row);
    proc.predecessorMatrix.push_back(row);
  }
  for (int i = 0; i < nodesNo; i++) {
    for (int j = 0; j < nodesNo; j++) {
      bool done = false;
      int distance = 0;
      int start = i;
      int target = j;
      int tempsave;
      while (start != target) {
	if (start > target) {
	  tempsave = start;	  
	  distance += proc.getEdgeWeight(proc.get_first_edge(start));
	  start = proc.getEdgeTarget(proc.get_first_edge(start));
	  if ((start == j) && !done) {
	    proc.predecessorMatrix[i][j] = tempsave;
	  }
	} else {
	  distance += proc.getEdgeWeight(proc.get_first_edge(target)); 
	  target = proc.getEdgeTarget(proc.get_first_edge(target));
	  if (!done) {
	    proc.predecessorMatrix[i][j] = proc.getEdgeTarget(proc.get_first_edge(j));
	    done = true;
	  }
	}
      }
      proc.distanceMatrix[i][j] = distance;
    }
  }    
  proc.refine = false;
}


/******************************************************/
/*                  build2DTreeRec                    */
/******************************************************/
void build2DTreeRec(int start, int total, int size, vector<int>* tree) {
   if (size == 1) {
    tree->push_back(start);
  } else {
    int half = size / 2;
    build2DTreeRec(start, total, half, tree);
    build2DTreeRec(start + half, total, half, tree);
    build2DTreeRec(start + (total * half), total, half, tree);
    build2DTreeRec(start + half + (total * half), total, half, tree);
  }
}

/******************************************************/
/*                  build3DTreeRec                    */
/******************************************************/
void build3DTreeRec(int start, int total, int size, vector<int>* tree) {
  if (size == 1) {
    tree->push_back(start);
  } else {
    int half = size / 2;
    build3DTreeRec(start, total, half, tree);
    build3DTreeRec(start + half, total, half, tree);
    build3DTreeRec(start + (total * half), total, half, tree);
    build3DTreeRec(start + half + (total * half), total, half, tree);
    build3DTreeRec(start + (total * total * half), total, half, tree);
    build3DTreeRec(start + half + (total * total * half), total, half, tree);
    build3DTreeRec(start + (total * half) + (total * total * half), total, half, tree);
    build3DTreeRec(start + half + (total * half) + (total * total * half), total, half, tree);
  }
}

/******************************************************/
/*                   buildLevDist                     */
/******************************************************/
vector< vector<int> > buildLevDist(int level) {
  vector< vector<int> > toRet;

  for (int i = 0; i < level; i++) {
    vector<int> row(level);
    toRet.push_back(row);
  }

  if (level == 4) {
    for (int i = 0; i < level; i++) {
      for (int j = 0; j < level; j++) {
	int x = 0;
	int y = 0;
	if (i % 2 != j % 2) {
	  x = 1;
	}
	if (i / 2 != j / 2) {
	  y = 1;
	}
        toRet[i][j] = x + y;
      }
    }
  }

  if (level == 8) {
    for (int i = 0; i < level; i++) {
      for (int j = 0; j < level; j++) {
	int x = 0;
	int y = 0;
	int z = 0;
	if (i % 2 != j % 2) {
	  x = 1;
	}
	if (((i / 2) % 2) != ((j / 2) % 2)) {
	  y = 1;
	}
	if (i / 4 != j / 4) {
	  z = 1;
	}
	toRet[i][j] = x + y + z;
      }
    }
  }

  return toRet;
}

/******************************************************/
/*              elementOfWeightVector                 */
/******************************************************/
int elementOfWeightVector(int size, int dest, int target) {

  int counter = 0;
  int diff = size / 2;
  while (true) {
    if (dest / diff != target / diff) {
      return counter;
    } else {
      diff /= 2;
      counter++;
    }
  }
}

/******************************************************/
/*                    build2DGrid                     */
/******************************************************/
void build2DGrid (processor_graph_access & proc, int size, vector<int> & bandwidths) {
  int power = getPower(size);

  //A Grid needs at least 2 nodes per dimension.
  if(size < 2) error("Size too small! Needs to be at least 2.",40);
  if(bandwidths.size() != 1 && ((int) (bandwidths.size())) != power) error("Wrong size of weight vector.", 12);

  int* nodes;
  nodes = new int[size * size + 1];
  int* edges;
  edges = new int[(4 * size * size) - 4 * size];
  int* edgeWeights;
  edgeWeights = new int[(4 * size * size) - 4 * size];
  int* nodeWeights;
  nodeWeights = new int [size * size];

  int counter = 0;
  nodes[0] = 0;

  for (int i = 0; i < size * size; i++) {
    nodeWeights[i] = 1;

    int startX = i % size;
    int startY = i / size;    

    //left
    if (i % size == 0) {
      //do nothing
    } else {
      edges[counter] = i - 1;
      if (bandwidths.size() == 1) {
	edgeWeights[counter] = bandwidths[0];
      } else {
	edgeWeights[counter] = bandwidths[elementOfWeightVector(size, startX, startX - 1)];
      }
      counter++;
    }
    //right
    if (i % size == size - 1) {
      //do nothing
    } else {
      edges[counter] = i + 1;
      if (bandwidths.size() == 1) {
	edgeWeights[counter] = bandwidths[0];
      } else {
	edgeWeights[counter] = bandwidths[elementOfWeightVector(size, startX, startX + 1)];
      }
      counter++;
    }
    //up
    if (i / size == 0) {
      //do nothing
    } else {
      edges[counter] = i - size;
      if (bandwidths.size() == 1) {
	edgeWeights[counter] = bandwidths[0];
      } else {
	edgeWeights[counter] = bandwidths[elementOfWeightVector(size, startY, startY - 1)];
      }
      counter++;
    }
    //down
    if (i / size == size - 1) {
      //do nothing
    } else {
      edges[counter] = i + size;
      if (bandwidths.size() == 1) {
	edgeWeights[counter] = bandwidths[0];
      } else {
	edgeWeights[counter] = bandwidths[elementOfWeightVector(size, startY, startY - 1)];
      }
      counter++;
    }

    nodes[i+1] = counter;

  }

  proc.build_from_metis_weighted(size * size, &nodes[0], &edges[0], &nodeWeights[0], &edgeWeights[0]);
  
  int nodesNo = proc.number_of_nodes();
  printf("proc has %d nodes\n", proc.number_of_nodes());

  for (int i = 0; i < nodesNo; i++) {
    vector<int> row(nodesNo);
    proc.distanceMatrix.push_back(row);
    proc.predecessorMatrix.push_back(row);
  }

  for (int i = 0; i < nodesNo; i++) {
    for (int j = 0; j < nodesNo; j++) {

      int startX = i % size;
      int startY = i / size;
      int targetX = j % size;
      int targetY = j / size;
      
      int distX, distY;
      
      if (bandwidths.size() == 1) {
	distX = std::max(startX, targetX) - std::min(startX, targetX);
	distY = std::max(startY, targetY) - std::min(startY, targetY);
	distX *= bandwidths[0];
	distY *= bandwidths[0];
      } else {

	int maxX = max(startX, targetX);
	int minX = min(startX, targetX);
	int maxY = max(startY, targetY);
	int minY = min(startY, targetY);
	distX = 0;
	distY = 0;
	for (int k = minX; k < maxX; k++) {
	  distX += bandwidths[elementOfWeightVector(size, k, k+1)];
	}
	for (int k = minY; k < maxY; k++) {
	  distY += bandwidths[elementOfWeightVector(size, k, k+1)];
	}
      }

      proc.distanceMatrix[i][j] = distX + distY;
      
      if (i == j) {
	proc.predecessorMatrix[i][j] = 0;
      } else {
	double partOfY = ((double) distY) / ((double) proc.distanceMatrix[i][j]);
	double randy = random01();
	if((randy <= partOfY) && (distY > 0)) {
	    if (startY > targetY) {
	      proc.predecessorMatrix[i][j] = j + size;
	    } else {
	      proc.predecessorMatrix[i][j] = j - size;
	    }
	  }
	else {// in particular distX > 0
	  if (startX > targetX) {
	    proc.predecessorMatrix[i][j] = j + 1;
	  } else {
	    proc.predecessorMatrix[i][j] = j - 1;
	  }
	}  
      }
    }
  }

  vector<int> split(power, 4);
  proc.splitDeg = split;

  vector<int> vec;
  proc.treeSort = vec;
  build2DTreeRec(0, size, size, &(proc.treeSort));
  
  proc.refine = true;
  proc.levelDist = buildLevDist(4);

}

/******************************************************/
/*                    build3DGrid                     */
/******************************************************/
void build3DGrid (processor_graph_access & proc, int size, vector<int> & bandwidths) {
  int power = getPower(size);

  //A Grid needs at least 2 nodes per dimension.
  if(size < 2) error("Size too small! Needs to be at least 2.",40);
  if(bandwidths.size() != 1 && ((int) (bandwidths.size())) != power) error("Wrong size of weight vector.", 12);

	
  int* nodes;
  nodes = new int[size * size * size + 1];
  int* edges;
  edges = new int[(6 * size * size * size) - (6 * size * size)];
  int* edgeWeights;
  edgeWeights = new int[(6 * size * size * size) - (6 * size * size)];
  int* nodeWeights;
  nodeWeights = new int[size * size * size];

  int counter = 0;
  nodes[0] = 0;

  for (int i = 0; i < size * size * size; i++) {
    nodeWeights[i] = 1;
    
    int startX = i % size;
    int startY = (i % (size * size)) / size;
    int startZ = i / (size * size);

    //left
    if (i % size == 0) {
      //do nothing                            
    } else {
      edges[counter] = i - 1;
      if (bandwidths.size() == 1) {
	edgeWeights[counter] = bandwidths[0];
      } else {
	edgeWeights[counter] = bandwidths[elementOfWeightVector(size, startX, startX - 1)];
      }      
      counter++;
    }
    //right
    if (i % size == size - 1) {
      //do nothing
    } else {
      edges[counter] = i + 1;
      if (bandwidths.size() == 1) {
	edgeWeights[counter] = bandwidths[0];
      } else {
	edgeWeights[counter] = bandwidths[elementOfWeightVector(size, startX, startX + 1)];
      }
      counter++;
    }
    //up
    if ((i % (size * size)) / size == 0) {
      //do nothing
    } else {
      edges[counter] = i - size;
      if (bandwidths.size() == 1) {
	edgeWeights[counter] = bandwidths[0];
      } else {
	edgeWeights[counter] = bandwidths[elementOfWeightVector(size, startY, startY - 1)];
      }
      counter++;
    }
    //down
    if ((i % (size * size)) / size == size - 1) {
      //do nothing
    } else {
      edges[counter] = i + size;
      if (bandwidths.size() == 1) {
	edgeWeights[counter] = bandwidths[0];
      } else {
	edgeWeights[counter] = bandwidths[elementOfWeightVector(size, startY, startY + 1)];
      }
      counter++;
    }

    //front
    if (i / (size * size) == 0) {
      //do nothing
    } else {
      edges[counter] = i - (size * size);
      if (bandwidths.size() == 1) {
	edgeWeights[counter] = bandwidths[0];
      } else {
	edgeWeights[counter] = bandwidths[elementOfWeightVector(size, startZ, startZ - 1)];
      }
      counter++;
    }
    //back
    if (i / (size * size) == size - 1) {
      //do nothing	
    } else {
      edges[counter] = i + (size * size);
      if (bandwidths.size() == 1) {
	edgeWeights[counter] = bandwidths[0];
      } else {
	edgeWeights[counter] = bandwidths[elementOfWeightVector(size, startZ, startZ - 1)];
      }
      counter++;
    }
    nodes[i+1] = counter;	
  }

  proc.build_from_metis_weighted(size * size * size, &nodes[0], &edges[0], &nodeWeights[0], &edgeWeights[0]);
  
  int nodesNo = proc.number_of_nodes();

  for (int i = 0; i < nodesNo; i++) {
    vector<int> row(nodesNo);
    proc.distanceMatrix.push_back(row);
    proc.predecessorMatrix.push_back(row);
  }

  for (int i = 0; i < nodesNo; i++) {
    for (int j = 0; j < nodesNo; j++) {

      int startX = i % size;
      int startY = (i % (size * size)) / size;
      int startZ = i / (size * size);
      int targetX = j % size;
      int targetY = (j % (size * size)) / size;
      int targetZ = j / (size * size);

      int distX, distY, distZ;

      if (bandwidths.size() == 1) {
	distX = bandwidths[0] * (std::max(startX, targetX) - std::min(startX, targetX));
	distY = bandwidths[0] * (std::max(startY, targetY) - std::min(startY, targetY));
	distZ = bandwidths[0] * (std::max(startZ, targetZ) - std::min(startZ, targetZ));
      } else {
	int maxX = max(startX, targetX);
	int minX = min(startX, targetX);
	int maxY = max(startY, targetY);
	int minY = min(startY, targetY);
	int maxZ = max(startZ, targetZ);
	int minZ = min(startZ, targetZ);

	distX = 0;
	distY = 0;
	distZ = 0;

	for (int k = minX; k < maxX; k++) {
	  distX += bandwidths[elementOfWeightVector(size, k, k+1)];
	}
	for (int k = minY; k < maxY; k++) {
	  distY += bandwidths[elementOfWeightVector(size, k, k+1)];
	}
	for (int k = minZ; k < maxZ; k++) {
	  distZ += bandwidths[elementOfWeightVector(size, k, k+1)];
	}
      }
      
      proc.distanceMatrix[i][j] = distX + distY + distZ;
      
      if (i == j) {
	proc.predecessorMatrix[i][j] = 0;
      } else {
	double partOfZ = ((double) distZ) / ((double)(proc.distanceMatrix[i][j]));
	double randy = random01();
	if((randy <= partOfZ) && (distZ > 0)) {
	  if (startZ > targetZ) {
	    proc.predecessorMatrix[i][j] = j + size * size;
	  } else {
	    proc.predecessorMatrix[i][j] = j - size * size;
	  }
	} else {//in particular distX + distY > 0
	  double partOfY = ((double) distY) / ((double) distX + distY);
	  randy = random01();
	  if((randy <= partOfY) && (distY > 0)) {
	    if (startY > targetY) {
	      proc.predecessorMatrix[i][j] = j + size;
	    } else {
	      proc.predecessorMatrix[i][j] = j - size;
	    }
	  }
	  else {//in particular distX > 0
	    if (startX > targetX) {
	      proc.predecessorMatrix[i][j] = j + 1;
	    } else {
	      proc.predecessorMatrix[i][j] = j - 1;
	    }
	  }
	}
      }
    }
  }
    

  vector<int> split(power, 8);
  proc.splitDeg = split;

  vector<int> vec;
  proc.treeSort = vec;
  build3DTreeRec(0, size, size, &(proc.treeSort));

  proc.refine = true;
  proc.levelDist = buildLevDist(8);

}

/************************************************/
/************ allPathsShortestPaths *************/
/************************************************/
void allPathsShortestPaths (vector< vector<int> >* predecessorMatrix, vector< vector<int> >* distanceMatrix, graph_access* graph) {


  int nodes = graph->number_of_nodes();

  //Build predecessorMatrix and distanceMatrix.
  for (int i = 0; i < nodes; i++) {
    vector<int> row(nodes, 100000);
    distanceMatrix->push_back(row);
    predecessorMatrix->push_back(row);
  }

  //Set trivial distanceMatrixs
  forall_nodes((*graph), node) {
    (*distanceMatrix)[node][node] = 0;
    forall_out_edges((*graph), edge, node) {
      (*distanceMatrix)[node][graph->getEdgeTarget(edge)] = graph->getEdgeWeight(edge);
      (*predecessorMatrix)[node][graph->getEdgeTarget(edge)] = node;
      } endfor
  } endfor

  
  //Floyd-Warshall-Algorithm
  for (int k = 0; k < nodes; k++) {
    for (int i = 0; i < nodes; i++) {
      for (int j = 0; j < nodes; j++) {
	if ((*distanceMatrix)[i][k] + (*distanceMatrix)[k][j] < (*distanceMatrix)[i][j]) {
	  (*distanceMatrix)[i][j] = (*distanceMatrix)[i][k] + (*distanceMatrix)[k][j];
	  (*predecessorMatrix)[i][j] = (*predecessorMatrix)[k][j];
	}
      }
    }
  }
  
  //Print matrix for debug purposes.
  /*
   for (int i = 0; i < nodes; i++) {
    for (int j = 0; j < nodes; j++) {
      cout << " " << (*distanceMatrix)[i][j];
    }
    cout << "\n";
    }*/
  
}

/******************************************************/
/*                   build2DTorus                     */
/******************************************************/
void build2DTorus (processor_graph_access & proc, int size, vector<int> & bandwidths) {
  int power = getPower(size);

  //A Torus needs at least 3 nodes per dimension.
  if(size < 3) error("Size too small! Needs to be at least 3.", 20);

  if(bandwidths.size() != 1 && ((int) (bandwidths.size())) != power) error("Wrong size of weight vector.", 12);

  int* nodes;
  nodes = new int[size * size + 1];
  int* edges;
  edges = new int[4 * size * size];
  int* edgeWeights;
  edgeWeights = new int[4 * size * size];
  int* nodeWeights;
  nodeWeights = new int[size * size];

  for (int i = 0; i < size * size; i++) {
    int startX = i % size;
    int startY = i / size;

    // 4 edges per node
    nodes[i] = 4 * i;
    nodeWeights[i] = 1;
    edges[4 * i] = 0;
    edges[4 * i + 1] = 0;
    edges[4 * i + 2] = 0;
    edges[4 * i + 3] = 0;
    if (bandwidths.size() == 1) {
      edgeWeights[4 * i] = bandwidths[0];
      edgeWeights[4 * i + 1] = bandwidths[0];		
      edgeWeights[4 * i + 2] = bandwidths[0];		
      edgeWeights[4 * i + 3] = bandwidths[0];
    } else {
      edgeWeights[4 * i] = bandwidths[elementOfWeightVector(size, startX, startX - 1)];
      edgeWeights[4 * i + 1] = bandwidths[elementOfWeightVector(size, startX, startX + 1)];		
      edgeWeights[4 * i + 2] = bandwidths[elementOfWeightVector(size, startY, startY - 1)];		
      edgeWeights[4 * i + 3] = bandwidths[elementOfWeightVector(size, startY, startY + 1)];
    }
  }
  nodes[size * size] = size * size * 4;

  printf("still alive 1\n");
  if(size % 2 != 0){
    printf("Extensions of torus must be even\n");
    exit(444);
  }

  //connects node i with his neighbours
  for (int i = 0; i < size * size; i++) {

    //left
    if (i % size == 0) {
      edges[4 * i] = i - 1 + size;
    } else {
      edges[4 * i] = i - 1;
    }
    //right
    if (i % size == size - 1) {
      edges[4 * i + 1] = i + 1 - size;
    } else {
      edges[4 * i + 1] = i + 1;
    }
    //up
    if (i / size == 0) {
      edges[4 * i + 2] = i - size + (size * size);
    } else {
      edges[4 * i + 2] = i - size;
    }
    //down
    if (i / size == size - 1) {
      edges[4 * i + 3] = i + size - (size * size);
    } else {
      edges[4 * i + 3] = i + size;
    }
  }
  proc.build_from_metis_weighted(size * size, &nodes[0], &edges[0], &nodeWeights[0], &edgeWeights[0]);
  proc.dim = 2;

  int nodesNo = proc.number_of_nodes();
  printf("proc has %d nodes\n", proc.number_of_nodes());

  for (int i = 0; i < nodesNo; i++) {
    vector<int> row(nodesNo);
    proc.distanceMatrix.push_back(row);
    proc.predecessorMatrix.push_back(row);
  }

  // vector< vector<int> > pred;
  // vector< vector<int> > dist;
  // int nodesNo = proc.number_of_nodes();

  // for (int i = 0; i < nodesNo; i++) {
  //   vector<int> row(nodesNo, 0);
  //   dist.push_back(row);
  //   pred.push_back(row);
  // }

  for (int i = 0; i < nodesNo; i++) {
    for (int j = 0; j < nodesNo; j++) {
      int startX = i % size;
      int startY = i / size;
      int targetX = j % size;
      int targetY = j / size;
      bool jumpX = false;
      bool jumpY = false;

      int distX; 
      int distY;

      if (bandwidths.size() == 1) {
	distY = std::max(startY, targetY) - std::min(startY, targetY);
	distX = std::max(startX, targetX) - std::min(startX, targetX);

     	if (distX * 2 >= size) {
	  distX = size - distX;
	  jumpX = true;
	}
	if (distY * 2 >= size) {
	  distY = size - distY;
	  jumpY = true;
	}
	distX *= bandwidths[0];
	distY *= bandwidths[0];
      } else {
	int maxX = max(startX, targetX);
	int minX = min(startX, targetX);
	int maxY = max(startY, targetY);
	int minY = min(startY, targetY);
	
	if ((maxY - minY) * 2 >= size) {
	  jumpY = true;
	} 
	if ((maxX - minX) * 2 >= size) {
	  jumpX = true;
	}

	int xPath = maxX;
	distX = 0;
	while (xPath != minX) {
	  if (!jumpX) {
	    distX += bandwidths[elementOfWeightVector(size, xPath, xPath - 1)];
	    xPath--;
	  } else {
	    distX += bandwidths[elementOfWeightVector(size, xPath, xPath + 1)];
	    if (xPath + 1 == size) {
	      xPath = 0;
	    } else {
	      xPath++;
	    }
	  }
	}

	int yPath = maxY;
	distY = 0;
	while (yPath != minY) {
	  if (!jumpY) {
	    distY += bandwidths[elementOfWeightVector(size, yPath, yPath - 1)];
	    yPath--;
	  } else {
	    distY += bandwidths[elementOfWeightVector(size, yPath, yPath + 1)];
	    if (yPath + 1 == size) {
	      yPath = 0;
	    } else {
	      yPath++;
	    }
	  }
	}
      }
      
      proc.distanceMatrix[i][j] = (distX + distY);
      
      if (i == j) {
	proc.predecessorMatrix[i][j] = 0;
      } else {
	double partOfY = ((double) distY) / ((double) proc.distanceMatrix[i][j]);
	double randy = random01();
	if((randy <= partOfY) && (distY > 0)) {
	  if (((startY > targetY) && !jumpY) || ((startY < targetY) && jumpY)) {
	    if (j + size < size * size) {
	      proc.predecessorMatrix[i][j] = j + size;
	    } else {
	      proc.predecessorMatrix[i][j] = j + size - (size * size);
	    }
	  } else {
	    if (j - size >= 0) {
	      proc.predecessorMatrix[i][j] = j - size;
	    } else {
	      proc.predecessorMatrix[i][j] = j - size + (size * size);
	    }
	  }
	} else {//in particular distX > 0
	  if (((startX > targetX) && !jumpX) || ((startX < targetX) && jumpX)) {
	    if ((j % size) != (size - 1)) {
	      proc.predecessorMatrix[i][j] = j + 1;
	    } else {
	      proc.predecessorMatrix[i][j] = j + 1 - size;
	    }
	  } else {
	    if ((j % size) > 0) {
	      
	      proc.predecessorMatrix[i][j] = j - 1;
	    } else {
	      proc.predecessorMatrix[i][j] = j - 1 + size;
	    }
	  }
	}
      }
     }
  }
  
  vector<int> split(power, 4);
  proc.splitDeg = split;
  
  vector<int> vec;
  proc.treeSort = vec;
  build2DTreeRec(0, size, size, &(proc.treeSort));
  proc.refine = true;
  proc.levelDist = buildLevDist(4);
}

/******************************************************/
/*                   build3DTorus                     */
/******************************************************/
void build3DTorus (processor_graph_access & proc, int size, vector<int> & bandwidths) {
  int power = getPower(size);

  //A Torus needs at least 3 nodes per dimension.
  if(size < 3)	error("Size too small! Needs to be at least 3.",30);
  if(bandwidths.size() != 1 && ((int) (bandwidths.size())) != power) error("Wrong size of weight vector.", 12);
	
  int* nodes;
  nodes = new int[size * size * size + 1];
  int* edges;
  edges = new int[6 * size * size * size];
  int* edgeWeights;
  edgeWeights = new int[6 * size * size * size];
  int* nodeWeights;
  nodeWeights = new int [size * size * size];

  for (int i = 0; i < size * size * size; i++) {

    // 6 edges per node
    nodes[i] = 6 * i;
    nodeWeights[i] = 1;

    int startX = i % size;
    int startY = (i % (size * size)) / size;
    int startZ = i / (size * size);

    edges[6 * i] = 0;
    edges[6 * i + 1] = 0;
    edges[6 * i + 2] = 0;
    edges[6 * i + 3] = 0;		
    edges[6 * i + 4] = 0;		
    edges[6 * i + 5] = 0;

    if (bandwidths.size() == 1) {
      edgeWeights[6 * i] = bandwidths[0];
      edgeWeights[6 * i + 1] = bandwidths[0];		
      edgeWeights[6 * i + 2] = bandwidths[0];		
      edgeWeights[6 * i + 3] = bandwidths[0];
      edgeWeights[6 * i + 4] = bandwidths[0];
      edgeWeights[6 * i + 5] = bandwidths[0];
    } else {
      edgeWeights[6 * i] = bandwidths[elementOfWeightVector(size, startX, startX - 1)];
      edgeWeights[6 * i + 1] = bandwidths[elementOfWeightVector(size, startX, startX + 1)];		
      edgeWeights[6 * i + 2] = bandwidths[elementOfWeightVector(size, startY, startY - 1)];		
      edgeWeights[6 * i + 3] = bandwidths[elementOfWeightVector(size, startY, startY + 1)];      		
      edgeWeights[6 * i + 4] = bandwidths[elementOfWeightVector(size, startZ, startZ - 1)];		
      edgeWeights[6 * i + 5] = bandwidths[elementOfWeightVector(size, startZ, startZ + 1)];
    }

  } 
  nodes[size * size * size] = size * size * size * 6;

  //connects node i with his neighbours
  for (int i = 0; i < size * size * size; i++) {
    //left
    if (i % size == 0) {
      edges[6 * i] = i - 1 + size;
    } else {
      edges[6 * i] = i - 1;
	  }
    //right
    if (i % size == size - 1) {
      edges[6 * i + 1] = i + 1 - size;
    } else {
      edges[6 * i + 1] = i + 1;
    }
    //up
    if ((i % (size * size)) / size == 0) {
      edges[6 * i + 2] = i - size + (size * size);
    } else {
      edges[6 * i + 2] = i - size;
    }
    //down
    if ((i % (size * size)) / size == size - 1) {
      edges[6 * i + 3] = i + size - (size * size);
    } else {
      edges[6 * i + 3] = i + size;
    }
    //front
    if (i / (size * size) == 0) {
      edges[6 * i + 4] = i - (size * size) + (size * size * size);
    } else {
      edges[6 * i + 4] = i - (size * size);
    }
    //back
    if (i / (size * size) == size - 1) {
      edges[6 * i + 5] = i + (size * size) - (size * size * size);
    } else {
      edges[6 * i + 5] = i + (size * size);
    }

  }
  
  vector<int> split(power, 8);
  proc.splitDeg = split;

  proc.build_from_metis_weighted(size * size * size, &nodes[0], &edges[0], &nodeWeights[0], &edgeWeights[0]);
  
  int nodesNo = proc.number_of_nodes();
  for (int i = 0; i < nodesNo; i++) {
    vector<int> row(nodesNo);
    proc.distanceMatrix.push_back(row);
    proc.predecessorMatrix.push_back(row);
  }

  for (int i = 0; i < nodesNo; i++) {
    for (int j = 0; j < nodesNo; j++) {
      int startX = i % size;
      int startY = (i % (size * size)) / size;
      int startZ = i / (size * size);
      int targetX = j % size;
      int targetY = (j % (size * size)) / size;
      int targetZ = j / (size * size);
      bool jumpX = false;
      bool jumpY = false;
      bool jumpZ = false;
      
      int distX, distY, distZ;

      if (bandwidths.size() == 1) {
	distX = std::max(startX, targetX) - std::min(startX, targetX);
	distY = std::max(startY, targetY) - std::min(startY, targetY);
	distZ = std::max(startZ, targetZ) - std::min(startZ, targetZ);

	if (distX * 2 >= size) {
	  distX = size - distX;
	  jumpX = true;
	}
	if (distY * 2 >= size) {
	  distY = size - distY;
	  jumpY = true;
	}
	if (distZ * 2 >= size) {
	  distZ = size - distZ;
	  jumpZ = true;
	}
      
      	distX *= bandwidths[0];
	distY *= bandwidths[0];
	distZ *= bandwidths[0];

      } else {

	int maxX = max(startX, targetX);
	int minX = min(startX, targetX);
	int maxY = max(startY, targetY);
	int minY = min(startY, targetY);
	int maxZ = max(startZ, targetZ);
	int minZ = min(startZ, targetZ);
	
	if ((maxY - minY) * 2 >= size) {
	  jumpY = true;
	} 
	if ((maxX - minX) * 2 >= size) {
	  jumpX = true;
	}
	if ((maxZ - minZ) * 2 >= size) {
	  jumpZ = true;
	}

	int xPath = maxX;
	distX = 0;
	while (xPath != minX) {
	  if (!jumpX) {
	    distX += bandwidths[elementOfWeightVector(size, xPath, xPath - 1)];
	    xPath--;
	  } else {
	    distX += bandwidths[elementOfWeightVector(size, xPath, xPath + 1)];
	    if (xPath + 1 == size) {
	      xPath = 0;
	    } else {
	      xPath++;
	    }
	  }
	}

	int yPath = maxY;
	distY = 0;
	while (yPath != minY) {
	  if (!jumpY) {
	    distY += bandwidths[elementOfWeightVector(size, yPath, yPath - 1)];
	    yPath--;
	  } else {
	    distY += bandwidths[elementOfWeightVector(size, yPath, yPath + 1)];
	    if (yPath + 1 == size) {
	      yPath = 0;
	    } else {
	      yPath++;
	    }
	  }
	}

	int zPath = maxZ;
	distZ = 0;
	while (zPath != minZ) {
	  if (!jumpZ) {
	    distZ += bandwidths[elementOfWeightVector(size, zPath, zPath - 1)];
	    zPath--;
	  } else {
	    distZ += bandwidths[elementOfWeightVector(size, zPath, zPath + 1)];
	    if (zPath + 1 == size) {
	      zPath = 0;
	    } else {
	      zPath++;
	    }
	  }
	}

      }
           
      proc.distanceMatrix[i][j] = distX + distY + distZ;
            
      if (i == j) {
	proc.predecessorMatrix[i][j] = 0;
      } else {
	double partOfZ = ((double) distZ) / ((double)(proc.distanceMatrix[i][j]));
	double randy = random01();
	if((randy <= partOfZ) && (distZ > 0)) {
	  if ((startZ > targetZ && !jumpZ) || (startZ < targetZ && jumpZ)) {
	    if (j + size * size < size * size * size) {
	      proc.predecessorMatrix[i][j] = j + (size * size);
	    } else {
	      proc.predecessorMatrix[i][j] = j + (size * size) - (size * size * size);
	    }
	  } else {
	    if (j >= size * size) {
	      proc.predecessorMatrix[i][j] = j - (size * size);
	    } else {
	      proc.predecessorMatrix[i][j] = j - (size * size) + (size * size * size);
	    }
	  }
	} else {//in particular distX + distY > 0  
	  double partOfY = ((double) distY) / ((double) distX + distY);
	  randy = random01();
	  if((randy <= partOfY) && (distY > 0)) {
	    if ((startY > targetY && !jumpY) || (startY < targetY && jumpY)) {
	      if ((j % (size * size)) + size < size * size) {
		proc.predecessorMatrix[i][j] = j + size;
	      } else {
		proc.predecessorMatrix[i][j] = j + size - size * size;
	      }
	    } else {
	      if (j % (size * size) >= size) {
		proc.predecessorMatrix[i][j] = j - size;
	      } else {
		proc.predecessorMatrix[i][j] = j - size + (size * size);
	      }
	    }
	  }
	  else {//in particular distX > 0
	    if ((startX > targetX && !jumpX) || (startX < targetX && jumpX)) {
	      if (j % size != size - 1) {
		proc.predecessorMatrix[i][j] = j + 1;
	      } else {
		proc.predecessorMatrix[i][j] = j + 1 - size;
	      }
	    } else {
	      if (j % size > 0) {
		proc.predecessorMatrix[i][j] = j - 1;
	      } else {
		proc.predecessorMatrix[i][j] = j - 1 + size;
	      }
	    }
	  }
	}
      }
    }
  }

  
  vector<int> vec;
  proc.treeSort = vec;
  build3DTreeRec(0, size, size, &(proc.treeSort));
  proc.refine = true;
  proc.levelDist = buildLevDist(8);
}

/******************************************************/
/*                   buildHypercube                   */
/******************************************************/
void buildHypercube (processor_graph_access & P, int dimension) {
  NodeID numberNodes = (NodeID) 1;
  for(int i = 0; i < dimension; i++) {
    numberNodes *= (NodeID) 2;
  }
  EdgeID numberEdges = (EdgeID)(numberNodes) * (EdgeID) dimension;

  P.start_construction((NodeID) numberNodes, (EdgeID) numberEdges);
  EdgeID ne = 0;
  for(unsigned i = 0; i < numberNodes; i++) {
    NodeID u = P.new_node();    
    P.setNodeWeight(u, 1);
    P.setPartitionIndex(u, 0);
    for(int pos = 0; pos < dimension; pos++) {
      NodeID v = (NodeID) (flipPosition((int64_t) pos, (int64_t) u));
      EdgeID e = P.new_edge(u, v);
      ne++,
      P.setEdgeWeight(e, 1);
    }
  }
  P.finish_construction();

  //build distance matrix
  P.distanceMatrix.resize(numberNodes);
  for(NodeID i = 0; i < numberNodes; i++) {
    ((P.distanceMatrix)[i]).resize(numberNodes, 0);
  }
  for(NodeID i = 0; i < numberNodes; i++) {
    for(NodeID j = 0; j < numberNodes; j++) {
      P.distanceMatrix[i][j] = hammingDistance(i, j);
    }
  }

  //build predecessor matrix
  P.predecessorMatrix.resize(numberNodes);
  for(NodeID i = 0; i < numberNodes; i++) {
    ((P.predecessorMatrix)[i]).resize(numberNodes, 0);
  }
  for(NodeID i = 0; i < numberNodes; i++) {
    for(NodeID j = 0; j < numberNodes; j++) {
      if (i == j) {
	P.predecessorMatrix[i][j] = 0;
      } else {
	int64_t exor = i ^ j;
	int64_t currNum1bits = 0;//current number of 1-bits in exor, starting the search from least significant position
	//NEXT: crucial number of 1-bits (note that P.distanceMatrix[i][j] == number of 1-bits in exor)
	int64_t crucialNum1bits = 1 + randomFunction(P.distanceMatrix[i][j]);
	int64_t posFromEnd = 0;//least significant position
	int64_t flippingPosition = posFromEnd;
	while(currNum1bits < crucialNum1bits) {
	  if(exor % 2 == 1) {
	    currNum1bits++;
	  }
	  if(currNum1bits == crucialNum1bits) {
	    flippingPosition = posFromEnd;
	  }
	  exor = exor >> 1;
	  posFromEnd++;
	}
	P.predecessorMatrix[i][j] = flipPosition(flippingPosition, j);
      }
    }
  }
}
