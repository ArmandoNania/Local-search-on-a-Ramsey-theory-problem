#include <stdio.h>		// printf
#include <stdlib.h>		// strtol, random
#include <time.h>		// time

int compute_dimension(const char *argv);

int **compute_vertices(int verticesCount, int dimension);

int **compute_4_coplanar_points(int *fourCoplanarVerticesCount, int **vertices, int verticesCount, int dimension);

int **compute_subgraph_edges_list(int **fourCoplanarVertices, int fourCoplanarVerticesCount, int verticesCount);

int compute_monocromy_count(int **subgraphEdgesList, int fourCoplanarVerticesCount, char coloredEdges[]);

void hill_climbing(char coloredEdges[], int edgesCount, int **subgraphEdgesList, int fourCoplanarVerticesCount);

int main(int argc, char *argv[]) { 
	// INPUT CHECK
	if (argc!=2) {
		fprintf(stderr, "Usage: %s <integer in [2, 7]>\n", argv[0]);
		return 1;
	}
	//hill_climbing(NULL, 2, NULL, 0);
	//return 0;
	int dimension = compute_dimension(argv[1]);
	if (dimension==-1) return 1;
	printf("Dimension: %i\n", dimension);
	int i;

	// VERTEX CALCULUS
	int verticesCount = 1;
	for (i=0; i<dimension; i++) {
		verticesCount *= 2;
	}
	int **vertices = compute_vertices(verticesCount, dimension);
	printf("Vertex number: %i\n", verticesCount);

	// COPLANAR VERTEX CALCULUS
	int fourCoplanarVerticesCount = 0;
	int **fourCoplanarVertices = compute_4_coplanar_points(&fourCoplanarVerticesCount, vertices, verticesCount, dimension);
	printf("Combinations of 4 coplanar vertex: %i\n", fourCoplanarVerticesCount);

	// LIST OF EDGES OF ALL COMPLETE SUBGRAPHS ON FOUR COPLANAR VERTICES CALCULUS
	int **subgraphEdgesList = compute_subgraph_edges_list(fourCoplanarVertices, fourCoplanarVerticesCount, verticesCount);
	printf("Number of constraints: %i\n", fourCoplanarVerticesCount);

	// UNNECESSARY DATA MEMORY RELEASE
    for (i = 0; i < verticesCount; ++i) {
        free(vertices[i]);
    }
    free(vertices);
    for (i = 0; i < fourCoplanarVerticesCount; ++i) {
        free(fourCoplanarVertices[i]);
    }
    free(fourCoplanarVertices);

    // EDGES CALCULUS
    int edgesCount = verticesCount*(verticesCount-1)/2;
    char coloredEdges[edgesCount];
    unsigned int seed = time(NULL);
    srand(seed);
    for (i=0; i<edgesCount; i++) {
    	coloredEdges[i] = rand() % 2;
    }
	printf("Number of edges: %i\n", edgesCount);
	printf("Seed: %u\n", seed);

	// MONOCROMY COUNT CALCULUS
	hill_climbing(coloredEdges, edgesCount, subgraphEdgesList, fourCoplanarVerticesCount);

    // END
    for (i = 0; i < fourCoplanarVerticesCount; ++i) {
        free(subgraphEdgesList[i]);
    }
    free(subgraphEdgesList);
    printf("\nDone.\n");
	return 0;
}

int compute_dimension(const char *argv) {
	long arg = strtol(argv, NULL, 10);
	if ((arg<2)||(arg>7)) {
		fprintf(stderr, "The number must be in the interval [2, 7]\n");
		return -1;
	}
	int dimension = (int) arg;
}

int **compute_vertices(int verticesCount, int dimension) {
	int i, decimalNumber = 0;
	int **vertices = (int **)malloc(verticesCount * sizeof(int *));
	for (i=0; i<verticesCount; i++) {
		vertices[i] = (int *)malloc(dimension * sizeof(int));
	}
	for (i=0; i<verticesCount; i++) {
		decimalNumber = i;
		for (int j=dimension-1; j>=0; j--) {
			vertices[i][j] = decimalNumber % 2;
			decimalNumber /= 2;
		}
	}
	return vertices;
}

void compute_4_points(int verticesCount,
					  int **fourVertices,
					  int vertex,
					  int *fourVerticesIndex1,
					  int fourVerticesIndex2,
					  int fourVerticesCount) {

	if (fourVerticesIndex2==4) {
		(*fourVerticesIndex1)++;
		if (*fourVerticesIndex1<fourVerticesCount) {
			for (int i=0; i<3; i++) {
				fourVertices[*fourVerticesIndex1][i] = fourVertices[*fourVerticesIndex1-1][i];
			}

		}
		return;
	}

	for (int i=vertex; i<=verticesCount-4+fourVerticesIndex2; i++) {
		fourVertices[*fourVerticesIndex1][fourVerticesIndex2] = i;
		compute_4_points(verticesCount, fourVertices, i+1, fourVerticesIndex1, fourVerticesIndex2+1, fourVerticesCount);
	}

}

int compute_coplanarity(char *flags,
						int **vertices,
						int dimension,
						int **fourVertices,
						int fourVerticesCount) {

	/*
	* Suppose we have the system of equations: Ax1 + By1 = z1; Ax2 + By2 = z2
	* in this case the solution is A = (y2z1 - y1z2)/D; B = (x1z2 - x2z1)/D;
	* with D = x1y2 - y1x2. So I check if ADx + BDy = Dz for each dimension in the matrix.
	* If that is the case, it means that the vectors matrix[0], matrix[1], matrix[2] are
	* linearly dependent, which means that the 4 points from which they derive are coplanar.
	*/

	int fourCoplanarVerticesCount = 0;
	int matrix[3][dimension];
	int i, j, k;
	int parameter_AD, parameter_BD, parameter_D;
	char flag;
	for (i=0; i<fourVerticesCount; i++) {
		// compute matrix of vectors
		for (j=0; j<dimension; j++) {
			matrix[0][j] = vertices[ fourVertices[i][0] ][ j ] - vertices[ fourVertices[i][3] ][ j ];
			matrix[1][j] = vertices[ fourVertices[i][1] ][ j ] - vertices[ fourVertices[i][3] ][ j ];
			matrix[2][j] = vertices[ fourVertices[i][2] ][ j ] - vertices[ fourVertices[i][3] ][ j ];
		}
		// find a non-zero parameter_D
		for (j=0; j<dimension-1; j++) {
			for (k=j+1; k<dimension; k++) {
				parameter_D = (matrix[0][j] * matrix[1][k]) - (matrix[0][k] * matrix[1][j]);
				if (parameter_D != 0) break;
			}
			if (parameter_D != 0) break;
		}
		if (parameter_D == 0) {
			// this should never be executed
			printf("---------------------------------------------------------------\n");
			printf("\tSome error occurred in the data generation!\n\tCheck the matrix variable, when it comes from vertices[%i]", i);
			printf("---------------------------------------------------------------\n");
			continue;
		}
		// compute other parameters
		parameter_AD = (matrix[1][k] * matrix[2][j]) - (matrix[1][j] * matrix[2][k]);
		parameter_BD = (matrix[0][j] * matrix[2][k]) - (matrix[0][k] * matrix[2][j]);
		// check linear dependence
		flag = 1;
		for (j=0; j<dimension; j++) {
			if ( ((parameter_AD*matrix[0][j]) + (parameter_BD*matrix[1][j])) != (parameter_D*matrix[2][j]) ) {
				flag = 0;
				break;
			}
		}
		// save detected coplanar points
		if (flag) {
			flags[i] = 1;
			fourCoplanarVerticesCount++;
		} else {
			flags[i] = 0;
		}
	}

	return fourCoplanarVerticesCount;
}

int **compute_4_coplanar_points(int *fourCoplanarVerticesCount, int **vertices, int verticesCount, int dimension) {
	// compute the number of combinations of 4 elements in a set of verticesCount elements
	int fourVerticesCount = 1;
	for (int i=verticesCount-3; i<=verticesCount; i++)
	{
		fourVerticesCount *= i;
	}
	fourVerticesCount /= 24;

	// compute the combinations of 4 elements in a set of verticesCount elements
	int **fourVertices = (int **)malloc(fourVerticesCount * sizeof(int *));
	for (int i=0; i<fourVerticesCount; i++) {
		fourVertices[i] = (int *)malloc(4 * sizeof(int));
	}
	int fourVerticesIndex1 = 0;
	compute_4_points(verticesCount, fourVertices, 0, &fourVerticesIndex1, 0, fourVerticesCount);
	int intCheck;
	char boolCheck = 0;
	for (int i=0; i<fourVerticesCount; i++) {
		intCheck = fourVertices[i][0];
		for (int j=1; j<4; j++) {
			if (fourVertices[i][j]<=intCheck) {
				boolCheck = 1;
				break;
			}
			intCheck = fourVertices[i][j];
		}
		if (boolCheck) {
			// this should never be executed
			printf("---------------------------------------------------------------\n");
			printf("\tSome error occurred in the data generation!\n\tCheck the fourVertices variable, at fourVertices[%i]\n", i);
			printf("---------------------------------------------------------------\n");
			break;
		}
	}

	// compute the combinations of 4 coplanar points
	char *flags = (char *)malloc(fourVerticesCount * sizeof(char));
	*fourCoplanarVerticesCount = compute_coplanarity(flags, vertices, dimension, fourVertices, fourVerticesCount);

	// store the result in a new variable
	int **result = (int **)malloc(*fourCoplanarVerticesCount * sizeof(int *));
	for (int i=0; i<*fourCoplanarVerticesCount; i++) {
		result[i] = (int *)malloc(4 * sizeof(int));
	}
	int resultIndex = 0;
	for (int i=0; i<fourVerticesCount; i++) {
		if (flags[i]) {
			for (int j=0; j<4; j++) {
				result[resultIndex][j] = fourVertices[i][j];
			}
			resultIndex++;
		}
	}

	// release the memory allocated & return
	free(flags);
    for (int i = 0; i < fourVerticesCount; ++i) {
        free(fourVertices[i]);
    }
    free(fourVertices);
    return result;
}

int **compute_subgraph_edges_list(int **fourCoplanarVertices, int fourCoplanarVerticesCount, int verticesCount) {
	int i, edgeIndex, subgraphEdgesIndex;
	int **subgraphEdgesList = (int **)malloc(fourCoplanarVerticesCount * sizeof(int *));
	for (i=0; i<fourCoplanarVerticesCount; i++) {
		subgraphEdgesList[i] = (int *)malloc(6 * sizeof(int));
	}
	for (i=0; i<fourCoplanarVerticesCount; i++){
		subgraphEdgesIndex = 0;
		for (int j=0; j<3; j++){
			for (int k=j+1; k<4; k++){
				edgeIndex = 0;
				for (int l=0; l<fourCoplanarVertices[i][j]; l++) {
					edgeIndex += (verticesCount-1-l);
				}
				edgeIndex += (fourCoplanarVertices[i][k]-1-fourCoplanarVertices[i][j]);
				subgraphEdgesList[i][subgraphEdgesIndex] = edgeIndex;
				subgraphEdgesIndex++;
			}
		}
	}
	return subgraphEdgesList;
}

int compute_monocromy_count(int **subgraphEdgesList, int fourCoplanarVerticesCount, char coloredEdges[]) {
	int monocromyCount = 0;
	char monocromyCheck, color;
	for (int i=0; i<fourCoplanarVerticesCount; i++) {
		color = coloredEdges[ subgraphEdgesList[i][0] ];
		monocromyCheck = 1;
		for (int j=1; j<6; j++) {
			if (coloredEdges[ subgraphEdgesList[i][j] ] != color) {
				monocromyCheck = 0;
				break;
			}
		}
		if (monocromyCheck) monocromyCount++;
	}
	return monocromyCount;
}

void swtich_color(char coloredEdges[], int index) {
	if (coloredEdges[index]==0) {
	    coloredEdges[index] = 1;
	} else {
	    coloredEdges[index] = 0;
	}
}

void hill_climbing(char coloredEdges[], int edgesCount, int **subgraphEdgesList, int fourCoplanarVerticesCount) {
	int prevChildIndex = rand() % edgesCount;
	int monocromyCount;
	int bestChildIndex;
	int bestChildScore;
	int absolutBest = fourCoplanarVerticesCount + 1;
    bestChildIndex = 0;
	for (int i=0; i<100; i++) {
		bestChildScore = fourCoplanarVerticesCount + 1;
	    for (int j=0; j<edgesCount; j++) {
	    	if (j==prevChildIndex) continue;
	    	swtich_color(coloredEdges, j);
		    monocromyCount = compute_monocromy_count(subgraphEdgesList, fourCoplanarVerticesCount, coloredEdges);
		    if (monocromyCount == 0) {
				printf("%i%%: Found a counter-exemple for the Graham Problem!!!\n", i+1);
		    	for (int j=0; j<edgesCount; j++) {
		    		if (coloredEdges[j]) {
		    			printf("red ");
		    		} else {
		    			printf("blue ");
		    		}
		    	}
		    	return;
		    } else if (monocromyCount<bestChildScore) {
		    	bestChildScore = monocromyCount;
		    	bestChildIndex = j;
		    }
			if (monocromyCount<absolutBest) {
		    	absolutBest = monocromyCount;
		    }
			swtich_color(coloredEdges, j);
		}
		if (i<10) printf(" ");
		printf("%i%%: %i\n", i+1, absolutBest);
		swtich_color(coloredEdges, bestChildIndex);
		prevChildIndex = bestChildIndex;
	}
  	printf("Counter-exemple not found. The coloring with the least number of\nsingle-coloured complete subgraphs had %i such subgraphs", absolutBest);
}