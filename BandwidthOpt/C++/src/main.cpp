#include <fstream>
#include <iostream>
#include <algorithm> // for finding
#include "sparse_graph.h"



int main() {
	SparseGraph sg("sparse_cpp.in");
	// construction phase
	sg.C1();
	sg.localSearch();
	return 0;
}