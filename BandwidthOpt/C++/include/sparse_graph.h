#ifndef SPARSE_GRAPH_H
#define SPARSE_GRAPH_H
#include "node.h"

class SparseGraph {
public:
	SparseGraph(std::string filename);

	void updateBandwidth() {
		bandwidth_ = 1;
		for (Node* node : nodes_) {
			node->update();
			if (node->g_ - node->min_ + 1 > bandwidth_) {
				bandwidth_ = node->g_ - node->min_ + 1;
				nbwdth_ = 1;
			}
			else if (node->g_ - node->min_ + 1 == bandwidth_) {
				nbwdth_++;
			}
		}
	}
	// GRASP algorithm
	// Construction phase
	void C1();

	// Improvement phase
	int checkMove(Node* u, Node* v); //return the change in the number of  which reach the maximum bandwidth while moving.
	void move(Node* u, Node* v) {
		int temp = u->g_;
		u->g_ = v->g_;
		v->g_ = temp;
		updateBandwidth();
	}

	void localSearch();

private:
	std::vector<Node*> nodes_;
	// 注意一下两者的值由第一次updateBandwidth的时候赋
	int bandwidth_;
	int nbwdth_; // number of  with maximum bandwidth
};
#endif
