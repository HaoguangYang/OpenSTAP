#ifndef NODE_CPP
#define NODE_CPP
#include "node.h"

// add neighbor_total_
void Node::addNeighborTotal(Node* a) {
	if (this != a && find(neighbors_total_.begin(), neighbors_total_.end(), a) == neighbors_total_.end()) {
		neighbors_total_.push_back(a);
		degree_total_++;
	}
}
// add neighbor_
void Node::addNeighbor(Node* a) {
	if (this != a && find(neighbors_.begin(), neighbors_.end(), a) == neighbors_.end()) {
		neighbors_.push_back(a);
	}
}
// update
void Node::update() {
	min_ = INT_MAX;
	for(auto node : neighbors_total_) {
		if (min_ > node->g_) {
			min_ = node->g_;
		}
	}
}

int Node::assignG(std::vector<bool>& assign) {
	for (auto node : neighbors_) {
		if (min_ > node->g_) {
			min_ = node->g_;
		}
		if (max_ < node->g_) {
			max_ = node->g_;
		}
	}
	int temp = (min_ + max_) / 2;
	if (!assign[temp]) {
		g_ = temp;
		assign[temp] = true;
		return temp;
	}
	int i = 1;
	while (true) {
		if (temp - i > 0 && !assign[temp - i]) {
			g_ = temp - i;
			assign[temp - i] = true;
			return temp;
		}
		if (temp + i < assign.size() && !assign[temp + i]) {
			g_ = temp + i;
			assign[temp + i] = true;
			return temp;
		}
		i++;
	}
}

#endif