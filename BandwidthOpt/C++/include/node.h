#ifndef NODE_H
#define NODE_H
#include <vector>
#include <string>
class Node{
public:
	// construction
	Node() :index_(0), g_(0), min_(INT_MAX), max_(0), degree_(0), degree_total_(0) {}
	// using the address to check for the equity.
	bool operator== (Node& a) { return  index_ == a.index_; }
	bool operator< (Node& a) { return this->degree_ < a.degree_; }

	// add neighbor_total_
	void addNeighborTotal(Node* a);
	// add neighbor_
	void addNeighbor(Node* a);
	// assign the proper g_ 
	int assignG(std::vector<bool>& assign);
	// update, could only be used after construction.
	void update();
	// ±Í ∂∫≈
	int index_;
	// the value of g(v)
	int g_;

	int min_;// the minimun index of the neigbor
	int max_;// the maximum index of the neigbor

	// the total degree
	int degree_total_;
	// the specific degree used in different phase
	int degree_;
	// N(v) for the whole graph
	std::vector<Node*> neighbors_total_;
	// N(v) for part of the graph
	std::vector<Node*> neighbors_;
};
#endif