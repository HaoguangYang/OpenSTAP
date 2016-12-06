#ifndef SPARSE_GRAPH_CPP
#define SPARSE_GRAPH_CPP
#include "sparse_graph.h"
#include <fstream>
#include <iostream>
#include <stdlib.h>
SparseGraph::SparseGraph(std::string filename){
	std::ifstream in(filename);
	if (in.is_open()) {
		//���벿��
		std::cout << "Input phase" << std::endl;
		int nnp;
		in >> nnp;
		for (int i = 0; i < nnp + 1; i++) {
			nodes_.push_back(new Node());
			nodes_[i]->index_ = i;
		}
		int nele;
		in >> nele;
		int npele; // node per element
		for (int i = 0; i < nele; i++) {
			in >> npele;
			std::vector<int> temp(npele);
			for (int j = 0; j < npele; j++) {
				in >> temp[j];
			}
			for (int j = 0; j < npele; j++) {
				for (int k = 0; k < npele; k++) {
					if (k != j) {
						nodes_[temp[j]]->addNeighborTotal(nodes_[temp[k]]);
					}
				}
			}
		}
	}
	else {
		std::cout << "Does not open file!";
	}
}

void SparseGraph::C1() {
	std::vector<bool> sign(nodes_.size(), true);    // �жϽڵ��Ƿ��ܸ�ֵ��
	std::vector<bool> assign(nodes_.size(), false); // ���ĸ�ֵ��������
	int i = rand()%(nodes_.size()-1) + 1;           // ������ɵ�һ��label
	sign[nodes_[1]->index_] = false;
	nodes_[1]->g_ = i;
	assign[i] = true;
	std::vector<Node*> CL;
	for (auto node : nodes_[1]->neighbors_total_) {
		CL.push_back(node);
		node->addNeighbor(nodes_[1]);
	}

	while (!CL.empty()) {
		int max_connectivity = 0;
		Node* max_node = CL[0];
		auto itr = CL.begin();
		for (; itr != CL.end(); itr++) {
			if (max_connectivity < (*itr)->degree_) {
				max_connectivity = (*itr)->degree_;
				max_node = *itr;
			}
		}
		for (auto node : max_node->neighbors_total_) {
			if (std::find(CL.begin(), CL.end(), node) == CL.end() && sign[node->index_]) {
				CL.push_back(node);
			}
				node->degree_++;
				node->addNeighbor(max_node);
		}
		//max_node->g_ = i++;
		max_node->assignG(assign);
		itr = std::find(CL.begin(), CL.end(), max_node);
		int tem = itr - CL.begin() + 1;
		sign[max_node->index_] = false;
		CL.erase(itr);
	}
	// ����֮�����Bandwidth
	updateBandwidth();
}

int SparseGraph::checkMove(Node* u, Node* v) {
	if (u->index_ == v->index_ || u->index_ == 0 || v->index_ == 0)
		return 0;
	// ע�⣬���ڽڵ����Ų����ظ������Դﵽ�ٽ�Ķ���ֻ������һ��neighbor����֮��Ϊ���
	// ���ֻ��Ҫ����������ڵ�����ڽڵ�֮��Ĺ�ϵ����
	int before = 0;
	int after  = 0;
	int gu = u->g_;
	int gv = v->g_;
	for (Node* node : u->neighbors_total_) {
		if (abs(node->g_ - gu) == bandwidth_-1) {
			before++;
		}
		if (abs(node->g_ - gv) >= bandwidth_-1) {
			after++;
		}
	}
	for (Node* node : v->neighbors_total_) {
		if (abs(node->g_ - gv) == bandwidth_-1) {
			before++;
		}
		if (abs(node->g_ - gu) >= bandwidth_-1) {
			after++;
		}
	}
	return after - before;
}

void SparseGraph::localSearch() {
	bool change = true;
	while (change) {
		change = false;
		for (Node* u : nodes_) {
			for (Node* v : nodes_) {
				if (checkMove(u, v) < 0) {
					move(u, v);
					change = true;
				}
			}
		}
	}
}
#endif