#include "grid.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>

fast_marcher::fast_marcher(size_t height, size_t width, double h, speed_func F):
	_nodes {new node[height*width]},
	_heap {static_cast<size_t>(log(height*width))}, // whatever
	_h {h},
	_F {F},
	_height {height},
	_width {width}
{
	std::cout << _h << std::endl;
	
	for (size_t i = 0; i < _height; ++i) {
		for (size_t j = 0; j < _width; ++j) {
			this->operator()(i, j).set_i(i);
			this->operator()(i, j).set_j(j);
		}
	}
}

fast_marcher::~fast_marcher() {
	delete[] _nodes;
}

node & fast_marcher::operator()(size_t i, size_t j) {
	return _nodes[_width*i + j];
}

node const & fast_marcher::operator()(size_t i, size_t j) const {
	return _nodes[_width*i + j];
}

void fast_marcher::add_boundary_node(size_t i, size_t j) {
	this->operator()(i, j) = node::make_boundary_node(i, j);
	update_neighbors(i, j);
}

void fast_marcher::update_neighbors(size_t i, size_t j) {
	static int offsets[4][2] = {{-1, 0}, {0, -1}, {0, 1}, {1, 0}};
	size_t a, b;

	// Change all 'far' neighbors to 'trial'
	node * n = nullptr;
	for (int k = 0; k < 4; ++k) {
		a = i + offsets[k][0], b = j + offsets[k][1];
		if (valid_index(a, b) && this->operator()(a, b).is_far()) {
			n = &_nodes[_width*a + b];
			n->set_trial();
			_heap.insert(n);
		}
	}

	for (int k = 0; k < 4; ++k) {
		a = i + offsets[k][0], b = j + offsets[k][1];
		if (valid_index(a, b) && !this->operator()(a, b).is_valid()) {
			update_node_value(a, b);
		}
	}
}

void fast_marcher::run() {
	node* n = nullptr;
	while (!_heap.empty()) {
		n = get_next_node();
		n->set_valid();
		update_neighbors(n->get_i(), n->get_j());
	}
}

double fast_marcher::get_value(size_t i, size_t j) const {
	return this->operator()(i, j).get_value();
}

void fast_marcher::update_node_value(size_t i, size_t j) {
	node* n = nullptr;
	node* nb[4] = {nullptr, nullptr, nullptr, nullptr}; // NESW
	get_valid_neighbors(i, j, nb);
	double rhs = _h/_F(_h*i, _h*j);
	double T = std::numeric_limits<double>::infinity();
	double tmp = 0, T1 = 0, T2 = 0, disc = 0;

	for (int k = 0, k1 = 1; k < 4; ++k, k1 = (k1 + 1) % 4) {
		if (nb[k] != nullptr && nb[k1] != nullptr) {
			T1 = nb[k]->get_value();
			T2 = nb[k1]->get_value();
			disc = 2*pow(rhs, 2) - pow(T1 - T2, 2);
			if (disc <= 0) continue;
			tmp = (T1 + T2 + std::sqrt(disc))/2;
			if (tmp >= T1 && tmp >= T2) T = std::min(T, tmp);
		} else if (nb[k] != nullptr || nb[k1] != nullptr) {
			n = nb[k] != nullptr ? nb[k] : nb[k1];
			T = std::min(T, n->get_value() + rhs);
		}
	}

	n = &_nodes[i*_width + j];
	assert(n->is_trial());
	n->set_value(T);
	_heap.adjust_entry(n);
}

void fast_marcher::get_valid_neighbors(size_t i, size_t j, node ** nb) const {
	auto const is_good = [this] (size_t a, size_t b) {
		return valid_index(a, b) && this->operator()(a, b).is_valid();
	};

	if (is_good(i - 1, j)) nb[0] = &_nodes[_width*(i - 1) + j];
	if (is_good(i, j + 1)) nb[1] = &_nodes[_width*i + j + 1];
	if (is_good(i + 1, j)) nb[2] = &_nodes[_width*(i + 1) + j];
	if (is_good(i, j - 1)) nb[3] = &_nodes[_width*i + j - 1];
}

int fast_marcher::get_far_neighbors(size_t i, size_t j, node ** nb) const {
	int nnb = 0;

	auto const is_good = [this] (size_t a, size_t b) {
		return valid_index(a, b) && this->operator()(a, b).is_far();
	};
	
	if (is_good(i - 1, j)) { // north
		nb[0] = &_nodes[_width*(i - 1) + j];
		++nnb;
	}
	if (is_good(i, j + 1)) { // east
		nb[1] = &_nodes[_width*i + j + 1];
		++nnb;
	}
	if (is_good(i + 1, j)) { // south
		nb[2] = &_nodes[_width*(i + 1) + j];
		++nnb;
	}
	if (is_good(i, j - 1)) { // west
		nb[3] = &_nodes[_width*i + j - 1];
		++nnb;
	}
	return nnb;
}

bool fast_marcher::valid_index(size_t i, size_t j) const {
	return i < _height && j < _width;
}

node* fast_marcher::get_next_node() {
	auto const elt = _heap.front();
	_heap.pop_front();
	return elt;
}

void fmm_mex(double * out, bool * in, size_t M, size_t N, double h,
			 speed_func F) {
	fast_marcher g(M, N, h, F);
	
	for (size_t i = 0; i < M; ++i) {
		for (size_t j = 0; j < N; ++j) {
			if (in[N*i + j]) {
				g.add_boundary_node(i, j);
			}
		}
	}

	g.run();

	for (size_t i = 0; i < M; ++i) {
		for (size_t j = 0; j < N; ++j) {
			out[M*i + j] = g.get_value(i, j);
		}
	}
}
