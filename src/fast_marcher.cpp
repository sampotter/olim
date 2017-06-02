#include "fast_marcher.hpp"

#include <cassert>
#include <cmath>

fast_marcher::fast_marcher(size_t height, size_t width, double h, speed_func F):
	_nodes {new node[height*width]},
	_heap {static_cast<size_t>(std::log(height*width))}, // whatever
	_h {h},
	_F {F},
	_height {height},
	_width {width}
{
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
	update_node_value_impl(i, j);
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

double fast_marcher::get_h() const {
	return _h;
}

double fast_marcher::F(double x, double y) const {
	return _F(x, y);
}

double fast_marcher::F(size_t i, size_t j) const {
	return _F(_h*i, _h*j);
}

void fast_marcher::adjust_heap_entry(node* n) {
	_heap.adjust_entry(n);
}
