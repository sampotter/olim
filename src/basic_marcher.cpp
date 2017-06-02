#include "basic_marcher.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>

void basic_marcher::update_node_value_impl(size_t i, size_t j) {
	node* n = nullptr;
	node* nb[4] = {nullptr, nullptr, nullptr, nullptr}; // NESW
	get_valid_neighbors(i, j, nb);
	double rhs = get_h()/F(i, j);
	double T = std::numeric_limits<double>::infinity();
	double tmp = 0, T1 = 0, T2 = 0, disc = 0;

	for (int k = 0, k1 = 1; k < 4; ++k, k1 = (k1 + 1) % 4) {
		if (nb[k] != nullptr && nb[k1] != nullptr) {
			T1 = nb[k]->get_value();
			T2 = nb[k1]->get_value();
			disc = 2*std::pow(rhs, 2) - std::pow(T1 - T2, 2);
			if (disc <= 0) continue;
			tmp = (T1 + T2 + std::sqrt(disc))/2;
			if (tmp >= T1 && tmp >= T2) T = std::min(T, tmp);
		} else if (nb[k] != nullptr || nb[k1] != nullptr) {
			n = nb[k] != nullptr ? nb[k] : nb[k1];
			T = std::min(T, n->get_value() + rhs);
		}
	}

	n = &this->operator()(i, j);
	assert(n->is_trial());
	n->set_value(T);
	adjust_heap_entry(n);
}
