#ifndef __NODE_HPP__
#define __NODE_HPP__

#include <limits>

enum class state {valid, trial, far};

struct node {
	static node make_boundary_node(size_t i, size_t j);

	double get_value() const;
	void set_value(double value);

	size_t get_i() const;
	void set_i(size_t i);
	
	size_t get_j() const;
	void set_j(size_t j);

	size_t get_heap_pos() const;
	void set_heap_pos(size_t pos);

	bool is_valid() const;
	bool is_trial() const;
	bool is_far() const;
	void set_valid();
	void set_trial();
	void set_far();
private:
	double _value {std::numeric_limits<double>::infinity()};
	state _state {state::far};
	size_t _i;
	size_t _j;
	size_t _heap_pos;
};

#endif // __NODE_HPP__
