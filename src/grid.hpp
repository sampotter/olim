#ifndef __GRID_HPP__
#define __GRID_HPP__

#include "heap.hpp"
#include "node.hpp"

using speed_func = double(*)(double, double);

struct fast_marcher
{
	fast_marcher(size_t height, size_t width, double h, speed_func F);
	~fast_marcher();

	node & operator()(size_t i, size_t j);
	node const & operator()(size_t i, size_t j) const;

	void add_boundary_node(size_t i, size_t j);
	void run();
	double get_value(size_t i, size_t j) const;
    
private:
	void update_node_value(size_t i, size_t j);
	void update_neighbors(size_t i, size_t j);
	int get_far_neighbors(size_t i, size_t j, node** nb) const;
	void get_valid_neighbors(size_t i, size_t j, node** nb) const;
	bool valid_index(size_t i, size_t j) const;
	node* get_next_node();
	
	node* _nodes;
	heap _heap;
	double _h;
	speed_func _F;
	size_t _height;
	size_t _width;
};

extern "C"
void fmm_mex(double * out, bool * in, size_t M, size_t N, double h,
			 speed_func F);

#endif // __GRID_HPP__

// Local Variables:
// indent-tabs-mode: nil
// End:
