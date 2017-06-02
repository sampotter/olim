#include "grid.hpp"

#include <cstdio>

double F(double x, double y) {
	(void) x;
	(void) y;
	return 1.0;
}

int main() {
	size_t height = 20;
	size_t width = 20;

	fast_marcher m(height, width, 1, F);
	m.add_boundary_node(0, 0);
	m.add_boundary_node(19, 10);
	m.add_boundary_node(10, 15);
	m.run();

	printf("[");
	for (size_t i = 0; i < height - 1; ++i) {
		for (size_t j = 0; j < width; ++j) {
			printf("%g ", m(j, i).get_value());
		}
		printf(";");
	}
	for (size_t j = 0; j < width; ++j) {
		printf("%g ", m(j, height - 1).get_value());
	}
	printf("];\n");
}
