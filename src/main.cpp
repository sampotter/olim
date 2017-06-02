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

	grid2d g(height, width, 1, F);
	g.add_boundary_node(0, 0);
	g.add_boundary_node(19, 10);
	g.add_boundary_node(10, 15);
	g.run();

	printf("[");
	for (size_t i = 0; i < height - 1; ++i) {
		for (size_t j = 0; j < width; ++j) {
			printf("%g ", g(j, i).get_value());
		}
		printf(";");
	}
	for (size_t j = 0; j < width; ++j) {
		printf("%g ", g(j, height - 1).get_value());
	}
	printf("];\n");
}
