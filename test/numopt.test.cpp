#include "numopt.hpp"
#include "test.hpp"

using namespace arma;

void test_qpez() {
  {
    mat G = {{2, 1}, {1, 2}};
    vec c = {3, 2};
    mat A; A.eye(2, 2);
    vec x = optim::qpez(G, c, A);
    IS_TRUE(all(abs(x) < 1e-15));
  }
  {
    mat G(2, 2); G.eye(2, 2);
    vec c = {1, 1};
    mat A = {{-1, 1}};
    vec x = optim::qpez(G, c, A);
    IS_APPROX_EQUAL(x(0), -1.0);
    IS_APPROX_EQUAL(x(1), -1.0);
  }
}

int main() {
  test_qpez();
}
