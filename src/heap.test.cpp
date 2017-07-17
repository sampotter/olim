#include "heap.hpp"
#include "node.hpp"
#include "test.hpp"

void empty_works() {
  heap<node> h {16};
  test::is_true(h.empty());
}

void size_works() {
  heap<node> h {16};
  node n[] = {{0, 0, 1}, {0, 0, 2}, {0, 0, 3}};
  test::is_true(h.size() == 0);
  h.insert(&n[1]);
  h.insert(&n[0]);
  h.insert(&n[2]);
  test::is_true(h.size() == 3);
  h.pop_front();
  test::is_true(h.size() == 2);
  h.pop_front();
  h.pop_front();
  test::is_true(h.size() == 0);
}

void insert_works() {
  heap<node> h {16};
  node n[] = {{0, 0, 1}, {0, 0, 2}, {0, 0, 3}};
  h.insert(&n[2]);
  h.insert(&n[1]);
  h.insert(&n[0]);
  node * front = h.front();
  test::is_true(front == &n[0]);
  test::is_true(front->get_value() == 1);
  h.pop_front();
  front = h.front();
  test::is_true(front == &n[1]);
  test::is_true(front->get_value() == 2);
  h.pop_front();
  front = h.front();
  test::is_true(front == &n[2]);
  test::is_true(front->get_value() == 3);
  h.pop_front();
}

void swim_works() {
  heap<node> h {16};
  node n[] = {{0, 0, 1}, {0, 0, 2}, {0, 0, 3}};
  h.insert(&n[0]);
  h.insert(&n[1]);
  h.insert(&n[2]);
  n[2].set_value(0);
  h.swim(&n[2]);
  node * front = h.front();
  test::is_true(front == &n[2]);
  test::is_true(front->get_value() == 0);
}

int main() {
  empty_works();
  size_works();
  insert_works();
  swim_works();
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
