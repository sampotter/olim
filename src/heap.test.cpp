#define BOOST_TEST_MODULE heap

#include <boost/test/included/unit_test.hpp>

#include "heap.hpp"

BOOST_AUTO_TEST_CASE (empty_works) {
  heap h {16};
  BOOST_TEST(h.empty());
}

BOOST_AUTO_TEST_CASE (size_works) {
  heap h {16};
  node n[] = {
    node::make_boundary_node(0, 0, 1),
    node::make_boundary_node(0, 0, 2),
    node::make_boundary_node(0, 0, 3)
  };
  BOOST_TEST(h.size() == 0);
  h.insert(&n[1]);
  h.insert(&n[0]);
  h.insert(&n[2]);
  BOOST_TEST(h.size() == 3);
  h.pop_front();
  BOOST_TEST(h.size() == 2);
  h.pop_front();
  h.pop_front();
  BOOST_TEST(h.size() == 0);
}

BOOST_AUTO_TEST_CASE (insert_works) {
  heap h {16};
  node n[] = {
    node::make_boundary_node(0, 0, 1),
    node::make_boundary_node(0, 0, 2),
    node::make_boundary_node(0, 0, 3)
  };
  h.insert(&n[2]);
  h.insert(&n[1]);
  h.insert(&n[0]);
  node * front = h.front();
  BOOST_TEST(front == &n[0]);
  BOOST_TEST(front->get_value() == 1);
  h.pop_front();
  front = h.front();
  BOOST_TEST(front == &n[1]);
  BOOST_TEST(front->get_value() == 2);
  h.pop_front();
  front = h.front();
  BOOST_TEST(front == &n[2]);
  BOOST_TEST(front->get_value() == 3);
  h.pop_front();
}

BOOST_AUTO_TEST_CASE (adjust_entry_works) {
  heap h {16};
  node n[] = {
    node::make_boundary_node(0, 0, 1),
    node::make_boundary_node(0, 0, 2),
    node::make_boundary_node(0, 0, 3)
  };
  h.insert(&n[0]);
  h.insert(&n[1]);
  h.insert(&n[2]);
  n[2].set_value(0);
  h.adjust_entry(&n[2]);
  node * front = h.front();
  BOOST_TEST(front == &n[2]);
  BOOST_TEST(front->get_value() == 0);
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
