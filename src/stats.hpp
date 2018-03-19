#ifndef __STATS_HPP_HPP__
#define __STATS_HPP_HPP__

enum class line_update: int {
  DEG_1,
  DEG_2,
  DEG_3,
  NUM
};

enum class tri_update: int {
  DEG_11,
  DEG_12,
  DEG_13,
  DEG_22,
  DEG_23,
  NUM
};

enum class tetra_update: int {
  DEG_111,
  DEG_112,
  DEG_113,
  DEG_122,
  DEG_123,
  DEG_222,
  DEG_223,
  NUM
};

struct olim3d_node_stats {
  olim3d_node_stats();
  int num_line_updates() const;
  int num_line_updates(int d) const;
  void inc_line_updates(int d);
  int num_tri_updates() const;
  int num_degenerate_tri_updates() const;
  void inc_tri_updates(int d1, int d2, bool degenerate);
  int num_tetra_updates() const;
  void inc_tetra_updates(int d1, int d2, int d3);
private:
  int _num_line_updates[static_cast<int>(line_update::NUM)];
  int _num_tri_updates[static_cast<int>(tri_update::NUM)];
  int _num_degenerate_tri_updates[static_cast<int>(tri_update::NUM)];
  int _num_tetra_updates[static_cast<int>(tetra_update::NUM)];
};

#endif // __STATS_HPP_HPP__
