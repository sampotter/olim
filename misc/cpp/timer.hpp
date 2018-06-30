#include <chrono>
#include <limits>

#ifndef REPS
#    define REPS 5
#endif

using namespace std::chrono;

system_clock::time_point t0;

void tic() {
  t0 = system_clock::now();
}

double toc() {
  return duration<double, std::milli>(system_clock::now() - t0).count()/1000;
}

double time(std::function<void()> thunk, int reps = REPS) {
  double t = std::numeric_limits<double>::infinity();
  for (int i = 0; i < reps; ++i) {
    tic();
    thunk();
    t = std::min(t, toc());
  }
  return t;
}
