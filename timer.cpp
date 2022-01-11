#include "timer.h"

timer::timer(const char* msg) : message(msg) {
  tstart = std::chrono::high_resolution_clock::now();
}

double timer::elapsed() {
  return std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::high_resolution_clock::now() - tstart).count();
}

void timer::reset() {
  tstart = std::chrono::high_resolution_clock::now();
}

timer::~timer() {
  if(!message.empty())
    std::cout << message << ' ' << elapsed() << " s\n";
}
