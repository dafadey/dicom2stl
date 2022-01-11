#pragma once

#include <iostream>
#include <ctime>
#include <ratio>
#include <chrono>
#include <string>

struct timer {
  timer() = delete;
  timer(const char* msg);

  ~timer();

  void reset();
  double elapsed();

  std::string message;
  
  std::chrono::high_resolution_clock::time_point tstart;

};
