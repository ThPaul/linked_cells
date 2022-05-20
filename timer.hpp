#pragma once


#include <chrono>


using Duration = std::chrono::duration<double, std::micro>;
using Timer = std::chrono::high_resolution_clock;

double to_ms(Duration d) {
  return std::chrono::duration_cast<std::chrono::microseconds>(d).count()/1000.;
}
