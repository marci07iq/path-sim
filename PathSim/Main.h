#pragma once

#include "Point.h"
#include <sstream>
#include <cmath>
#include <algorithm>
#include <chrono>

struct movement {
  dVec3 pos;
  dVec3 vel;
  dVec3 acc;
  double t;
};

dVec3 B(double time, dVec3 pos);

dVec3 E(double time, dVec3 pos);

dVec3 acc(double time, dVec3 pos, dVec3 vel);

void rundekutta(double time, movement& path);

int main();