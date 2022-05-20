#pragma once

#include "soa/soa.h"
#include "utils/Vector.hpp"
enum class PPROPS: std::size_t {
   POS=0, V=1, FORCE=2, MASS=3, EXTRA_DATA };

constexpr const std::size_t prop_idx(const PPROPS p) { return std::size_t(p); };


struct ExtraData { char x[1000]; };

using ParticlePool = SoA<
    Utils::Vector3d, // pos
    Utils::Vector3d, // vel
    Utils::Vector3d, // force
    double,
    ExtraData>; // mass

