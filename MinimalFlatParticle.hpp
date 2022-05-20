#include "utils/Vector.hpp"

template <int extra_size=0> 
struct MinimalFlatParticle {
private:
  Utils::Vector3d _pos{};
  Utils::Vector3d _v{};
  Utils::Vector3d _force{};
 double _mass;
 char _waste[extra_size];
public:
  auto const& pos() const { return _pos; };
  auto const& v()  const{ return _v; };
  auto const& force() const { return _force; };
  auto const& mass() const { return _mass; };
  
  auto & pos()  { return _pos; };
  auto & v() { return _v; };
  auto & force() { return _force; };
  auto & mass() { return _mass; };
};
