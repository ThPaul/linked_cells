#pragma once

#include "ParticlePool.hpp"
#include "utils/Vector.hpp"

ParticlePool global_pool;
class SoABackedParticle {
   private: 
     ParticlePool& m_pool;
     static constexpr const Utils::Vector3d vec_zero{};
     int m_idx;
   public:
     SoABackedParticle(ParticlePool& pool=global_pool) : m_pool(pool) {
       m_pool.push_back(vec_zero, vec_zero, vec_zero, 0,ExtraData{});
       m_idx = m_pool.size()-1;
     };
     Utils::Vector3d& pos() { return m_pool.get<prop_idx(PPROPS::POS)>(m_idx); };
     Utils::Vector3d& v() { return m_pool.get<prop_idx(PPROPS::V)>(m_idx); };
     Utils::Vector3d& force() { return m_pool.get<prop_idx(PPROPS::FORCE)>(m_idx); };
     double& mass() { return m_pool.get<prop_idx(PPROPS::MASS)>(m_idx); };
     
     const Utils::Vector3d& pos() const { return m_pool.get<prop_idx(PPROPS::POS)>(m_idx); };
     const Utils::Vector3d& v() const { return m_pool.get<prop_idx(PPROPS::V)>(m_idx); };
     const Utils::Vector3d& force() const { return m_pool.get<prop_idx(PPROPS::FORCE)>(m_idx); };
     const double& mass() const { return m_pool.get<prop_idx(PPROPS::MASS)>(m_idx); };
  };
