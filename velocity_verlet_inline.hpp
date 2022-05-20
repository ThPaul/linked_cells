/*
 * Copyright (C) 2010-2019 The ESPResSo project
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef INTEGRATORS_VELOCITY_VERLET_HPP
#define INTEGRATORS_VELOCITY_VERLET_HPP

#include "config.hpp"


/** Propagate the velocities and positions. Integration steps before force
 *  calculation of the Velocity Verlet integrator: <br> \f[ v(t+0.5 \Delta t) =
 *  v(t) + 0.5 \Delta t f(t)/m \f] <br> \f[ p(t+\Delta t) = p(t) + \Delta t
 *  v(t+0.5 \Delta t) \f]
 */
template <typename ParticleIterable>
void velocity_verlet_propagate_vel_pos(ParticleIterable &particles,
                                              double time_step) {

  for (auto &p : particles) {

    for (int j = 0; j < 3; j++) {
        /* Propagate velocities: v(t+0.5*dt) = v(t) + 0.5 * dt * a(t) */
        p.v()[j] += 0.5 * time_step * p.force()[j] / p.mass();

        /* Propagate positions (only NVT): p(t + dt)   = p(t) + dt *
         * v(t+0.5*dt) */
        p.pos()[j] += time_step * p.v()[j];
      }
    }

}

/** Final integration step of the Velocity Verlet integrator
 *  \f[ v(t+\Delta t) = v(t+0.5 \Delta t) + 0.5 \Delta t f(t+\Delta t)/m \f]
 */
template <typename ParticleIterable> 
void velocity_verlet_propagate_vel_final(ParticleIterable &particles,
                                                double time_step) {

  for (auto &p : particles) {
    // Virtual sites are not propagated during integration
    for (int j = 0; j < 3; j++) {
        /* Propagate velocity: v(t+dt) = v(t+0.5*dt) + 0.5*dt * a(t+dt) */
        p.v()[j] += 0.5 * time_step * p.force()[j] / p.mass();
      }
    }
}

template <typename ParticleIterable>
inline void velocity_verlet_step_1(ParticleIterable &particles,
                                   double time_step) {
  velocity_verlet_propagate_vel_pos(particles, time_step);
}

template <typename ParticleIterable>
inline void velocity_verlet_step_2(ParticleIterable &particles,
                                   double time_step) {
  velocity_verlet_propagate_vel_final(particles, time_step);
}

#endif
