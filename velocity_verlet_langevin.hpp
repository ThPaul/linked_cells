#include "random.hpp"

double random_uniform(int counter,int seed, int key, int step){
	auto a = Random::noise_uniform<RNGSalt::LANGEVIN>(counter,seed,key,step);
	return a[0]; //random uniform -0.5 to 0.5
}


template  <typename Particle>
void half_velocity(Particle& p,double temp,double friction,double dt,int counter, int seed,int step){
	double variance = 2*p.mass()*friction*temp/dt;
	double randomforce=random_uniform(counter,seed,p.index(),step)*sqrt(12*variance);
	p.vel()=p.vel()+0.5*dt*(p.force()/p.mass()-friction*p.vel()+randomforce/p.mass());
}

template  <typename Particle>
void new_positions(Particle& p,double dt){
	p.pos()+=dt*p.vel();
}

template  <typename Particle>
void verlet_step_1(Box<Particle>& box){
	for (auto& cell:box.all()){
		for (auto& p:cell.particles()){
			half_velocity(p,box.temperature(),box.frictionCoef(),box.dt(),box.currenttimestep(),box.seed(),1);
			new_position(p,box.dt());
		}
	}
}


template  <typename Particle>
void verlet_step_2(Box<Particle>& box){
	for (auto& cell:box.all()){
		for (auto& p:cell.particles()){
			half_velocity(p,box.temperature(),box.frictionCoef(),box.dt(),box.currenttimestep(),box.seed(),2);
		}
	}
}
