#include "random.hpp"

Utils::Vector3d random_uniform(int counter,int seed, int key, int step){
	auto a = Random::noise_uniform<RNGSalt::LANGEVIN>(counter,seed,key,step);
	Utils::Vector3d result={a[0],a[1],a[2]};
	return result; //random uniform -0.5 to 0.5
}

Utils::Vector3d maxf;
Utils::Vector3d maxrf;

template  <typename Particle>
void new_velocity(Particle& p,double temp,double friction,double dt,int counter, int seed,int step){
	double variance = 2*p.mass()*friction*temp/dt;
	auto  randomforce=random_uniform(counter,seed,p.index(),step)*sqrt(12*variance);
	p.v()=p.v()+dt*(p.force()/p.mass()-friction*p.v()+randomforce/p.mass());
	for (auto a : {0,1,2}) {
	maxf[a]=std::max(maxf[a],std::abs(p.force()[a]));
	maxrf[a]=std::max(maxrf[a],std::abs(randomforce[a]));

	}
}

template  <typename Particle>
void new_position(Particle& p,double dt){
	p.pos()+=dt*p.v(); //check for outside the box in move_to_new_cell()
	}



template  <typename Particle>
void verlet_step_1(Box<Particle>& box){
	maxf={0,0,0};
	maxrf={0,0,0};
	for (auto& cell:box.all()){
		for (auto& p:cell.particles()){
			new_velocity(p,box.temperature(),box.frictionCoef(),box.dt()/2,box.currenttimestep(),box.seed(),1);
			new_position(p,box.dt());
		}
	}

	std::cout<<"f: "<<maxf<<"     rf: "<<maxrf<<std::endl;
	move_to_new_cells(box);

}


template  <typename Particle>
void verlet_step_2(Box<Particle>& box){
	for (auto& cell:box.all()){
		for (auto& p:cell.particles()){
			new_velocity(p,box.temperature(),box.frictionCoef(),box.dt()/2,box.currenttimestep(),box.seed(),2);
		}
	}
}


template  <typename Particle>
void move_to_new_cells(Box<Particle>& box){
	auto nrCells=box.nrOfCells();
	for (int n=0; n< box.all().size(); n++){
		auto it =box.all()[n].particles().begin();
		while (it!= box.all()[n].particles().end()){
			Utils::Vector3d d;
			for(int i:{0,1,2}){
				if (it->pos()[i]<0) it->pos()[i]+=box.boxSize()[i];
				else if (it->pos()[i]>=box.boxSize()[i]) it->pos()[i]-=box.boxSize()[i];
			       	d[i]=it->pos()[i]/box.cellSize()[i];
			}
			std::array<int,3>cellcoords ={(int)d[0],(int)d[1],(int)d[2]};
			int cellindex=Map3dto1d(cellcoords,nrCells);
			if(cellindex!=n){
				//std::cout<<"old "<< box.all()[n].particles().size()<<" "<<box.all()[cellindex].particles().size()<<std::endl;
				box.all()[cellindex].particles().insert(*it);

				it=box.all()[n].particles().erase(it);
				//std::cout<<"moved "<<it->pos()<<"from "<<Map1dto3d(n,nrCells)<< " to " <<Map1dto3d(cellindex,nrCells)<<std::endl;

				//std::cout<<"new "<< box.all()[n].particles().size()<<" "<<box.all()[cellindex].particles().size()<<std::endl;
			}
			else it++;
		}
	}
}
