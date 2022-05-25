#include "Cell.hpp"
#include "MinimalFlatParticle.hpp"
#include "utils/Vector.hpp"
#include "utils/for_each_pair.hpp"
#include "maps_between_1d_3d.hpp"
#include "Box.hpp"
#include <vector>
#include <iostream>




template  <typename Particle>
void fill_Cell(Box<Particle>& cells, std::vector<Utils::Vector3d>particleCoords){
	for (auto par : particleCoords){
		Utils::Vector3d d;
		for(int i:{0,1,2}) d[i]=par[i]/cells.cellSize()[i];
		std::array<int,3>cellcoords ={(int)d[0],(int)d[1],(int)d[2]};
		int cellindex=Map3dto1d(cellcoords,cells.nrOfCells());
		Particle particle;
		particle.pos()=par;
		cells[cellindex].particles().insert(particle);
		//std::cout<<"put particle " << particle.pos() << "    in cell "<<cellindex<<std::endl;
	}
};

template <typename Particle, typename BinaryOp>
void op_on_pairs_within_cutoff(Box<Particle>& cells,BinaryOp op){
	auto op_within_cutoff=[cutoff2=cells.cutoff2(),op](auto& par1, auto& par2){
		auto dist=par2.pos()-par1.pos();
		auto norm2 = dist.norm2();
		if(norm2<=cutoff2){
			op(par1,par2);

		}
	};
	for (auto& cell : cells.all()){
		Utils::for_each_pair(cell.particles().begin(), cell.particles().end(), op_within_cutoff);
		for  (auto& redNeighbor : cell.neighbors().red()){
			if (redNeighbor!=NULL){ //BC
				for(auto&  par1 : cell.particles()){
					for (auto& par2 : redNeighbor->particles()){
						op_within_cutoff(par1,par2);
					}
				}
			}
		}
	}
}

template <typename Particle>
void print_forces(Box<Particle> cells){

	std::cout<<"Number of Cells  " << cells.nrOfCells() << "  size of one cell: "<<cells.cellSize()<<std::endl;
	int counter=0;
	for (auto cell : cells.all()){
		for (auto particle : cell.particles()){
			auto pos=particle.pos();
			auto force=particle.force();
			std::cout<<pos<<" "<<force<<"    in Cell "<<counter<<std::endl;

		}
		counter++;
	}
}

int main(int argc, char** argv) {
	double cutoff= 1;
	Utils::Vector3d BoxSize={15,25,35};
	auto op=[](auto& par1,auto& par2){
		par1.force()[0]+=1;
		par2.force()[1]-=1;
		par1.force()[2]+=1;
		par2.force()[2]+=1;
	};



	auto run=[cutoff,BoxSize,op](std::vector<Utils::Vector3d> particles){
		Box<MinimalFlatParticle<0>>cells(cutoff,BoxSize);
		fill_Cell(cells,particles);
		op_on_pairs_within_cutoff(cells,op);
		print_forces(cells);
	};
	using Particles = std::vector<Utils::Vector3d>;
	Particles a = {{1.1,1.1,1.1},{1.1,1.2,1.2}};
	run(a);
	Particles b = {{1,0,0},{2,0,0},{3,0,0},{4,0,0}};
	run(b);
	Particles c;
	for (double i=0.1;i<5.0;i+=0.6){
		for  (double j=0.1;j<5.0;j+=0.6){
			for(double k=0.1;k<5.0;k+=0.6){
				c.push_back({i,j,k});
			}
		}
	}
	run(c);
	return 0; }
