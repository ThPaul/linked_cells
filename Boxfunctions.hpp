#include "Cell.hpp"
#include "MinimalFlatParticle.hpp"
#include "utils/Vector.hpp"
#include "utils/for_each_pair.hpp"
#include "maps_between_1d_3d.hpp"
#include "Box.hpp"
#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <hpx/future.hpp>



template  <typename Particle>
void fill_Cell(Box<Particle>& box, std::vector<Utils::Vector3d>particleCoords){
	for (auto par : particleCoords){
		Utils::Vector3d d;
		for(int i:{0,1,2}) d[i]=par[i]/box.cellSize()[i];
		std::array<int,3>cellcoords ={(int)d[0],(int)d[1],(int)d[2]};
		int cellindex=Map3dto1d(cellcoords,box.nrOfCells());
		Particle particle;
		particle.pos()=par;
		box[cellindex].particles().insert(particle);
	}
};

template <typename Particle>
void fill_Cell(Box<Particle>& box, std::string path_to_file){
	std::ifstream particleFile(path_to_file);
	std::vector<Utils::Vector3d>particleCoords;
	if(particleFile){
		std::string line;
		while(getline(particleFile,line)){
			std::vector<std::string> split;
			boost::split(split,line,[](char c){return c==',';});
			Utils::Vector3d coords;
			for (int i : {0,1,2}){
				coords[i]=std::stod(split[i]);

			}
			particleCoords.push_back(coords);
		}
		particleFile.close();
		fill_Cell(box,particleCoords);
	}
	else std::cout<< "FILE DOESNT EXIST" <<std::endl;


}	


template <typename Particle, typename BinaryOp>
void op_on_pairs_within_cutoff(Box<Particle>& box,BinaryOp op){
	
	auto op_within_cutoff=[cutoff2=box.cutoff2(),op](auto& par1, auto& par2){
		auto dist=par2.pos()-par1.pos();
		auto norm2 = dist.norm2();
		if(norm2<=cutoff2){
			op(par1,par2);

		}
	};
	for (auto& cell : box.all()){
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



template <typename Particle, typename BinaryOp>
void op_on_pairs_within_cutoff_hpx(Box<Particle>& box,BinaryOp op){
		auto op_within_cutoff=[cutoff2=box.cutoff2(),op](auto& par1, auto& par2){
		auto dist=par2.pos()-par1.pos();
		auto norm2 = dist.norm2();
		if(norm2<=cutoff2){
			op(par1,par2);

		}
	};
	auto op_on_cell=[op_within_cutoff](auto& cell){
	
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
	
	
	};
	std::vector<hpx::future<void>> fut;
	//fut.push_back(hpx::async(op_on_pairs_within_cutoff<Particle,BinaryOp>,std::ref(box),op));
	
	for (auto& cell : box.all()){
		fut.push_back(hpx::async(op_on_cell,std::ref(cell)));


	}
	hpx::wait_all(fut);
}


template <typename Particle>
void print_forces(Box<Particle> box){

	std::cout<<"Number of Cells  " << box.nrOfCells() << "  size of one cell: "<<box.cellSize()<<std::endl;
	for (auto cell : box.all()){
		for (auto particle : cell.particles()){
			auto pos=particle.pos();
			auto force=particle.force();
			std::cout<<pos<<" "<<force<<std::endl;

		}
	}
}


template <typename Particle>
void print_forces_sorted(Box<Particle> box){

	std::cout<<"Number of Cells  " << box.nrOfCells() << "  size of one cell: "<<box.cellSize()<<std::endl;
	std::vector<Particle> allPart;
	for (auto cell : box.all()){
		for (auto particle : cell.particles()){
			allPart.push_back(particle);

		}
	}
	sort(allPart.begin(),allPart.end(),[](auto a, auto b){return a.pos().norm2() < b.pos().norm2();});
	for (auto p : allPart){
		std::cout<<p.force()<<std::endl;}
}

template <typename Particle>
void print_forces_sorted(Box<Particle> box, std::string file){
	std::vector<Particle> allPart;
	for (auto cell : box.all()){
		for (auto particle : cell.particles()){
			allPart.push_back(particle);

		}
	}
	sort(allPart.begin(),allPart.end(),[](auto a, auto b){return a.pos().norm2() < b.pos().norm2();});
	std::ofstream myfile(file);
	for (auto p : allPart){
		auto a =[](auto f){return boost::lexical_cast<std::string>(f);};
		myfile<<a(p.force()[0])<<","<<a(p.force()[1])<<","<<a(p.force()[2])<<std::endl;

	}

	myfile.close();
}
