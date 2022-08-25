#pragma once
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


enum class HPX_PROTOCOL : int {
	NONE,
	ASYNC
};





template  <typename Particle>
void fill_Cell(Box<Particle>& box, std::vector<Utils::Vector3d>particleCoords, std::vector<double> particleCharges){
	for (int i=0; i< particleCoords.size();i++){
		auto par = particleCoords[i];
		Utils::Vector3d d;
		for(int i:{0,1,2}) d[i]=par[i]/box.cellSize()[i];
		std::array<int,3>cellcoords ={(int)d[0],(int)d[1],(int)d[2]};
		int cellindex=Map3dto1d(cellcoords,box.nrOfCells());
		Particle particle;
		particle.pos()=par;
		particle.charge()=particleCharges[i];
		particle.index()=i;
		box[cellindex].particles().insert(particle);
	}
};

template <typename Particle>
void fill_Cell(Box<Particle>& box, std::string file_pos, std::string file_charge=""){
	std::ifstream particleFile(file_pos);
	std::vector<Utils::Vector3d>particleCoords;
	std::vector<double>particleCharges;
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



	}
	else{ std::cout<< "POS FILE DOESNT EXIST" <<std::endl;}
	
	if (file_charge!=""){
		std::ifstream particleFile(file_charge);
        	if(particleFile){
                	std::string line;
                	while(getline(particleFile,line)){
                        	double charge=std::stod(line);
                        	particleCharges.push_back(charge);
               		 }
                	particleFile.close();
		}
        	else std::cout<< "CHARGE FILE DOESNT EXIST" <<std::endl;
	}
	else  particleCharges.resize(particleCoords.size(),0.0);


	fill_Cell(box,particleCoords,particleCharges);

}	

template <HPX_PROTOCOL hpx_protocol, typename Particle, typename BinaryOp1, typename BinaryOp2>
void op_on_pairs_within_cutoff(Box<Particle>& box,BinaryOp1 inner, BinaryOp2 neighbor){


	auto op_on_cell=[](auto& cell,auto inner, auto neighbor){
		{
			std::scoped_lock lock(cell.cellMutex);
			Utils::for_each_pair(cell.particles().begin(), cell.particles().end(), inner);
		}
		for  (auto& redNeighbor : cell.neighbors().red()){
			auto offset=redNeighbor.offset();
			std::scoped_lock lock(cell.cellMutex,redNeighbor.cellRef()->cellMutex);
			for(auto&  par1 : cell.particles()){
				for (auto& par2 : redNeighbor.cellRef()->particles()){
					neighbor(par1,par2,offset);

					}
				}
			}
		};

	auto final_op=[op_on_cell,inner,neighbor](auto& cell){
	op_on_cell(cell,inner,neighbor);};

	if constexpr (hpx_protocol == HPX_PROTOCOL::ASYNC) {
		std::vector<hpx::future<void>> fut;
		for (auto& cell : box.all()){
			fut.push_back(hpx::async(final_op,std::ref(cell)));
		}
		hpx::wait_all(fut);
	} else if constexpr (hpx_protocol == HPX_PROTOCOL::NONE) {
		for (auto& cell : box.all()){
			final_op(cell);
		}
	}
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
void print_forces_sorted(Box<Particle>& box){

	std::cout<<"Number of Cells  " << box.nrOfCells() << "  size of one cell: "<<box.cellSize()<<std::endl;
	std::vector<Particle> allPart;
	for (auto& cell : box.all()){
		for (auto& particle : cell.particles()){
			allPart.push_back(particle);

		}
	}
	sort(allPart.begin(),allPart.end(),[](auto a, auto b){return a.pos().norm2() < b.pos().norm2();});
	for (auto p : allPart){
		std::cout<<p.force()<<std::endl;}
}

template <typename Particle>
void print_forces_sorted(Box<Particle>& box, std::string file){
	std::vector<Particle> allPart;
	for (auto& cell : box.all()){
		for (auto& particle : cell.particles()){
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
