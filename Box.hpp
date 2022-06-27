#pragma once
#include "maps_between_1d_3d.hpp" 
#include <iostream>
#include "NeighborSetup.hpp"




template<typename Particle>
struct Box{

	private:
		BC boundaryCon_;
		double cutoff_;
		double cutoff2_;
		double boxSize1d_;
		double eps_;
		double sigma_;
		std::vector<Cell<Particle>> listofCells_;
		Utils::Vector3d boxSize_;
		Utils::Vector3d cellSize_;
		Utils::Vector3i nrOfCells_;
	public:
		Cell<Particle>& operator[](int idx) { return listofCells_[idx]; }
		int size() {return listofCells_.size();}
		std::vector<Cell<Particle>>& all() {return listofCells_;}
		Utils::Vector3d boxSize() {return boxSize_;}
		Utils::Vector3d cellSize() {return cellSize_;}
		Utils::Vector3i nrOfCells() {return nrOfCells_;}
		double cutoff() {return cutoff_;}
		double cutoff2() {return cutoff2_;}
		double eps() {return eps_;}
		double sigma() {return sigma_;}

		Box(Utils::Vector3d boxSize, double cutoff, BC boundaryCondition, double eps = 0, double sigma = 0) : boxSize_(boxSize),cutoff_(cutoff),cutoff2_(cutoff*cutoff),eps_(eps),sigma_(sigma), boundaryCon_(boundaryCondition) {
			// use a cell size smaller than the cutoff to reduce spurious pairs
			for (int i:{0,1,2}){
				nrOfCells_[i]=(int)(boxSize_[i]/cutoff_*20./20.);//32./20.
				cellSize_[i]=boxSize_[i]/(double)nrOfCells_[i];
			}
			int nrOfCells=nrOfCells_[0]*nrOfCells_[1]*nrOfCells_[2];
			std::cout <<"nrOfCells="<<nrOfCells<<" ("<<nrOfCells_[0]<<", "<<nrOfCells_[1]<<", "<<nrOfCells_[2]<<")\n";
			listofCells_=std::vector<Cell<Particle>>(nrOfCells);
			NeighborSetup(listofCells_,nrOfCells_,cellSize_,cutoff2_,boundaryCon_);
		}
};
