#pragma once
#include "maps_between_1d_3d.hpp" 
#include <iostream>
#include "NeighborSetup.hpp"


enum class FORCES : int {
	LJ,
	EWALDRS,
	EWALDKS
};


template<typename Particle>
struct Box{

	private:
		BC boundaryCon_;
		double ljcutoff_;
		double ljcutoff2_;
		double boxSize1d_;
		double eps_;
		double sigma_;
		double coulombC_;
		std::vector<Cell<Particle>> listofCells_;
		Utils::Vector3d boxSize_;
		Utils::Vector3d cellSize_;
		Utils::Vector3i nrOfCells_;
		std::vector<FORCES> forces_;

		double dt_=0.01;
		int currenttimestep_=0;
		int seed_=12;
		double ewaldcutoff_;
		double ewaldcutoff2_;
		int kmax_;
		double alpha_;
		double temperature_=1;
		double frictionCoef_=1;


	public:
		Cell<Particle>& operator[](int idx) { return listofCells_[idx]; }
		int size() {return listofCells_.size();}
		std::vector<Cell<Particle>>& all() {return listofCells_;}
		Utils::Vector3d boxSize() {return boxSize_;}
		Utils::Vector3d cellSize() {return cellSize_;}
		Utils::Vector3i nrOfCells() {return nrOfCells_;}
		std::vector<FORCES> forces() {return forces_;}
		double ljcutoff() {return ljcutoff_;}
		double ljcutoff2() {return ljcutoff2_;}
		double eps() {return eps_;}
		double sigma() {return sigma_;}
		double coulombC() {return coulombC_;}

		double dt() {return dt_;}

		double ewaldcutoff() {return ewaldcutoff_;}
		double ewaldcutoff2() {return ewaldcutoff2_;}
		double alpha() {return alpha_;}
		int kmax() {return kmax_;}

		double temperature() {return temperature_;}
		double frictionCoef() {return frictionCoef_;}

		int currenttimestep() {return currenttimestep_;}
		int settimestep(int a){currenttimestep_=a; return currenttimestep_;}
		int increasetimestep() {currenttimestep_++; return currenttimestep_;}
		int seed() {return seed_;}

		Box(Utils::Vector3d boxSize, double ljcutoff, BC boundaryCondition, double eps, double sigma, double coulombC, double ewaldcutoff, double alpha, int kmax,std::vector<FORCES> forces) : boxSize_(boxSize),ljcutoff_(ljcutoff),ljcutoff2_(ljcutoff*ljcutoff),eps_(eps),sigma_(sigma), boundaryCon_(boundaryCondition), coulombC_(coulombC),ewaldcutoff_(ewaldcutoff),ewaldcutoff2_(ewaldcutoff*ewaldcutoff),alpha_(alpha),kmax_(kmax),forces_(forces) {
			// use a cell size smaller than the ljcutoff to reduce spurious pairs
			double cutoff=std::max(ljcutoff_,ewaldcutoff_);
			for (int i:{0,1,2}){
				nrOfCells_[i]=(int)(boxSize_[i]/cutoff*32./20.);//32./20.
				cellSize_[i]=boxSize_[i]/(double)nrOfCells_[i];
			}
			int nrOfCells=nrOfCells_[0]*nrOfCells_[1]*nrOfCells_[2];
			std::cout <<"nrOfCells="<<nrOfCells<<" ("<<nrOfCells_[0]<<", "<<nrOfCells_[1]<<", "<<nrOfCells_[2]<<")\n";
			listofCells_=std::vector<Cell<Particle>>(nrOfCells);
			NeighborSetup(listofCells_,nrOfCells_,cellSize_,ljcutoff2_,boundaryCon_);
		}
};
