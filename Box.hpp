#pragma once

#include "maps_between_1d_3d.hpp"
template<typename Particle>
struct Box{

	private:
		double cutoff_;
		double cutoff2_;
		double boxSize1d_;
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

		Box(double cutoff, Utils::Vector3d boxSize) : cutoff_(cutoff),boxSize_(boxSize),cutoff2_(cutoff*cutoff) {
			for (int i:{0,1,2}){
				nrOfCells_[i]=(int)(boxSize_[i]/cutoff_);
				cellSize_[i]=boxSize_[i]/(double)nrOfCells_[i];
			}
			int nrOfCells=nrOfCells_[0]*nrOfCells_[1]*nrOfCells_[2];
			listofCells_.resize(nrOfCells);
			std::vector<std::array<int,3>>relPosRed;
			std::vector<std::array<int,3>>relPosBlack;
			int count=0;
			for (int x: {-1,0,1}){
				for (int y :{-1,0,1}){
					for (int z: {-1,0,1}) {
						count++;
						std::array<int,3> temp={x,y,z};
						if(x==-1||(y==-1&&x!=1)||(z==-1&&x!=1&&y!=1)){
							relPosRed.push_back(temp);
						}
						else{
							if(x!=0||y!=0||z!=0)	
								relPosBlack.push_back(temp);

						}
					}
				}
			}

			for (int i=0; i<this->size(); ++i){
				std::vector<Cell<Particle>*> red;
				std::vector<Cell<Particle>*> black;
				for (auto relPos : relPosRed){
					int index = findShifted(relPos,i,nrOfCells_);
					if(index!=-1){
						red.push_back(&listofCells_[index]);
					}
					else{
						red.push_back(NULL); //BC
					}
				}

				for (auto relPos : relPosBlack){
					int index = findShifted(relPos,i,nrOfCells_);
					if(index!=-1){
						black.push_back(&listofCells_[index]);
					}
					else{
						black.push_back(NULL); //BC
					}
				}
				//? 
				Utils::Span<Cell<Particle>*> r(red);
				Utils::Span<Cell<Particle>*> b(black);
				listofCells_[i].neighbors()=Neighbors<Cell<Particle>*>(r,b);
			
			}
		}
};
