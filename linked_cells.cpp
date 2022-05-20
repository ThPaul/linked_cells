#include "Cell.hpp"
#include "MinimalFlatParticle.hpp"
#include "utils/Vector.hpp"
#include "utils/for_each_pair.hpp"
#include <vector>
#include <iostream>


int Map3dto1d (int x, int y, int z, int NrOfCells1d){
	int result=x+(y*NrOfCells1d)+(z*NrOfCells1d*NrOfCells1d);
	return result;
}	

int Map3dto1d(std::array<int,3> arr, int NrOfCells1d){
	return Map3dto1d (arr[0] ,arr[1],arr[2],NrOfCells1d);
}
std::array<int,3> Map1dto3d(int i, int NrOfCells1d){
	std::array<int,3> result;
	result[2]=i/(NrOfCells1d*NrOfCells1d);
	i-=result[2]*(NrOfCells1d*NrOfCells1d);
	result[1]=i/NrOfCells1d;
	i-=result[1]*NrOfCells1d;
	result[0]=i;
	return result;

}
// find index of Cell with given relative Postion, returns -1 if outside of Box
int findShifted(int xMove, int yMove, int zMove, int index, int NrOfCells1d){
	auto pos = Map1dto3d(index, NrOfCells1d);
	pos[0]+=xMove;
	pos[1]+=yMove;
	pos[2]+=zMove;
	bool inRange=true;
	for (int i:{0,1,2}){
		inRange*=(pos[i]>=0&&pos[i]<NrOfCells1d);
	}
	if(inRange)
		return Map3dto1d(pos[0],pos[1],pos[2],NrOfCells1d);
	else{
		return -1;}
}
int findShifted(std::array<int,3> arr, int index, int NrOfCells1d){
	return findShifted(arr[0],arr[1],arr[2],index,NrOfCells1d);
}


template  <typename Particle>
std::vector<Cell<Particle>> create_Cells(double cutoff, double boxSize){
	int nrOfCells1d = boxSize/cutoff;
	if (boxSize-nrOfCells1d*cutoff>0.000001) nrOfCells1d+=1;
	double cellSize=boxSize/nrOfCells1d;
	std::cout<<"Number of Cells in 1D " << nrOfCells1d << "  size of one cell: "<<cellSize<<std::endl;
	int nrOfCells=nrOfCells1d*nrOfCells1d*nrOfCells1d;
	std::vector<Cell<Particle>> allCells(nrOfCells);
	std::vector<std::array<int,3>>relPosRed;
	std::vector<std::array<int,3>>relPosBlack;
	std::array<int,3> temp;
	int count=0;
	for (int x: {-1,0,1}){
		for (int y :{-1,0,1}){
			for (int z: {-1,0,1}) {
				count++;
				temp={x,y,z};
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

	for (int i=0; i<allCells.size(); ++i){
		std::vector<Cell<Particle>*> red;
		std::vector<Cell<Particle>*> black;
		for (auto relPos : relPosRed){
			int index = findShifted(relPos,i,nrOfCells1d);
			if(index!=-1){
				red.push_back(&allCells[index]);
			}
			else{
				red.push_back(NULL); //BC
			}
		}

		for (auto relPos : relPosBlack){
			int index = findShifted(relPos,i,nrOfCells1d);
			if(index!=-1){
				black.push_back(&allCells[index]);
			}
			else{
				black.push_back(NULL); //BC
			}
		}
		//? 
		Utils::Span<Cell<Particle>*> r(red);
		Utils::Span<Cell<Particle>*> b(black);
		allCells[i].neighbors()=Neighbors<Cell<Particle>*>(r,b);
		/*auto p=Map1dto3d(i,nrOfCells1d);
		  int counter=0;
		  for (auto a : allCells[i].neighbors().all()){
		  if (a==NULL) counter++;}
		  std::cout << p[0] << " "<<p[1]<<" "<<p[2]<<"    "<<counter<< std::endl;
		  */	
	}

	return allCells;


}
template  <typename Particle>
void fill_Cell(std::vector<Cell<Particle>>& cells, int boxSize, std::vector<Utils::Vector3d>particleCoords){
	int numberOfCells1d=std::cbrt(cells.size());
	double cellSize=boxSize/numberOfCells1d;
	for (auto par : particleCoords){
		auto d=par/cellSize;
		std::array<int,3>cellcoords ={(int)d[0],(int)d[1],(int)d[2]};
		int cellindex=Map3dto1d(cellcoords,boxSize);
		Particle particle;
		particle.pos()=par;
		cells[cellindex].particles().insert(particle);
	}
};

template <typename Particle>
void calc_force(std::vector<Cell<Particle>>& cells, int boxSize,double sig, double eps, double cutoff){
	auto force=[sig,eps,cutoff](auto& par1, auto& par2){
		auto dist=par2.pos()-par1.pos();
		auto norm2 = dist.norm2();
		if(norm2<=cutoff*cutoff){
			par1.force()[0]+=1;
			par2.force()[1]-=1;
			par1.force()[2]+=1;
			par2.force()[2]+=1;
		}
	};
	for (auto& cell : cells){
		Utils::for_each_pair(cell.particles().begin(), cell.particles().end(), force);
		for  (auto& redNeighbor : cell.neighbors().red()){
			if (redNeighbor!=NULL){ //BC
				for(auto&  par1 : cell.particles()){
					for (auto& par2 : redNeighbor->particles()){
						force(par1,par2);
					}
				}
			}
		}
	}
}

template <typename Particle>
void print_forces(std::vector<Cell<Particle>> cells){
	for (auto cell : cells){
		for (auto particle : cell.particles()){
			auto pos=particle.pos();
			auto force=particle.force();
			std::cout<<pos<<" "<<force<<std::endl;
		}
	}
}
void do_all(std::vector<Utils::Vector3d> a,double cutoff,double BoxSize,double sig,double eps){
	
	auto cells=create_Cells<MinimalFlatParticle<0>>(cutoff,BoxSize);
	fill_Cell(cells,BoxSize,a);
	calc_force(cells,BoxSize,sig,eps,cutoff);
	print_forces(cells);

}
int main(int argc, char** argv) {
	double cutoff= 1;
	double BoxSize=5;
	double sig=1;
	double eps=1;


	auto run=[cutoff,BoxSize,sig,eps](std::vector<Utils::Vector3d> particles){
do_all(particles,cutoff,BoxSize,sig,eps);};
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
