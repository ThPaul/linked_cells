#include "Cell.hpp"


enum class BC : int {
	NONE,
	PERIODIC
};



std::vector<Utils::Vector3i> calcrelPosRed(Utils::Vector3d cellSize_,double cutoff2){
	auto distance2=[cellSize_](Utils::Vector3i shift){
		for (auto& a : shift){
			if(a!=0) a=abs(a)-1; 
		}
		double distance2=0;
		for (auto i : {0,1,2}) distance2+=cellSize_[i]*cellSize_[i]*shift[i]*shift[i];
		return distance2;	
	};
	Utils::Vector3d maxValue;
	for (auto i : {0,1,2}){
		Utils::Vector3i shift;
		shift.fill(0);


		while(distance2(shift)<cutoff2)
		{

			maxValue[i]=shift[i];
			shift[i]+=1;
		}

	}

	std::vector<Utils::Vector3i> relPosQuadrant;
	for (int x=maxValue[0];x>=0;x--){
		for (int y=maxValue[1];y>=0;y--){
			for (int z=maxValue[2];z>=0;z--){
				if (distance2({x,y,z})<cutoff2) relPosQuadrant.push_back({x,y,z});
			}

		}

	}
	relPosQuadrant.pop_back();
	// x>0 || y>0 && x=0 || z>0 &&x,y=0
	std::vector<Utils::Vector3i> relPosRed;
	for (auto pos : relPosQuadrant){
		relPosRed.push_back(pos);
		if(pos[0]!=0){
			if (pos[1]!=0) relPosRed.push_back({pos[0],-pos[1],pos[2]});
			if (pos[2]!=0) relPosRed.push_back({pos[0],pos[1],-pos[2]});
			if (pos[1]!=0&&pos[2]!=0) relPosRed.push_back({pos[0],-pos[1],-pos[2]});
		}
		else if (pos[1]!=0){
			if (pos[2]!=0) relPosRed.push_back({pos[0],pos[1],-pos[2]});
		}



	}
	return relPosRed;

}

template<typename Particle>
void NeighborSetup(std::vector<Cell<Particle>>& listofCells_, Utils::Vector3i nrOfCells_, Utils::Vector3d cellSize_ , double cutoff2_, BC boundaryCond_
		){
	std::vector<Utils::Vector3i>relPosRed=calcrelPosRed(cellSize_,cutoff2_);
	std::vector<Utils::Vector3i>relPosBlack;
	for (auto pos : relPosRed) relPosBlack.push_back(pos*-1);





	for (int i=0; i<listofCells_.size(); ++i){
		std::vector<Neighbor<Cell<Particle>*>> red;
		std::vector<Neighbor<Cell<Particle>*>> black;
		Utils::Vector3i nullV={0,0,0};
		for (auto relPos : relPosRed){
			auto NeighborPos = findShifted(relPos,i,nrOfCells_);
			bool isImage=(NeighborPos.second!=nullV);
			Utils::Vector3d offset;
			for (auto i : {0,1,2}) offset[i]=NeighborPos.second[i]*cellSize_[i]*nrOfCells_[i];
			if(boundaryCond_== BC::NONE){
				auto a = NeighborPos.first;
				if(!isImage) red.push_back(Neighbor(&listofCells_[NeighborPos.first],offset,isImage));}
			if(boundaryCond_==BC::PERIODIC){
				red.push_back(Neighbor(&listofCells_[NeighborPos.first],offset,isImage));
			}
		} 

		for (auto relPos : relPosBlack){

			auto NeighborPos = findShifted(relPos,i,nrOfCells_);
			Utils::Vector3d offset;
			bool isImage=(NeighborPos.second==nullV);
			for (auto i : {0,1,2}) offset[i]=NeighborPos.second[i]*cellSize_[i];

			if(boundaryCond_== BC::NONE){
				if(!isImage)black.push_back(Neighbor(&listofCells_[NeighborPos.first],offset,isImage));}
			if(boundaryCond_==BC::PERIODIC){
				black.push_back(Neighbor(&listofCells_[NeighborPos.first],offset,isImage));
			}
		} 
		//? 
		Utils::Span<Neighbor<Cell<Particle>*>> r(red);
		Utils::Span<Neighbor<Cell<Particle>*>> b(black);
		listofCells_[i].neighbors()=Neighbors<Cell<Particle>*>(r,b);

	}
}
