#include <iostream>
#pragma once
int Map3dto1d (int x, int y, int z, Utils::Vector3i nrOfCells){
	int result=x+(y*nrOfCells[0])+(z*nrOfCells[0]*nrOfCells[1]);
	return result;
}	

int Map3dto1d(Utils::Vector3i arr, Utils::Vector3i nrOfCells){
	return Map3dto1d (arr[0] ,arr[1],arr[2],nrOfCells);
}

int Map3dto1d(std::array<int,3> arr, Utils::Vector3i nrOfCells){
	return Map3dto1d (arr[0] ,arr[1],arr[2],nrOfCells);
}
Utils::Vector3i Map1dto3d(int i, Utils::Vector3i nrOfCells){
	Utils::Vector3i result;
	result[2]=i/(nrOfCells[0]*nrOfCells[1]);
	i-=result[2]*(nrOfCells[0]*nrOfCells[1]);
	result[1]=i/nrOfCells[0];
	i-=result[1]*nrOfCells[0];
	result[0]=i;
	return result;

}

std::pair<int,Utils::Vector3i> findShifted(int xMove, int yMove, int zMove, int index, Utils::Vector3i nrOfCells){
	auto pos = Map1dto3d(index, nrOfCells);
	pos[0]+=xMove;
	pos[1]+=yMove;
	pos[2]+=zMove;
	std::pair<int,Utils::Vector3i> result;


	for (int i:{0,1,2}){
		if(pos[i]<0) {pos[i]+=nrOfCells[i];
			result.second[i]=-1;}
		else if (pos[i]>=nrOfCells[i]) {pos[i]-=nrOfCells[i];
			result.second[i]=+1;
		}
	}
	result.first=Map3dto1d(pos[0],pos[1],pos[2],nrOfCells);
	return result;

}
std::pair<int,Utils::Vector3i> findShifted(Utils::Vector3i arr, int index, Utils::Vector3i nrOfCells){
	return findShifted(arr[0],arr[1],arr[2],index,nrOfCells);
}

