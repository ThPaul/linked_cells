#pragma once
int Map3dto1d (int x, int y, int z, Utils::Vector3i nrOfCells){
	int result=x+(y*nrOfCells[0])+(z*nrOfCells[0]*nrOfCells[1]);
	return result;
}	

int Map3dto1d(std::array<int,3> arr, Utils::Vector3i nrOfCells){
	return Map3dto1d (arr[0] ,arr[1],arr[2],nrOfCells);
}
std::array<int,3> Map1dto3d(int i, Utils::Vector3i nrOfCells){
	std::array<int,3> result;
	result[2]=i/(nrOfCells[0]*nrOfCells[1]);
	i-=result[2]*(nrOfCells[0]*nrOfCells[1]);
	result[1]=i/nrOfCells[0];
	i-=result[1]*nrOfCells[0];
	result[0]=i;
	return result;

}
// find index of Cell with given relative Postion, returns -1 if outside of Box
int findShifted(int xMove, int yMove, int zMove, int index, Utils::Vector3i nrOfCells){
	auto pos = Map1dto3d(index, nrOfCells);
	pos[0]+=xMove;
	pos[1]+=yMove;
	pos[2]+=zMove;
	bool inRange=true;
	for (int i:{0,1,2}){
		inRange*=(pos[i]>=0&&pos[i]<nrOfCells[i]);
	}
	if(inRange)
		return Map3dto1d(pos[0],pos[1],pos[2],nrOfCells);
	else{
		return -1;}
}
int findShifted(std::array<int,3> arr, int index, Utils::Vector3i nrOfCells){
	return findShifted(arr[0],arr[1],arr[2],index,nrOfCells);
}

