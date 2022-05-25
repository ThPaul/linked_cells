#include <atomic>
#include <vector>
#include <iostream>

#define ATOMIC 1 // 0 for int, 1 for std::atomic<int>
#if ATOMIC==1

struct Particle{
	private: std::atomic<int> force_;
	public:
		 std::atomic<int>& force() {return force_;}
};
#endif
#if ATOMIC==0

struct Particle{
	private: int force_;
	public:
		 int& force() {return force_;}
};
#endif

struct Cell{
	private:
		std::vector<Particle> listOfParticles_;
	public:
		Cell(int size):listOfParticles_(size){}
		std::vector<Particle>& particles(){return listOfParticles_;}

};

void works(){
	Particle a;
	a.force()=0;
	a.force()+=1;
	Cell myCell(5);
	myCell.particles()[2].force()=5;
	std::cout<<a.force()<<" "<<myCell.particles()[2].force()<<std::endl;
}
void doesntwork(){
	Cell myCell(3);
	Particle myParticle;
	myCell.particles().push_back(myParticle); //adding Particle
	myCell.particles().erase(myCell.particles().begin()+1); //deleting 2nd Particle
}

int main(int argc, char** argv) {
	works();
	doesntwork();
	return 0;
}
