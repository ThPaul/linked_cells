#include "Boxfunctions.hpp"
#include "MinimalFlatParticle.hpp"
#include "Box.hpp"
#include <vector>
#include <iostream>
#include <chrono>
#include <hpx/hpx_main.hpp>

template <typename Particle>
void calc_lj_force(Box<Particle>& box){
	auto lj_force=[eps=box.eps(),sigma=box.sigma()](auto& par1,auto& par2){
		double distance2 =(par2.pos()-par1.pos()).norm2();
		double sig_r6 = pow((sigma*sigma/distance2),3.0);
		double sig_r12= sig_r6*sig_r6;
		double force_norm = 24*eps/distance2*(2*sig_r12-sig_r6); // =force/distance
		Utils::Vector3d forcevector = -(par2.pos()-par1.pos())*force_norm;
		par1.force()+= forcevector;
		par2.force()-= forcevector;
	};

	op_on_pairs_within_cutoff_hpx(box,lj_force);
}

int main(int argc, char** argv) {
	double cutoff= 2.5;
	Utils::Vector3d boxSize={50,50,50};
	double eps=1.0;
	double sigma=1.0;
	Box<MinimalFlatParticle<0>>box(boxSize,cutoff,eps,sigma);
	fill_Cell(box,"/tikhome/tpaul/Documents/linked_cells/build/partPos");
	auto  begin= std::chrono::steady_clock::now();
	calc_lj_force(box);
	auto end = std::chrono::steady_clock::now();
	std::cout << "time elapsed " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " [ms]"<<std::endl;
	print_forces_sorted(box,"forces_lc");

	return 0; }
