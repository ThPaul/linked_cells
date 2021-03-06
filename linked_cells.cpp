#include "Forces.hpp"
#include "Boxfunctions.hpp"
#include "MinimalFlatParticle.hpp"
#include "Box.hpp"
#include <vector>
#include <iostream>
#include <chrono>
#include <hpx/hpx_main.hpp>


#include <vector>
#include <algorithm>
#include <numeric>
#include <utility>

namespace detail {
std::pair<double, double> statistics(std::vector<double> &vec) {
	if (vec.empty()) {
		return {0.0, 0.0};
	}
	auto const n = static_cast<double>(vec.size());
	auto const mean = std::accumulate(vec.begin(), vec.end(), 0.0) / n;
	auto var_func = [mean, n](double acc, const double val) {
		auto const diff = val - mean;
		return acc + diff*diff/(n-1.);
	};
	auto const var = std::accumulate(vec.begin(), vec.end(), 0.0, var_func);
	return {mean, var};
}
} // namespace detail

template <typename Particle>
void calc_ewald_rs_force(Box<Particle>& box){
	auto ewald_rs_force=[coulombC=box.coulombC()](auto& par1,auto& par2, double dist2, auto dist){
		double alpha=0.666;
		double distnorm = sqrt(dist2);
		const double pi = M_PI;
		double part1 =2.0*alpha/sqrt(pi)*exp(-alpha*alpha*dist2);
		double part2 =erfc(alpha*distnorm)/distnorm;
		double charge = par1.charge()*par2.charge();
		double force_norm = coulombC*charge*(part1+part2)/dist2; // =force/distance
		Utils::Vector3d forcevector = dist*force_norm;
		par1.force()+= forcevector;
		par2.force()-= forcevector;
	};

short_range_forces<HPX_PROTOCOL::ASYNC>(box,ewald_rs_force);
}

template <typename Particle>
void calc_lj_force(Box<Particle>& box){
	auto lj_force=[eps=box.eps(),sigma=box.sigma()](auto& par1,auto& par2, double dist2, auto dist){
		auto const sig_r2 = sigma * sigma / dist2;
		auto const sig_r6 = sig_r2 * sig_r2 * sig_r2;
		auto const sig_r12= sig_r6 * sig_r6;
		double force_norm = 24*eps/dist2*(2*sig_r12-sig_r6); // =force/distance
		Utils::Vector3d forcevector = -dist*force_norm;
		par1.force()+= forcevector;
		par2.force()-= forcevector;
	};

short_range_forces<HPX_PROTOCOL::ASYNC>(box,lj_force);
}

auto kernel(Box<MinimalFlatParticle<0>>& box) {
	auto  begin= std::chrono::steady_clock::now();
	//calc_lj_force(box);
	//calc_ewald_rs_force(box);

	Utils::Vector3i a={30,30,30};
	kSpaceForces(a,box);
	auto end = std::chrono::steady_clock::now();
	auto res = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
	std::cout << "time elapsed " << res << " [ms]"<<std::endl;
	return static_cast<double>(res);
}

int main(int argc, char** argv) {
	double cutoff= 2.5;
	Utils::Vector3d boxSize={50,50,50};
	double eps=1.0;
	double sigma=1.0;
	BC bc= BC::PERIODIC;
	Box<MinimalFlatParticle<0>>box(boxSize,cutoff,bc,eps,sigma);
	fill_Cell(box,"/tikhome/tpaul/Documents/linked_cells/build/partPos","/tikhome/tpaul/Documents/linked_cells/build/charges");


	std::vector<double> timings = {};
	for (int i = 0; i < 1; ++i) {
	    for (auto& c : box.all()){
		    for (auto& p : c.particles()){
			    p.force()={0,0,0};
			    }}
	    timings.push_back(kernel(box));
	}
	auto const stats = detail::statistics(timings);
	auto const timings_avg = stats.first;
	auto const timings_sem = stats.second / std::sqrt(static_cast<double>(timings.size() - 1));
	auto const timings_c95 = static_cast<int>(1.96 * timings_sem);
	std::cout << "stats: " << timings_avg << " +/- " << timings_c95 << " [ms]"<<std::endl;
	print_forces_sorted(box,"forces_lc");

	return 0; }
