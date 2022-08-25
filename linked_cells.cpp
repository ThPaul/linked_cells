
#include "Forces.hpp"
#include "Boxfunctions.hpp"
#include "MinimalFlatParticle.hpp"
#include "Box.hpp"
#include "InputData.hpp"
#include <vector>
#include <iostream>
#include <chrono>
#include <hpx/hpx_main.hpp>
#include "velocity_verlet_langevin.hpp"


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



auto kernel(Box<MinimalFlatParticle<0>>& box) {
	auto  begin= std::chrono::steady_clock::now();
	calc_forces<HPX_PROTOCOL::NONE>(box);
	//calc_ewald_rs_force(box);
	Utils::Vector3i a={30,30,30};
	//kSpaceForces(a,box);
	auto end = std::chrono::steady_clock::now();
	auto res = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
	std::cout << "time elapsed " << res << " [ms]"<<std::endl;
	return static_cast<double>(res);
}

int main(int argc, char** argv) {
	auto settings = readInput("/tikhome/tpaul/Documents/linked_cells/build/settings");
	auto contains=[forces=settings.forces](FORCES f){return std::find(forces.begin(),forces.end(),f)!=forces.end();};
//	if(settings.ewaldacc!=-1&&(contains(FORCES::EWALDKS)||contains(FORCES::EWALDRS))) ewaldTuning(settings);
	Box<MinimalFlatParticle<0>>box(settings.boxSize,settings.LJcutoff,settings.boundaryCondition,settings.eps,settings.sigma,settings.coulombC,settings.EwaldCutoff,settings.alpha,settings.kmax,settings.forces);
	fill_Cell(box,settings.posFile,settings.chargeFile);


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
