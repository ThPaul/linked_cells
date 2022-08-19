#pragma once
#include "Boxfunctions.hpp"
#include "EwaldKspace.hpp"



template <HPX_PROTOCOL hpx_protocol, typename Particle, typename BinaryOp>
void calc_short_range_forces(Box<Particle>& box, BinaryOp op){

	auto inner=[op](auto& par1, auto& par2){
		auto dist=par2.pos()-par1.pos();
		auto norm2 = dist.norm2();
			op(par1,par2,norm2,dist);

		
	};

	auto neighbor=[op](auto& par1, auto& par2, auto& offset){
		auto dist=par2.pos()-par1.pos()+offset;
		auto norm2 = dist.norm2();
			op(par1,par2,norm2,dist);

		
	};


	op_on_pairs_within_cutoff<hpx_protocol>(box,inner,neighbor);



}


template <HPX_PROTOCOL hpx_protocol, typename Particle>
void calc_forces(Box<Particle>& box){

	auto ewald_rs_force=[alpha=box.alpha(),coulombC=box.coulombC(),cutoff2=box.ewaldcutoff2()](auto& par1,auto& par2, double dist2, auto dist){
		if(dist2<=cutoff2){
			double distnorm = sqrt(dist2);
			const double pi = M_PI;
			double part1 =2.0*alpha/sqrt(pi)*exp(-alpha*alpha*dist2);
			double part2 =erfc(alpha*distnorm)/distnorm;
			double charge = par1.charge()*par2.charge();
			double force_norm =- coulombC*charge*(part1+part2)/dist2; // =force/distance
			Utils::Vector3d forcevector = dist*force_norm;
			par1.force()+= forcevector;
			par2.force()-= forcevector;
		}
	};


	auto lj_force=[eps=box.eps(),sigma=box.sigma(),cutoff2=box.ljcutoff2()](auto& par1,auto& par2, double dist2, auto dist){
		if(dist2<=cutoff2){
			auto const sig_r2 = sigma * sigma / dist2;
			auto const sig_r6 = sig_r2 * sig_r2 * sig_r2;
			auto const sig_r12= sig_r6 * sig_r6;
			double force_norm = - 24*eps/dist2*(2*sig_r12-sig_r6); // =force/distance
			Utils::Vector3d forcevector = dist*force_norm;
			par1.force()+= forcevector;
			par2.force()-= forcevector;
		}
	};

	auto lj_ewald_force=[ewald_rs_force,lj_force](auto& par1,auto& par2, double dist2, auto dist){
		//could save some dist <= cutoff) coparisons
		lj_force(par1,par2,dist2,dist);
		ewald_rs_force(par1,par2,dist2,dist);
	};
	auto contains=[forces=box.forces()](FORCES f){return std::find(forces.begin(),forces.end(),f)!=forces.end();};
	if (contains(FORCES::LJ)&&contains(FORCES::EWALDRS))calc_short_range_forces<hpx_protocol>(box,lj_ewald_force);
	else if (contains(FORCES::LJ))calc_short_range_forces<hpx_protocol>(box,lj_force);
	else if (contains(FORCES::EWALDRS)) calc_short_range_forces<hpx_protocol>(box,ewald_rs_force);
	if (contains(FORCES::EWALDKS)) kSpaceForces<hpx_protocol>(box);

}

