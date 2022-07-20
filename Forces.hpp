#include "Boxfunctions.hpp"
#include "EwaldKspace.hpp"

template <HPX_PROTOCOL hpx_protocol, typename Particle, typename BinaryOp>
void short_range_forces(Box<Particle>& box, BinaryOp op){


	auto inner=[op,cutoff2=box.cutoff2()](auto& par1, auto& par2){
		auto dist=par2.pos()-par1.pos();
		auto norm2 = dist.norm2();
		if(norm2<=cutoff2){
			op(par1,par2,norm2,dist);

		}
	};

	auto neighbor=[op,cutoff2=box.cutoff2()](auto& par1, auto& par2, auto& offset){
		auto dist=par2.pos()-par1.pos()+offset;
		auto norm2 = dist.norm2();
		if(norm2<=cutoff2){
			op(par1,par2,norm2,dist);

		}
	};


op_on_pairs_within_cutoff<hpx_protocol>(box,inner,neighbor);



}

