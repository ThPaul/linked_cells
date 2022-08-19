struct kValue{
	Utils::Vector3d k;
	double sin=0;
	double cos=0;


};



template <typename Particle>
std::vector<kValue> kSpaceSetup(Box<Particle>& box){
	std::vector<kValue> kValues;
	auto kVector=box.kSpaceSize();
	for (int i=0; i<=kVector[0]; i++){
		for (int j=0; j<=kVector[1]; j++){
			for (int l=0; l<=kVector[2]; l++){
				if(i!=0||j!=0||l!=0){
					kValue a;
					a.k={i/box.boxSize()[0],j/box.boxSize()[1],l/box.boxSize()[2]};
					a.k*=2.0*M_PI;
					kValues.push_back(a);
				}




			}
		}
	}
	return kValues;
}


template <HPX_PROTOCOL hpx_protocol, typename Particle>
std::vector<kValue>& calcGlobalValues(std::vector<kValue>&kValues,Box<Particle>& box){

	auto calcKValues=[&box](auto& kValue){

		kValue.sin=0.0;
		kValue.cos=0.0;
		for (auto& cell:box.all()){
			for (auto& p:cell.particles()){
				kValue.sin+=p.charge()*sin(kValue.k*p.pos());
				kValue.cos+=p.charge()*cos(kValue.k*p.pos());
			}
		}
	};



	if constexpr (hpx_protocol == HPX_PROTOCOL::ASYNC) {
		std::vector<hpx::future<void>> fut;
		for (auto& kValue : kValues){
			fut.push_back(hpx::async(calcKValues,std::ref(kValue)));
		}
		hpx::wait_all(fut);
	} else if constexpr (hpx_protocol == HPX_PROTOCOL::NONE) {
		for (auto& kValue : kValues){
			calcKValues(kValue);
		}
	}
	return kValues;

}

template <HPX_PROTOCOL hpx_protocol, typename Particle>
void calckSpaceForces(std::vector<kValue>&kValues,Box<Particle>& box) {
	double alpha= box.alpha();
	double Vol = box.boxSize()[0]*box.boxSize()[1]*box.boxSize()[2];
	auto forcecalc=[alpha,Vol,kValues,coulombC=box.coulombC()](auto& p){

		for (auto kValue : kValues){
			double k2=kValue.k.norm2();
			double force1d=coulombC* p.charge()/Vol*4.0*M_PI/k2*exp(-k2/4.0/alpha/alpha);
			force1d*=-cos(kValue.k*p.pos())*kValue.sin+kValue.cos*sin(kValue.k*p.pos());
			p.force()+=force1d*kValue.k;
		}
	};

	if constexpr (hpx_protocol == HPX_PROTOCOL::ASYNC) {
		std::vector<hpx::future<void>> fut;

		for (auto& cell:box.all()){
			for (auto& p:cell.particles()){
				fut.push_back(hpx::async(forcecalc,std::ref(p)));
			}
		}
		hpx::wait_all(fut);
	} else if constexpr (hpx_protocol == HPX_PROTOCOL::NONE) {

		for (auto& cell:box.all()){
			for (auto& p:cell.particles()){
				forcecalc(p);
			}
		}
	}




}






template <typename Particle>
void calckSpaceForcesTest(std::vector<kValue>&kValues,Box<Particle>& box){
	Utils::Vector3d force;

	Utils::Vector3d sum;
	double alpha=box.alpha();
	bool test=false;
	double Vol = box.boxSize()[0]*box.boxSize()[1]*box.boxSize()[2];
	for (auto& cell:box.all()){
		for (auto& p:cell.particles()){
			force={0,0,0};
			test=true;
			for (auto& cell :  box.all()){
				for (auto& q:cell.particles()){
					sum={0,0,0};
					for (auto& kV:kValues){
						auto k =kV.k;
						sum+=4.0*M_PI*k/k.norm2()*exp(-k.norm2()/(4.0*alpha*alpha))*sin(k*(p.pos()-q.pos()));

					}
					force+=q.charge()*sum;
				}
			}
			force*=p.charge()/Vol;
			p.force()=force;			


			//if (test) break;
		}
		//if(test) break;
	}
	std::cout<<force<<std::endl;

}


template <HPX_PROTOCOL hpx_protocol, typename Particle>
void kSpaceForces(Box<Particle>& box){
	auto k = kSpaceSetup(box);
	k= calcGlobalValues<hpx_protocol>(k,box);
	calckSpaceForces<hpx_protocol>(k,box);
}

