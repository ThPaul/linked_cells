struct kValue{
	Utils::Vector3d k;
	double sin=0;
	double cos=0;


};

template <typename Particle>
void kSpaceForces(Utils::Vector3i kVector ,Box<Particle>& box){
	auto k = kSpaceSetup(kVector,box);
	k= calcGlobalValues(k,box);
	calckSpaceForces(k,box);

}


template <typename Particle>
std::vector<kValue> kSpaceSetup(Utils::Vector3i kVector ,Box<Particle>& box){
	std::vector<kValue> kValues;
	for (int i=0; i<=kVector[0]; i++){
		for (int j=0; j<=kVector[0]; j++){
			for (int l=0; l<=kVector[0]; l++){
				if(i!=0||j!=0||l!=0){
					kValue a;
					a.k={i/box.boxSize()[0],j/box.boxSize()[1],l/box.boxSize()[2]};
					a.k*=2.0*M_PI;
					kValues.push_back(a);
				}




			}
		}
	}
	std::cout<<kValues.size()<<std::endl;
	return kValues;
}


template <typename Particle>
std::vector<kValue>& calcGlobalValues(std::vector<kValue>&kValues,Box<Particle>& box){
	for (auto &kValue:kValues){
		kValue.sin=0.0;
		kValue.cos=0.0;
		for (auto& cell:box.all()){
			for (auto& p:cell.particles()){
				kValue.sin+=p.charge()*sin(kValue.k*p.pos());
				kValue.cos+=p.charge()*cos(kValue.k*p.pos());
			}
		}
	}
	return kValues;

}

template <typename Particle>
void calckSpaceForces(std::vector<kValue>&kValues,Box<Particle>& box){
	double alpha=0.666;
	double Vol = box.boxSize()[0]*box.boxSize()[1]*box.boxSize()[2];
		for (auto& cell:box.all()){
			for (auto& p:cell.particles()){
				for (auto& kValue : kValues){
					double k2=kValue.k.norm2();
					double force1d=box.coulombC()* p.charge()/Vol*4*M_PI/k2*exp(-k2/4.0/alpha/alpha);
					force1d*=cos(kValue.k*p.pos())*kValue.sin-kValue.cos*sin(kValue.k*p.pos());
					p.force()+=force1d*kValue.k;
					

					
				}
			}
		}
}
