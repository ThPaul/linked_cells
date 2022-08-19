struct InputData{
       Utils::Vector3d boxSize;
       double LJcutoff;
       double EwaldCutoff;
       BC boundaryCondition;
       double eps = 0;
       double sigma = 0;
       double coulombC =1;
       double ewaldacc = -1;
       std::string posFile;
       std::string chargeFile;
       double alpha;
       Utils::Vector3i kSpaceSize;
       std::vector<FORCES> forces;
};


InputData readInput(std::string file){
	InputData dat;


	std::ifstream settings (file);
	std::vector<Utils::Vector3d>particleCoords;
	std::vector<double>particleCharges;
	if(settings){
		std::string line;
		while(getline(settings,line)){
			if (line[0]!='#'){
			std::vector<std::string> split;
			boost::split(split,line,[](char c){return c=='=';});
			boost::trim(split[0]);
			boost::to_lower(split[0]);
			boost::trim(split[1]);
			if (split[0]=="ljcutoff") dat.LJcutoff=std::stod(split[1]);
			else if (split[0]=="ewaldcutoff") dat.EwaldCutoff=std::stod(split[1]);
			else if(split[0]=="epsilon") dat.eps=std::stod(split[1]);
			else if(split[0]=="sigma") dat.sigma=std::stod(split[1]);
			else if(split[0]=="coulombconstant") dat.coulombC=std::stod(split[1]);
			else if(split[0]=="positionfile") dat.posFile=split[1];
			else if(split[0]=="chargefile") dat.chargeFile=split[1];
			else if(split[0]=="alpha") dat.alpha=std::stod(split[1]);
			else if(split[0]=="ewaldaccurracy") dat.ewaldacc=std::stod(split[1]);
			else if(split[0]=="boundarycondition"){
				boost::to_lower(split[1]);
				if(split[1]=="none") dat.boundaryCondition=BC::NONE;
				else  dat.boundaryCondition=BC::PERIODIC; 
			}
			else if (split[0]=="boxsize"){

			std::vector<std::string> dim;
			boost::split(dim,split[1],[](char c){return c==',';});
			Utils::Vector3d coords;
			for (int i : {0,1,2}){
				coords[i]=std::stod(dim[i]);
				}
			dat.boxSize=coords;
			}

			else if (split[0]=="kspacesize"){

			std::vector<std::string> dim;
			boost::split(dim,split[1],[](char c){return c==',';});
			Utils::Vector3i coords;
			for (int i : {0,1,2}){
				coords[i]=std::stoi(dim[i]);
				}
			dat.kSpaceSize=coords;
			}

			else if (split[0]=="forces"){

			std::vector<std::string> forces;
			boost::split(forces,split[1],[](char c){return c==',';});
			for (auto force:forces){
				if (force=="LJ") dat.forces.push_back(FORCES::LJ);
				else if (force=="EWALDRS") dat.forces.push_back(FORCES::EWALDRS);
				else if (force=="EWALDKS") dat.forces.push_back(FORCES::EWALDKS);
				else std::cout << force << " is NOT an available force" << std::endl;
			}
			}

			else std::cout<<"INVALID :" <<line<<std::endl;
			}
		}
		settings.close();



	}
	return dat;
}
