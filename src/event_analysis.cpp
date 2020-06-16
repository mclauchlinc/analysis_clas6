#include "event_analysis.hpp"
		  
Analysis::Analysis(std::shared_ptr<Branches> data_, std::shared_ptr<Histogram> hist_, std::shared_ptr<Environment> envi_, int run_type_){
	switch(run_type_){
		case 1:
			_set = 1;
			_recon = true;
		break;
		case 2:
			_set = 0;
			_recon = true;
		break;
		case 3:
			_set = 1;
			_sim = true;
			_thrown = true;
			_recon = true;
		break;
		case 4:
			_set = 0;
			_sim = true;
			_thrown = true;
			_recon = true;
		break;
		case 5:
			_set = 1; 
			_sim = true;
			_thrown = true;
		break;
		case 6:
			_set = 0;
			_sim = true;
			_thrown = true;
		break;
	}
	//std::cout<<"Set: " <<_set <<"	| Sim: " <<_sim <<"	| Thrown: " <<_thrown <<"	| Recon: " <<_recon <<std::endl;

	if(_sim){
		_npart = data_->Branches::npart();
		//std::cout<<"confirm is sim" <<std::endl;
	}else{
		_npart = data_->Branches::gpart();
		//std::cout<<"confirm not sim" <<std::endl;
	}

	//Not sure if I actually need this, because I think I'm doing everything already inside the events? 
	if(_thrown){
		//std::cout<<"Making Thrown W and Q2: " ;
		_W[1] = physics::WP(run_type_,data_,1);
		_Q2[1] = physics::Qsquared(run_type_,data_,1);
		//std::cout<<_W[1] <<" " <<_Q2[1] <<std::endl;
	}
	if(_recon){
		_W[0] = physics::WP(run_type_,data_,0);
		_Q2[0] = physics::Qsquared(run_type_,data_,0);
	}

	if(_sim){//Thrown analysis
		//std::cout<<"Filling sim particles" <<std::endl;
		//Particle Loop
		Particle tPart[4]; 
		Event tEvent[4];
		for(int j = 0; j< 4; j++){
			tPart[j].Particle::Fill_Particle(data_,j,_set,_sim,true);
			tPart[j].Particle::PID(data_,envi_,hist_,_W[1]);
			//tPart[j].Particle::Check_Particle();
		}
		//Topology Assignment and Particle Histogram Filling
		for(int k = 0; k< 4; k++){
			if(tPart[k].Particle::Is_Elec()){
				for(int l = 0; l < 4; l++){
					if(l!= k && tPart[l].Particle::Is_Pip()){//Pro Miss
						for(int m = 0; m< 4; m++){
							if(m!=l && m!=k && tPart[m].Particle::Is_Pim()){
								tEvent[0].Event::Fill_Event(envi_,hist_,0,_W[1],_Q2[1],tPart[k],tPart[l],tPart[m]);
							}
						}
					}
					if(l!=k && tPart[l].Particle::Is_Pro()){//
						for(int n = 0; n< 4; n++){
							if(n!=l && n!=k && tPart[n].Particle::Is_Pim()){//Pip Miss
								tEvent[1].Event::Fill_Event(envi_,hist_,1,_W[1],_Q2[1],tPart[k],tPart[l],tPart[n]);
							}else if(n!=l && n!=k && tPart[n].Particle::Is_Pip()){//Pim && Zero Miss
								tEvent[2].Event::Fill_Event(envi_,hist_,2,_W[1],_Q2[1],tPart[k],tPart[l],tPart[n]);
								for(int o = 0; o<4; o++){
									if(o!=l && o!=k && o!=n && tPart[o].Particle::Is_Pim()){//Zero Miss
										tEvent[3].Event::Fill_Event(envi_,hist_,3,_W[1],_Q2[1],tPart[k],tPart[l],tPart[n],tPart[o]);
									}
								}
							}
						}
					}
				}
			}
		}
		for(int x = 0; x<4; x++){
			tEvent[x].Event::Assign_Weight(data_->Branches::weight()/4.0);
		}
		//*filling event Histogram here*
	}

	if(_recon){//Experimental or reconstructed data
		//std::cout<<"Entered Non-Sim Event" <<std::endl;
		//Particle Loop
		int _nevnt_ = -1; 
		int _evnt_idx_ = 0;
		int mix = 0; 
		Particle ePart[_npart]; 
		for(int a = 0; a< _npart; a++){
			ePart[a].Particle::Fill_Particle(data_,a,_set,_sim);
			ePart[a].Particle::PID(data_,envi_,hist_,_W[0]);
			//ePart[a].Particle::Check_Particle();
		}
		//Figure out how many possible events we have
		for(int w = 1; w<_npart; w++){
			if(ePart[w].Particle::Is_Elec()){
				_gpart[0]+=1;
				//std::cout<<"	Particle " <<w <<" is Elec" <<std::endl;
			}
			if(ePart[w].Particle::Is_Pro()){
				_gpart[1]+=1;

				//std::cout<<"	Particle " <<w <<" is Pro" <<std::endl;
			}
			if(ePart[w].Particle::Is_Pip()){
				_gpart[2]+=1;
				//std::cout<<"	Particle " <<w <<" is PIP" <<std::endl;
			}
			if(ePart[w].Particle::Is_Pim()){
				_gpart[3]+=1;
				//std::cout<<"	Particle " <<w <<" is PIM" <<std::endl;
			}
			if(ePart[w].Particle::ID_crisis() == 1){
				mix+=1;
			}
		}
		_nevnt_ = (_gpart[1]*_gpart[3]+_gpart[2]*_gpart[3]+(_gpart[1]*_gpart[2]-mix)*(1+_gpart[3]))*ePart[0].Particle::Is_Elec();
		//if(ePart[0].Particle::Is_Elec()){
			//std::cout<<"Ele: " <<ePart[0].Particle::Is_Elec() <<"	| Pro: " <<_gpart[1] <<"	| Pip: " <<_gpart[2] <<"	| Pim: " <<_gpart[3] <<"	| Mix: " <<mix <<std::endl;
			//std::cout<<"Possible events: " <<_nevnt_ <<std::endl;
		//}
		if(_nevnt_ > 0){
			Event eEvent[_nevnt_];

			//Event Assignment into Topologies for Thrown
			//I know this part is inefficient, but it does work, and I assume I'm not going to have an absurd number of events like this
			//std::cout<<"	is there a good electron? " <<ePart[0].Particle::Is_Elec() <<std::endl <<std::endl;
			if(ePart[0].Particle::Is_Elec()){//Assume that only the first particle can be our event electron
				for(int g = 1; g < _npart; g++){
					if(ePart[g].Particle::Is_Pip()){//Pro Miss
						for(int h = 1; h< _npart; h++){
							if(h!=g && ePart[h].Particle::Is_Pim()){
								//std::cout<<"		Event Combo: "<<_evnt_idx_ <<" top: Pro Miss"<<std::endl;
								//std::cout<<"		Pre Event Fill for event Combo: " <<_evnt_idx_ <<std::endl;
								//ePart[0].Particle::Check_Particle();
								//ePart[g].Particle::Check_Particle();
								//ePart[h].Particle::Check_Particle();
								eEvent[_evnt_idx_].Event::Fill_Event(envi_,hist_,0,_W[0],_Q2[0],ePart[0],ePart[g],ePart[h]);
								//std::cout<<"		Post Event Fill for event Combo: " <<_evnt_idx_ <<std::endl;
								//std::cout<<"			Elec " <<ePart[0].Particle::Get_idx() <<" pid: " <<ePart[0].Particle::Is_Elec() <<"| Pip " <<ePart[g].Particle::Get_idx() <<" pid: " <<ePart[g].Particle::Is_Pip() <<"| Pim " <<ePart[h].Particle::Get_idx() <<" pid: " <<ePart[h].Particle::Is_Pim() <<std::endl;
								//ePart[0].Particle::Check_Particle();
								//ePart[g].Particle::Check_Particle();
								//ePart[h].Particle::Check_Particle();
								_evnt_idx_+=1;
								_mtop[0]+=1;
							}
						}
					}
					if(ePart[g].Particle::Is_Pro()){//
						for(int s = 1; s< _npart; s++){
							if(s!=g && ePart[s].Particle::Is_Pim()){//Pip Miss
								//std::cout<<"		Event Combo: " <<_evnt_idx_ <<" top: PIP Miss"<<std::endl;
								//std::cout<<"		Pre Event Fill for event Combo: " <<_evnt_idx_ <<std::endl;
								//std::cout<<"			Elec " <<ePart[0].Particle::Get_idx() <<" pid: " <<ePart[0].Particle::Is_Sim() <<"| Pro " <<ePart[g].Particle::Get_idx() <<" pid: " <<ePart[g].Particle::Is_Pro() <<"| Pim " <<ePart[s].Particle::Get_idx() <<" pid: " <<ePart[s].Particle::Is_Pim() <<std::endl;
								//ePart[0].Particle::Check_Particle();
								//ePart[g].Particle::Check_Particle();
								//ePart[s].Particle::Check_Particle();
								eEvent[_evnt_idx_].Event::Fill_Event(envi_,hist_,1,_W[0],_Q2[0],ePart[0],ePart[g],ePart[s]);
								//std::cout<<"		Post Event Fill for event Combo: " <<_evnt_idx_ <<std::endl;
								//std::cout<<"			Elec " <<ePart[0].Particle::Get_idx() <<" pid: " <<ePart[0].Particle::Is_Sim() <<"| Pro " <<ePart[g].Particle::Get_idx() <<" pid: " <<ePart[g].Particle::Is_Pro() <<"| Pim " <<ePart[s].Particle::Get_idx() <<" pid: " <<ePart[s].Particle::Is_Pim() <<std::endl;
								//ePart[0].Particle::Check_Particle();
								//ePart[g].Particle::Check_Particle();
								//ePart[s].Particle::Check_Particle();
								_evnt_idx_+=1;
								_mtop[1]+=1;
							}else if(s!=g && ePart[s].Particle::Is_Pip()){//Pim
								//std::cout<<"		Event Combo: "<<_evnt_idx_ <<" top: PIM Miss"<<std::endl;
								//std::cout<<"		Pre Event Fill for event Combo: " <<_evnt_idx_ <<std::endl;
								//std::cout<<"			Elec " <<ePart[0].Particle::Get_idx() <<" pid: " <<ePart[0].Particle::Is_Sim() <<"| Pro " <<ePart[g].Particle::Get_idx() <<" pid: " <<ePart[g].Particle::Is_Pro() <<"| Pip " <<ePart[s].Particle::Get_idx() <<" pid: " <<ePart[s].Particle::Is_Pip() <<std::endl;
								//ePart[0].Particle::Check_Particle();
								//ePart[g].Particle::Check_Particle();
								//ePart[s].Particle::Check_Particle();
									 							    //(envi_,hist_,top_,W_,Q2_,&p1,&p2,&p3,hel_)
								eEvent[_evnt_idx_].Event::Fill_Event(envi_,hist_,2,_W[0],_Q2[0],ePart[0],ePart[g],ePart[s]);
								//std::cout<<"		Post Event Fill for event Combo: " <<_evnt_idx_ <<std::endl;
								//std::cout<<"			Elec " <<ePart[0].Particle::Get_idx() <<" pid: " <<ePart[0].Particle::Is_Sim() <<"| Pro " <<ePart[g].Particle::Get_idx() <<" pid: " <<ePart[g].Particle::Is_Pro() <<"| Pip " <<ePart[s].Particle::Get_idx() <<" pid: " <<ePart[s].Particle::Is_Pip() <<std::endl;
								//ePart[g].Particle::Check_Particle();
								//ePart[0].Particle::Check_Particle();
								//ePart[s].Particle::Check_Particle();
								_evnt_idx_+=1;
								_mtop[2]+=1;
								for(int u = 1; u<_npart; u++){
									if(u!=g && u!=s && ePart[u].Particle::Is_Pim()){//Zero Miss
										//std::cout<<"		Event Combo: "<<_evnt_idx_ <<" top:zero"<<std::endl;
										//std::cout<<"		Pre Event Fill for event Combo: " <<_evnt_idx_ <<std::endl;
										//std::cout<<"			Elec " <<ePart[0].Particle::Get_idx() <<" pid: " <<ePart[0].Particle::Is_Sim() <<"| Pro " <<ePart[g].Particle::Get_idx() <<" pid: " <<ePart[g].Particle::Is_Pro() <<"| Pip " <<ePart[s].Particle::Get_idx() <<" pid: " <<ePart[s].Particle::Is_Pip() <<"| Pim " <<ePart[u].Particle::Get_idx() <<" pid: " <<ePart[u].Particle::Is_Pim() <<std::endl;									
										//ePart[0].Particle::Check_Particle();
										//ePart[g].Particle::Check_Particle();
										//ePart[s].Particle::Check_Particle();
										//ePart[u].Particle::Check_Particle();
										eEvent[_evnt_idx_].Event::Fill_Event(envi_,hist_,3,_W[0],_Q2[0],ePart[0],ePart[g],ePart[s],ePart[u]);
										//std::cout<<"		Post Event Fill for event Combo: " <<_evnt_idx_ <<std::endl;
										//std::cout<<"			Elec " <<ePart[0].Particle::Get_idx() <<" pid: " <<ePart[0].Particle::Is_Sim() <<"| Pro " <<ePart[g].Particle::Get_idx() <<" pid: " <<ePart[g].Particle::Is_Pro() <<"| Pip " <<ePart[s].Particle::Get_idx() <<" pid: " <<ePart[s].Particle::Is_Pip() <<"| Pim " <<ePart[u].Particle::Get_idx() <<" pid: " <<ePart[u].Particle::Is_Pim() <<std::endl;
										//ePart[0].Particle::Check_Particle();
										//ePart[g].Particle::Check_Particle();
										//ePart[s].Particle::Check_Particle();
										//ePart[u].Particle::Check_Particle();
										_evnt_idx_+=1;
										_mtop[3]+=1;
									}
								}
							}
						}
					}
				}
			}
			//Event Counting!
			for(int ev = 0; ev < _evnt_idx_; ev++){

				if(eEvent[ev].Event::Gevnt()){
					_gevts = _gevts + 1.0;
					eEvent[ev].Event::COM_4Vec();//Get 4-momenta into COM frame
					eEvent[ev].Event::Vars();//Get alpha, theta, and MM for given expressions
					if(eEvent[ev].Event::Top(0)){
						_ntop[0]+=1;
					}else if(eEvent[ev].Event::Top(1)){
						_ntop[1]+=1;
					}else if(eEvent[ev].Event::Top(2)){
						_ntop[2]+=1;
					}else if(eEvent[ev].Event::Top(3)){
						_ntop[3]+=1;
					}
				}
			}
			//Assign Event Weights
			for(int ev2 = 0; ev2 < _evnt_idx_; ev2++){
				if(eEvent[ev2].Event::Gevnt()){
					if(_sim){
						eEvent[ev2].Event::Assign_Weight(data_->Branches::weight()/_gevts);
					}else{
						eEvent[ev2].Event::Assign_Weight(1.0/_gevts);
					}
				}
				//*Fill Event Histogram here*
			}
		}	
	}
}


int Analysis::Gevts(){
	return _gevts;
}
int Analysis::Ntop(int i){
	return _ntop[i];
}

float Analysis::Get_W(int i){
	return _W[i];
}
float Analysis::Get_Q2(int i){
	return _Q2[i];
}
int Analysis::Get_set(){
	return _set;
}
bool Analysis::Is_Sim(){
	return _sim;
}
bool Analysis::Is_Thrown(){
	return _thrown;
}