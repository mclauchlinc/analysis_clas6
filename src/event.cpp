#include "event.hpp"


Event::Event(std::shared_ptr<Branches> data_, std::shared_ptr<Histogram> hist_, int data_set_, bool sim_){
	int num_par = data_->Branches::gpart();
	//Create particle candidates
	Particle elec_candidate; 
	Particle had_candidate[num_par-1];
	//Particle indexes of the friends we have made
	int pro_idx[20];
	int pip_idx[20];
	int pim_idx[20];
	int pro_pip_idx[20];

	bool top_pos[4] = {false,false,false,false};
	bool top_good[4] = {false,false,false,false};

	switch(data_set_){
		case 0:
			_k_mu = k_mu_e16;
		break;
		case 1:
			_k_mu = k_mu_e1f;
		break;
		default:
			_k_mu = k_mu_e16;
		break;
	}


	for(int i = 0; i< 20; i++){
		pro_idx[i] = -99;
		pip_idx[i] = -99;
		pim_idx[i] = -99;
		pro_pip_idx[i] = -99;
	}

	//Look into Electron candidate first
	elec_candidate.Particle::Fill_Particle(data_,0);

	//Particle Identification
	elec_candidate.Particle::EID_Cuts();//hist_); //EID Cuts and Plottting 
	if(elec_candidate.Particle::Is_Elec()){//without a good electron delta t is unreliable so hadron ID is subsequently unreliable
		_elec = elec_candidate;//Assign the electron to the event as good
		for(int h = 1; h<num_par; h++){//First Hadron Loop for initial tagging
			had_candidate[h-1].Particle::Fill_Particle(data_,h);
			had_candidate[h-1].Particle::HID_Cuts();//hist_);
			if(had_candidate[h-1].Particle::Is_Pro()){
				_meas_had[0]++;
				pro_idx[_meas_had[0]-1]= h-1;
			}else if(had_candidate[h-1].Particle::Is_Pip()){
				_meas_had[1]++;
				pip_idx[_meas_had[1]-1]= h-1;
			}else if(had_candidate[h-1].Particle::Is_Pim()){
				_meas_had[2]++;
				pim_idx[_meas_had[2]-1]= h-1;
			}else if(had_candidate[h-1].Particle::Pro_Pip_Prob()){
				_meas_had[3]++;
				pro_pip_idx[_meas_had[3]-1]= h-1;
			}
		}
		//Fixing if a particle has been double identified
		if(_meas_had[3] > 0){
			for(int i = 0; i< _meas_had[3] ; i++){
				if(_meas_had[0] == 1 && _meas_had[2] == 1){//If only one proton and pim then we know it has to be a proton
					had_candidate[pro_pip_idx[i]].Particle::Pro_Pip_Resolved(1);//Assign that hadron to be a pip
					_meas_had[3] = _meas_had[3]-1;
				}else if(_meas_had[1] == 1 && _meas_had[2] == 1){//If only one pip and pim then we know it has to be a proton
					had_candidate[pro_pip_idx[i]].Particle::Pro_Pip_Resolved(0);//Assign that hadron to be a proton
					_meas_had[3] = _meas_had[3]-1;
				}else{//More nuanced case where we have multiple particles doing their thing
					//Pro Missing top 
					if(_meas_had[2]> 0){//Pim measured
						for(int j = 0; j<_meas_had[2]; j++){
							if(_meas_had[1]>0){//Pip Measured
								for(int k = 0; k<_meas_had[1]; k++){
									if(Selection::Event_Selection(_k_mu,_elec,had_candidate[pip_idx[k]],had_candidate[pim_idx[j]]) == false){//If we didn't find another PIP to satisfy the event criteria
										if(Selection::Event_Selection(_k_mu,_elec,had_candidate[pro_pip_idx[i]],had_candidate[pim_idx[j]])){//If the double IDed particle satisfies the event as a pip
											had_candidate[pro_pip_idx[i]].Particle::Pro_Pip_Resolved(1);//Then let it be a pip
											_meas_had[3] = _meas_had[3]-1;
											top_good[0] = true;
											_pip = had_candidate[pro_pip_idx[i]];
											_pim = had_candidate[pim_idx[j]];
										}
									}
								}
							}else{//No measured Pip
								if(Selection::Event_Selection(_k_mu,_elec,had_candidate[pro_pip_idx[i]],had_candidate[pim_idx[j]])){
									had_candidate[pro_pip_idx[i]].Particle::Pro_Pip_Resolved(1);//Then let it be a pip
									_meas_had[3] = _meas_had[3]-1;
									top_good[0] = true;
									_pip = had_candidate[pro_pip_idx[i]];
									_pim = had_candidate[pim_idx[j]];
								}
							}
						}
					}
				}
			}
		}else{//No double ID of particles! Yay! 

		}
	} 
	
}

Event::~Event(){ }