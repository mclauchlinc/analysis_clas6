#include "event_class.hpp"


Event_Class::Event_Class(std::shared_ptr<Branches> data, std::shared_ptr<Histogram> _hists, int run_type, int plate_info, std::shared_ptr<Environment> envi ){ 
						 
	bool fid_e_pass = false;
	bool sf_e_pass = false;
	bool cc_e_pass = false;
	int num_parts = data->Branches::gpart();
	float theta[num_parts];
	float phi[num_parts];
	float sector[num_parts];
	float p_[num_parts];
	bool in_range = false;

	_helicity = physics::event_helicity(data,plate_info);
	_run_type = run_type;

	switch(run_type%2){
		case 1:
			_beam = k_mu_e16;
		break;
		case 0:
			_beam = k_mu_e1f;
		break;
		default:
			_beam = k_mu_e16;
		break;
	}
	
	_target = p_mu;

	//std::cout<<std::endl <<"gpart: " <<num_parts <<std::endl;
	for(int i=0; i < num_parts; i++){ 
		//hadron[i]->Particle::Fill_Particle(data,i);
					theta[i] = physics::get_theta(data->Branches::cz(i));
					phi[i] = physics::get_phi(data->Branches::cx(i),data->Branches::cy(i));
					sector[i] = physics::get_sector(phi[i]);
					p_[i] = data->Branches::p(i);
		//std::cout<<"theta " <<i <<": " <<theta[i] <<std::endl;
		//std::cout<<"phi " <<i <<": " <<phi[i]<<std::endl;
		//std::cout<<"sector " <<i <<": " <<sector[i]<<std::endl;
		//std::cout<<"momentum " <<i <<": " <<_p[i]<<std::endl;
	}	


	_W = physics::WP(0,data);//The first index {0,1} -> {e16, e1f}
	_Q2 = physics::Qsquared(0,data);//The first index {0,1} -> {e16, e1f}

	/*if(_W >= Wmin && _W <= Wmax ){//Checking to see if the particle is in the relevant W Q2 region 
		if(_Q2 >= Q2min && _Q2 <= Q2max){
			in_range = true;
		}
	}*/
	//Pre ID Filling
		//Event-Wide
		//Electrons
	_hists->Histogram::WQ2_Fill(envi,0,0,_W,_Q2);
	_hists->Histogram::Fid_Fill(envi,0,physics::get_theta(data->Branches::cz(0)),physics::get_phi(data->Branches::cx(0),data->Branches::cy(0)),0,0,0,_W,data->Branches::p(0));
	_hists->Histogram::SF_Fill(envi,0,data->Branches::p(0),data->Branches::etot(0),0,0,_W,sector[0]);
	_hists->Histogram::CC_Fill(envi,0,data->Branches::cc_sect(0),data->Branches::cc_segm(0),data->Branches::nphe(0),0,0);
	
	//electron ID
	//if(in_range && data->Branches::q(0)==-1 && data->Branches::cc(0)!=0 && data->Branches::dc(0)!=0 && data->Branches::sc(0)!=0 && data->Branches::ec(0)!=0){//Sanity are q = -1, cc, dc, ec, and sc hits
	if(cuts::in_range(_W,_Q2,envi) && cuts::e_sanity(data,envi)){
		_hists->Histogram::WQ2_Fill(envi,0,1,_W,_Q2);
		_hists->Histogram::Fid_Fill(envi,0,physics::get_theta(data->Branches::cz(0)),physics::get_phi(data->Branches::cx(0),data->Branches::cy(0)),0,1,0,_W,data->Branches::p(0));
		_hists->Histogram::SF_Fill(envi,0,data->Branches::p(0),data->Branches::etot(0),1,0,_W,sector[0]);
		_hists->Histogram::CC_Fill(envi,0,data->Branches::cc_sect(0),data->Branches::cc_segm(0),data->Branches::nphe(0),1,0);
		//if(cuts::fid_cut(0,data->Branches::p(0),data->Branches::cx(0),data->Branches::cy(0),data->Branches::cz(0))){
		if(cuts::e_fid(data,envi)){
			fid_e_pass = true;
			_hists->Histogram::WQ2_Fill(envi,0,2,_W,_Q2);
			_hists->Histogram::Fid_Fill(envi,0,physics::get_theta(data->Branches::cz(0)),physics::get_phi(data->Branches::cx(0),data->Branches::cy(0)),0,2,0,_W,data->Branches::p(0));
			_hists->Histogram::SF_Fill(envi,0,data->Branches::p(0),data->Branches::etot(0),2,0,_W,sector[0]);
			_hists->Histogram::CC_Fill(envi,0,data->Branches::cc_sect(0),data->Branches::cc_segm(0),data->Branches::nphe(0),2,0);
		}else{
			_hists->Histogram::Fid_Fill(envi,0,physics::get_theta(data->Branches::cz(0)),physics::get_phi(data->Branches::cx(0),data->Branches::cy(0)),0,2,1,_W,data->Branches::p(0));
			_hists->Histogram::SF_Fill(envi,0,data->Branches::p(0),data->Branches::etot(0),2,1,_W,sector[0]);
			_hists->Histogram::CC_Fill(envi,0,data->Branches::cc_sect(0),data->Branches::cc_segm(0),data->Branches::nphe(0),2,1);
		}
		//if(cuts::sf_cut(data->Branches::p(0),data->Branches::etot(0),data->Branches::cx(0),data->Branches::cy(0))){
		if(cuts::e_sf(data,envi)){
			sf_e_pass = true;
			_hists->Histogram::WQ2_Fill(envi,0,3,_W,_Q2);
			_hists->Histogram::Fid_Fill(envi,0,physics::get_theta(data->Branches::cz(0)),physics::get_phi(data->Branches::cx(0),data->Branches::cy(0)),0,3,0,_W,data->Branches::p(0));
			_hists->Histogram::SF_Fill(envi,0,data->Branches::p(0),data->Branches::etot(0),3,0,_W,sector[0]);
			_hists->Histogram::CC_Fill(envi,0,data->Branches::cc_sect(0),data->Branches::cc_segm(0),data->Branches::nphe(0),3,0);
		}else{
			_hists->Histogram::Fid_Fill(envi,0,physics::get_theta(data->Branches::cz(0)),physics::get_phi(data->Branches::cx(0),data->Branches::cy(0)),0,3,1,_W,data->Branches::p(0));
			_hists->Histogram::SF_Fill(envi,0,data->Branches::p(0),data->Branches::etot(0),3,1,_W,sector[0]);
			_hists->Histogram::CC_Fill(envi,0,data->Branches::cc_sect(0),data->Branches::cc_segm(0),data->Branches::nphe(0),3,1);
		}
		//if(cuts::min_cc(data->Branches::cc_segm(0),data->Branches::cc_sect(0),data->Branches::nphe(0))){
		if(cuts::e_cc(data,envi)){
			cc_e_pass = true;
			_hists->Histogram::WQ2_Fill(envi,0,4,_W,_Q2);
			_hists->Histogram::Fid_Fill(envi,0,physics::get_theta(data->Branches::cz(0)),physics::get_phi(data->Branches::cx(0),data->Branches::cy(0)),0,4,0,_W,data->Branches::p(0));
			_hists->Histogram::SF_Fill(envi,0,data->Branches::p(0),data->Branches::etot(0),4,0,_W,sector[0]);
			_hists->Histogram::CC_Fill(envi,0,data->Branches::cc_sect(0),data->Branches::cc_segm(0),data->Branches::nphe(0),4,0);
		}else{
			_hists->Histogram::Fid_Fill(envi,0,physics::get_theta(data->Branches::cz(0)),physics::get_phi(data->Branches::cx(0),data->Branches::cy(0)),0,4,1,_W,data->Branches::p(0));
			_hists->Histogram::SF_Fill(envi,0,data->Branches::p(0),data->Branches::etot(0),4,1,_W,sector[0]);
			_hists->Histogram::CC_Fill(envi,0,data->Branches::cc_sect(0),data->Branches::cc_segm(0),data->Branches::nphe(0),4,1);
		}
		if(fid_e_pass && sf_e_pass){
			_hists->Histogram::WQ2_Fill(envi,0,5,_W,_Q2);
			_hists->Histogram::Fid_Fill(envi,0,physics::get_theta(data->Branches::cz(0)),physics::get_phi(data->Branches::cx(0),data->Branches::cy(0)),0,5,0,_W,data->Branches::p(0));
			_hists->Histogram::SF_Fill(envi,0,data->Branches::p(0),data->Branches::etot(0),5,0,_W,sector[0]);
			_hists->Histogram::CC_Fill(envi,0,data->Branches::cc_sect(0),data->Branches::cc_segm(0),data->Branches::nphe(0),5,0);
		}else{
			_hists->Histogram::Fid_Fill(envi,0,physics::get_theta(data->Branches::cz(0)),physics::get_phi(data->Branches::cx(0),data->Branches::cy(0)),0,5,1,_W,data->Branches::p(0));
			_hists->Histogram::SF_Fill(envi,0,data->Branches::p(0),data->Branches::etot(0),5,1,_W,sector[0]);
			_hists->Histogram::CC_Fill(envi,0,data->Branches::cc_sect(0),data->Branches::cc_segm(0),data->Branches::nphe(0),5,1);
		}
		if(fid_e_pass && cc_e_pass){
			_hists->Histogram::WQ2_Fill(envi,0,6,_W,_Q2);
			_hists->Histogram::Fid_Fill(envi,0,physics::get_theta(data->Branches::cz(0)),physics::get_phi(data->Branches::cx(0),data->Branches::cy(0)),0,6,0,_W,data->Branches::p(0));
			_hists->Histogram::SF_Fill(envi,0,data->Branches::p(0),data->Branches::etot(0),6,0,_W,sector[0]);
			_hists->Histogram::CC_Fill(envi,0,data->Branches::cc_sect(0),data->Branches::cc_segm(0),data->Branches::nphe(0),6,0);
		}else{
			_hists->Histogram::Fid_Fill(envi,0,physics::get_theta(data->Branches::cz(0)),physics::get_phi(data->Branches::cx(0),data->Branches::cy(0)),0,6,1,_W,data->Branches::p(0));
			_hists->Histogram::SF_Fill(envi,0,data->Branches::p(0),data->Branches::etot(0),6,1,_W,sector[0]);
			_hists->Histogram::CC_Fill(envi,0,data->Branches::cc_sect(0),data->Branches::cc_segm(0),data->Branches::nphe(0),6,1);
		}
		if(sf_e_pass && cc_e_pass){//SF and Min CC
			_hists->Histogram::WQ2_Fill(envi,0,7,_W,_Q2);
			_hists->Histogram::Fid_Fill(envi,0,physics::get_theta(data->Branches::cz(0)),physics::get_phi(data->Branches::cx(0),data->Branches::cy(0)),0,7,0,_W,data->Branches::p(0));
			_hists->Histogram::SF_Fill(envi,0,data->Branches::p(0),data->Branches::etot(0),7,0,_W,sector[0]);
			_hists->Histogram::CC_Fill(envi,0,data->Branches::cc_sect(0),data->Branches::cc_segm(0),data->Branches::nphe(0),7,0);
		}else{
			_hists->Histogram::Fid_Fill(envi,0,physics::get_theta(data->Branches::cz(0)),physics::get_phi(data->Branches::cx(0),data->Branches::cy(0)),0,7,1,_W,data->Branches::p(0));
			_hists->Histogram::SF_Fill(envi,0,data->Branches::p(0),data->Branches::etot(0),7,1,_W,sector[0]);
			_hists->Histogram::CC_Fill(envi,0,data->Branches::cc_sect(0),data->Branches::cc_segm(0),data->Branches::nphe(0),7,1);
		}
		if(cuts::eid(data,envi)){//EID
			_hists->Histogram::WQ2_Fill(envi,0,8,_W,_Q2);
			_hists->Histogram::Fid_Fill(envi,0,physics::get_theta(data->Branches::cz(0)),physics::get_phi(data->Branches::cx(0),data->Branches::cy(0)),0,8,0,_W,data->Branches::p(0));
			_hists->Histogram::SF_Fill(envi,0,data->Branches::p(0),data->Branches::etot(0),8,0,_W,sector[0]);
			_hists->Histogram::CC_Fill(envi,0,data->Branches::cc_sect(0),data->Branches::cc_segm(0),data->Branches::nphe(0),8,0);
			_elec = physics::Make_4Vector(data->Branches::p(0),data->Branches::cx(0),data->Branches::cy(0),data->Branches::cz(0),me);
			good_electron++;
		}else{
			_hists->Histogram::Fid_Fill(envi,0,physics::get_theta(data->Branches::cz(0)),physics::get_phi(data->Branches::cx(0),data->Branches::cy(0)),0,8,1,_W,data->Branches::p(0));
			_hists->Histogram::SF_Fill(envi,0,data->Branches::p(0),data->Branches::etot(0),8,1,_W,sector[0]);
			_hists->Histogram::CC_Fill(envi,0,data->Branches::cc_sect(0),data->Branches::cc_segm(0),data->Branches::nphe(0),8,1);
		}
		if((data->Branches::id(0))==ELECTRON){
			_hists->Histogram::WQ2_Fill(envi,0,9,_W,_Q2);
			_hists->Histogram::Fid_Fill(envi,0,physics::get_theta(data->Branches::cz(0)),physics::get_phi(data->Branches::cx(0),data->Branches::cy(0)),0,9,0,_W,data->Branches::p(0));
			_hists->Histogram::SF_Fill(envi,0,data->Branches::p(0),data->Branches::etot(0),9,0,_W,sector[0]);
			_hists->Histogram::CC_Fill(envi,0,data->Branches::cc_sect(0),data->Branches::cc_segm(0),data->Branches::nphe(0),9,0);
		}
		else{
			_hists->Histogram::Fid_Fill(envi,0,physics::get_theta(data->Branches::cz(0)),physics::get_phi(data->Branches::cx(0),data->Branches::cy(0)),0,9,1,_W,data->Branches::p(0));
			_hists->Histogram::SF_Fill(envi,0,data->Branches::p(0),data->Branches::etot(0),9,1,_W,sector[0]);
			_hists->Histogram::CC_Fill(envi,0,data->Branches::cc_sect(0),data->Branches::cc_segm(0),data->Branches::nphe(0),9,1);
		}
	}else{
		_hists->Histogram::Fid_Fill(envi,0,physics::get_theta(data->Branches::cz(0)),physics::get_phi(data->Branches::cx(0),data->Branches::cy(0)),0,1,1,_W,data->Branches::p(0));
		_hists->Histogram::SF_Fill(envi,0,data->Branches::p(0),data->Branches::etot(0),1,1,_W,sector[0]);
		_hists->Histogram::CC_Fill(envi,0,data->Branches::cc_sect(0),data->Branches::cc_segm(0),data->Branches::nphe(0),1,1);
	}
	
	//Hadron Loop
	int q_h = 0; //desired charge for hadron
	bool fid_h = false; 
	bool dt_h = false;
	int hid_num = -99; 
	int p_loop = 0;
	int pip_loop = 0;
	int pim_loop = 0; 

	if(cuts::in_range(_W,_Q2,envi)){
		for(int h = 1; h < num_parts; h++){//Particle Loop
			for(int part = 0; part < 3; part++){//Hadron loop
				fid_h = false; 
				dt_h = false; 
				hid_num = -99;
				//Assign desired charge
				switch(part){
					case 0:
					q_h = 1;
					hid_num = PROTON;
					break;
					case 1:
					q_h = 1; 
					hid_num = PION; 
					break;
					case 2:
					q_h = -1;
					hid_num = -PION;
					break;
				}
				//std::cout<<std::endl <<"for part = "<<part <<" q_h = " <<q_h <<" and hid_num = " <<hid_num;
				//Pre Cut
				_hists->Histogram::Fid_Fill(envi,0,theta[h],phi[h],part+1,0,0,_W,data->Branches::p(h));
				_hists->Histogram::DT_Fill(envi,0,part,data->Branches::p(h), data->Branches::sc_r(h), data->Branches::sc_t(h), data->Branches::sc_r(0), data->Branches::sc_t(0),0,0,_W,sector[h]);
				//Sanity Cut
				//if(data->Branches::q(h)==q_h && data->Branches::dc(h)!=0 && data->Branches::sc(h)!=0){
				if(data->Branches::q(h)==q_h && cuts::h_sanity(data,envi,h)){
					_hists->Histogram::Fid_Fill(envi,0,theta[h],phi[h],part+1,1,0,_W,data->Branches::p(h));
					_hists->Histogram::DT_Fill(envi,0,part,data->Branches::p(h), data->Branches::sc_r(h), data->Branches::sc_t(h), data->Branches::sc_r(0), data->Branches::sc_t(0),1,0,_W,sector[h]);
					//if(cuts::fid_cut(part+1,data->Branches::p(h),data->Branches::cx(h),data->Branches::cy(h),data->Branches::cz(h))){
					if(cuts::h_fid(data,envi,h,part)){
						fid_h = true; 
						_hists->Histogram::Fid_Fill(envi,0,theta[h],phi[h],part+1,2,0,_W,data->Branches::p(h));
						_hists->Histogram::DT_Fill(envi,0,part,data->Branches::p(h), data->Branches::sc_r(h), data->Branches::sc_t(h), data->Branches::sc_r(0), data->Branches::sc_t(0),2,0,_W,sector[h]);
					}else{
						_hists->Histogram::Fid_Fill(envi,0,theta[h],phi[h],part+1,2,1,_W,data->Branches::p(h));
						_hists->Histogram::DT_Fill(envi,0,part,data->Branches::p(h), data->Branches::sc_r(h), data->Branches::sc_t(h), data->Branches::sc_r(0), data->Branches::sc_t(0),2,1,_W,sector[h]);
					}	
					//if(cuts::delta_t_cut(part, data->Branches::p(h), data->Branches::sc_r(0), data->Branches::sc_r(h), data->Branches::sc_t(0), data->Branches::sc_t(h))){
					if(cuts::h_dt(data,envi,h,part)){
						dt_h = true; 
						_hists->Histogram::Fid_Fill(envi,0,theta[h],phi[h],part+1,3,0,_W,data->Branches::p(h));
						_hists->Histogram::DT_Fill(envi,0,part,data->Branches::p(h), data->Branches::sc_r(h), data->Branches::sc_t(h), data->Branches::sc_r(0), data->Branches::sc_t(0),3,0,_W,sector[h]);
					}else{
						_hists->Histogram::Fid_Fill(envi,0,theta[h],phi[h],part+1,3,1,_W,data->Branches::p(h));
						_hists->Histogram::DT_Fill(envi,0,part,data->Branches::p(h), data->Branches::sc_r(h), data->Branches::sc_t(h), data->Branches::sc_r(0), data->Branches::sc_t(0),3,1,_W,sector[h]);
					}
					if(cuts::hid(data,envi,h,part)){//Contains some additional cuts for pim to try and cut out some electrons
						//std::cout<<std::endl <<"have measured a hadron for species number: " <<part;
						_hists->Histogram::Fid_Fill(envi,0,theta[h],phi[h],part+1,4,0,_W,data->Branches::p(h));
						_hists->Histogram::DT_Fill(envi,0,part,data->Branches::p(h), data->Branches::sc_r(h), data->Branches::sc_t(h), data->Branches::sc_r(0), data->Branches::sc_t(0),4,0,_W,sector[h]);
						switch(part){
							case 0:
								//std::cout<<std::endl <<"have measured a proton ";
								//_prot = physics::Make_4Vector(data->Branches::p(h),data->Branches::cx(h),data->Branches::cy(h),data->Branches::cz(h),mp);
								good_pro++;//Number of identified protons
								//std::cout<<" assign pro idx:" <<h <<std::endl;
								dpro[p_loop] = data->Branches::sc_r(h);
								tpro[p_loop] = data->Branches::sc_t(h);
								cxpro[p_loop] = data->Branches::cx(h);
								cypro[p_loop] = data->Branches::cy(h);
								czpro[p_loop] = data->Branches::cz(h);
								ppro[p_loop] = data->Branches::p(h);
								h_secpro[p_loop] = sector[h];
								check_idx[p_loop] = h; 
								pro_idx[p_loop] = h;
								//std::cout<<"Sector Comparison for Proton: " <<sector[h] <<" " <<h_secpro[p_loop] <<std::endl; 
								p_loop++;
								//std::cout<<"| num good proton = " <<good_pro;
							break;
							case 1:
								//std::cout<<std::endl <<"have measured a pip ";
								//_pip = physics::Make_4Vector(data->Branches::p(h),data->Branches::cx(h),data->Branches::cy(h),data->Branches::cz(h),mpi);
								good_pip++;//number of Identified pip
								//std::cout<<" assign pip idx:" <<h <<std::endl;
								dpip[pip_loop] = data->Branches::sc_r(h);
								tpip[pip_loop] = data->Branches::sc_t(h);
								cxpip[pip_loop] = data->Branches::cx(h);
								cypip[pip_loop] = data->Branches::cy(h);
								czpip[pip_loop] = data->Branches::cz(h);
								ppip[pip_loop] = data->Branches::p(h);
								h_secpip[pip_loop] = sector[h];
								pip_idx[pip_loop]=h;
								check_idx[1] = h; 
								pip_loop++;
								//std::cout<<"| num good pip = " <<good_pip;
							break;
							case 2:
								//std::cout<<std::endl <<"have measured a pim ";
								//_pim = physics::Make_4Vector(data->Branches::p(h),data->Branches::cx(h),data->Branches::cy(h),data->Branches::cz(h),mpi);
								good_pim++;//number of identified pim
								check_idx[2] = h; 
								//std::cout<<" assign pim idx:" <<h <<std::endl;
								dpim[pim_loop] = data->Branches::sc_r(h);
								tpim[pim_loop] = data->Branches::sc_t(h);
								cxpim[pim_loop] = data->Branches::cx(h);
								cypim[pim_loop] = data->Branches::cy(h);
								czpim[pim_loop] = data->Branches::cz(h);
								ppim[pim_loop] = data->Branches::p(h);
								h_secpim[pim_loop] = sector[h];
								pim_idx[pim_loop]=h;
								pim_loop++;
								//std::cout<<"| num good pim = " <<good_pim;
							break;
						}
					}else{
						_hists->Histogram::Fid_Fill(envi,0,theta[h],phi[h],part+1,4,1,_W,data->Branches::p(h));
						_hists->Histogram::DT_Fill(envi,0,part,data->Branches::p(h), data->Branches::sc_r(h), data->Branches::sc_t(h), data->Branches::sc_r(0), data->Branches::sc_t(0),4,1,_W,sector[h]);
					}
					if(data->Branches::id(h) == hid_num){
						_hists->Histogram::Fid_Fill(envi,0,theta[h],phi[h],part+1,5,0,_W,data->Branches::p(h));
						_hists->Histogram::DT_Fill(envi,0,part,data->Branches::p(h), data->Branches::sc_r(h), data->Branches::sc_t(h), data->Branches::sc_r(0), data->Branches::sc_t(0),5,0,_W,sector[h]);
					}else{
						_hists->Histogram::Fid_Fill(envi,0,theta[h],phi[h],part+1,5,1,_W,data->Branches::p(h));
						_hists->Histogram::DT_Fill(envi,0,part,data->Branches::p(h), data->Branches::sc_r(h), data->Branches::sc_t(h), data->Branches::sc_r(0), data->Branches::sc_t(0),5,1,_W,sector[h]);
					}
				}else{
					_hists->Histogram::Fid_Fill(envi,0,theta[h],phi[h],part+1,1,1,_W,data->Branches::p(h));
					_hists->Histogram::DT_Fill(envi,0,part,data->Branches::p(h), data->Branches::sc_r(h), data->Branches::sc_t(h), data->Branches::sc_r(0), data->Branches::sc_t(0),1,1,_W,sector[h]);
				}
			}
		}
	}

	bool double_id = false; 
	int double_idx = -99;
	int num_double_id = 0;

	bool have_prot = false;
	bool have_pip = false;
	bool have_pim = false;

	//This just goes with the first identified particle of that type
	if(check_idx[2]==pim_idx[0]){
		have_pim = true; 
	}
	if(check_idx[1]==pip_idx[0]){
		have_pip = true; 
	}
	if(check_idx[0]==pro_idx[0]){
		have_prot = true; 
	}

	bool ideal_top[4] = {false,false,false,false}; //For isolating events in which only one of each particle was detected for the relevant topology

	//std::cout<<"pro_idx:" <<pro_idx[0];
	//std::cout<<" pip_idx:" <<pip_idx[0];
	//std::cout<<" pim_idx:" <<pim_idx[0] <<std::endl;

	//std::cout<<std::endl <<"p_loop: " <<p_loop << "   pip_loop: " <<pip_loop << "   pim_loop: " <<pim_loop <<std::endl;
	//std::cout<<"good_proton: " <<good_pro << "   good_pip: " <<good_pip <<"   good_pim: " <<good_pim <<std::endl;
	
	//Resolve Multiple Identification Issues as well as double measurement for proton/pip
	TLorentzVector test_pro,test_pip,test_pim; // Four vectors for measuring our fun Missing Masses
	//Proton Missing Topology Look 
	if(good_electron == 1){//Need to have a good electron for the event
		//if(cuts::p_miss(envi)){
		if(envi->was_top(0)){
				//std::cout<<std::endl <<"I am in the p missing piece"; 
			if(good_pim >0 && good_pip >0){
				//std::cout<<std::endl <<"	good pim and pip  ";
				if( good_pip > 1){
					if(good_pim > 1){//Multiple pi+ and pi-
						for(int i = 0; i < good_pip; i++){
							for( int j = 0; j < good_pim; j++){
								test_pip = physics::Make_4Vector(ppip[i],cxpip[i],cypip[i],czpip[i],mpi);
								test_pim = physics::Make_4Vector(ppim[j],cxpim[j],cypim[j],czpim[j],mpi);
								//_hists->Histogram::MM_Fill(envi,0,physics::MM_event(0,0,_elec,test_pip,test_pim),0,0);
								//_hists->Histogram::MM_Fill(envi,0,physics::MM_event(0,1,_elec,test_pip,test_pim),0,1);
								if((physics::MM_event(0,0,_elec,test_pip,test_pim) > (p_center-p_sig))&&(physics::MM_event(0,0,_elec,test_pip,test_pim) < (p_center+p_sig))){
									_pim = test_pim;
									_pip = test_pip;
									top_possible[0]=true;
									d[2] = dpim[j];
									t[2] = tpim[j];
									cx[2] = cxpim[j]; 
									cy[2]  = cypim[j];
									cz[2] = czpim[j];
									_p[2] = ppim[j];
									h_sec[2] = h_secpim[j];
									d[1] = dpip[j];
									t[1] = tpip[j];
									cx[1] = cxpip[j]; 
									cy[1]  = cypip[j];
									cz[1] = czpip[j];
									_p[1] = ppip[j];
									h_sec[1] = h_secpip[j];
								}
							}
						}
					}else{//Multiple pi+, but only 1 good pi-
						_pim = physics::Make_4Vector(ppim[0],cxpim[0],cypim[0],czpim[0],mpi);
						for(int i = 0; i < good_pip; i++){
							test_pip = physics::Make_4Vector(ppip[i],cxpip[i],cypip[i],czpip[i],mpi);
						//	_hists->Histogram::MM_Fill(envi,0,physics::MM_event(0,0,_elec,test_pip,_pim),0,0);
						//	_hists->Histogram::MM_Fill(envi,0,physics::MM_event(0,1,_elec,test_pip,_pim),0,1);
							if((physics::MM_event(0,0,_elec,test_pip,_pim) > (p_center-p_sig))&&(physics::MM_event(0,0,_elec,test_pip,_pim) < (p_center+p_sig))){
								_pip = test_pip;
								top_possible[0]=true;
								d[2] = dpim[0];
								t[2] = tpim[0];
								cx[2] = cxpim[0]; 
								cy[2]  = cypim[0];
								cz[2] = czpim[0];
								_p[2] = ppim[0];
								h_sec[2] = h_secpim[0];
								d[1] = dpip[i];
								t[1] = tpip[i];
								cx[1] = cxpip[i]; 
								cy[1]  = cypip[i];
								cz[1] = czpip[i];
								_p[1] = ppip[i];
								h_sec[1] = h_secpip[i];
							}
						}
					}
				}else{
					if(good_pim > 1){//Multiple pi-, but only 1 good pi+
						_pip = physics::Make_4Vector(ppip[0],cxpip[0],cypip[0],czpip[0],mpi);
						for(int i = 0; i < good_pim; i++){
							test_pim = physics::Make_4Vector(ppim[i],cxpim[i],cypim[i],czpim[i],mpi);
							//_hists->Histogram::MM_Fill(envi,0,physics::MM_event(0,0,_elec,_pip,test_pim),0,0);
							//_hists->Histogram::MM_Fill(envi,0,physics::MM_event(0,1,_elec,_pip,test_pim),0,1);
							if((physics::MM_event(0,0,_elec,_pip,test_pim) > (p_center-p_sig))&&(physics::MM_event(0,0,_elec,_pip,test_pim) < (p_center+p_sig))){
								_pim = test_pim;
								top_possible[0]=true;
								d[2] = dpim[i];
								t[2] = tpim[i];
								cx[2] = cxpim[i]; 
								cy[2]  = cypim[i];
								cz[2] = czpim[i];
								_p[2] = ppim[i];
								h_sec[2] = h_secpim[i];
								d[1] = dpip[0];
								t[1] = tpip[0];
								cx[1] = cxpip[0]; 
								cy[1]  = cypip[0];
								cz[1] = czpip[0];
								_p[1] = ppip[0];
								h_sec[1] = h_secpip[0];
							}
						}
					}else{//We have just one pi+ candidate and one pi- candidate so no need to do MM cuts to legitimize them
						_pim = physics::Make_4Vector(ppim[0],cxpim[0],cypim[0],czpim[0],mpi);
						_pip = physics::Make_4Vector(ppip[0],cxpip[0],cypip[0],czpip[0],mpi);
						ideal_top[0] = true;
						//_hists->Histogram::MM_Fill(envi,0,physics::MM_event(0,0,_elec,_pip,_pim),0,0);
						//_hists->Histogram::MM_Fill(envi,0,physics::MM_event(0,1,_elec,_pip,_pim),0,1);
						top_possible[0]=true;
						d[2] = dpim[0];
						t[2] = tpim[0];
						cx[2] = cxpim[0]; 
						cy[2]  = cypim[0];
						cz[2] = czpim[0];
						_p[2] = ppim[0];
						h_sec[2] = h_secpim[0];
						d[1] = dpip[0];
						t[1] = tpip[0];
						cx[1] = cxpip[0]; 
						cy[1]  = cypip[0];
						cz[1] = czpip[0];
						_p[1] = ppip[0];
						h_sec[1] = h_secpip[0];
					}
				}
			}
		}
		//Pip Missing Topology Look 
		//if(cuts::pip_miss(envi)){	
		if(envi->was_top(1)){
			//std::cout<<std::endl <<"I am in the pip missing piece  "; 
			if(good_pro > 0 && good_pim >0){
				if( good_pro > 1){
					if(good_pim > 1){//Multiple proton and pi-
						for(int i = 0; i < good_pro; i++){
							for( int j = 0; j < good_pim; j++){
								test_pro = physics::Make_4Vector(ppro[i],cxpro[i],cypro[i],czpro[i],mp);
								test_pim = physics::Make_4Vector(ppim[j],cxpim[j],cypim[j],czpim[j],mpi);
							//	_hists->Histogram::MM_Fill(envi,1,physics::MM_event(0,0,_elec,test_pro,test_pim),0,0);
							//	_hists->Histogram::MM_Fill(envi,1,physics::MM_event(0,1,_elec,test_pro,test_pim),0,1);
								if((physics::MM_event(0,0,_elec,test_pro,test_pim) > (pip_center-pip_sig))&&(physics::MM_event(0,0,_elec,test_pro,test_pim) < (pip_center+pip_sig))){
									_pim = test_pim;
									_pro = test_pro;
								//	top_possible[1]=true;
									d[2] = dpim[j];
									t[2] = tpim[j];
									cx[2] = cxpim[j]; 
									cy[2]  = cypim[j];
									cz[2] = czpim[j];
									_p[2] = ppim[j];
									h_sec[2] = h_secpim[j];
									d[0] = dpro[i];
									t[0] = tpro[i];
									cx[0] = cxpro[i]; 
									cy[0]  = cypro[i];
									cz[0] = czpro[i];
									_p[0] = ppro[i];
									h_sec[0] = h_secpro[i];
								}
							}
						}
					}else{//Multiple proton, but only 1 good pi-
						_pim = physics::Make_4Vector(ppim[0],cxpim[0],cypim[0],czpim[0],mpi);
						for(int i = 0; i < good_pro; i++){
							test_pro = physics::Make_4Vector(ppro[i],cxpro[i],cypro[i],czpro[i],mp);
						//	_hists->Histogram::MM_Fill(envi,1,physics::MM_event(0,0,_elec,test_pro,_pim),0,0);
						//	_hists->Histogram::MM_Fill(envi,1,physics::MM_event(0,1,_elec,test_pro,_pim),0,1);
							if((physics::MM_event(0,0,_elec,test_pro,_pim) > (pip_center-pip_sig))&&(physics::MM_event(0,0,_elec,test_pro,_pim) < (pip_center+pip_sig))){
								_pro = test_pro;
							//	top_possible[1]=true;
								d[2] = dpim[0];
								t[2] = tpim[0];
								cx[2] = cxpim[0]; 
								cy[2]  = cypim[0];
								cz[2] = czpim[0];
								_p[2] = ppim[0];
								h_sec[2] = h_secpim[0];
								d[0] = dpro[i];
								t[0] = tpro[i];
								cx[0] = cxpro[i]; 
								cy[0]  = cypro[i];
								cz[0] = czpro[i];
								_p[0] = ppro[i];
								h_sec[0] = h_secpro[i];
							}
						}
					}
				}else{
					if(good_pim > 1){//Multiple pi-, but only 1 good proton
						_pro = physics::Make_4Vector(ppro[0],cxpro[0],cypro[0],czpro[0],mp);
						for(int i = 0; i < good_pim; i++){
							test_pim = physics::Make_4Vector(ppim[i],cxpim[i],cypim[i],czpim[i],mpi);
						//	_hists->Histogram::MM_Fill(envi,1,physics::MM_event(0,0,_elec,_pro,test_pim),0,0);
						//	_hists->Histogram::MM_Fill(envi,1,physics::MM_event(0,1,_elec,_pro,test_pim),0,1);
							if((physics::MM_event(0,0,_elec,_pro,test_pim) > (pip_center-pip_sig))&&(physics::MM_event(0,0,_elec,_pro,test_pim) < (pip_center+pip_sig))){
								_pim = test_pim;
							//	top_possible[1]=true;
								d[2] = dpim[i];
								t[2] = tpim[i];
								cx[2] = cxpim[i]; 
								cy[2]  = cypim[i];
								cz[2] = czpim[i];
								_p[2] = ppim[i];
								h_sec[2] = h_secpim[i];
								d[0] = dpro[0];
								t[0] = tpro[0];
								cx[0] = cxpro[0]; 
								cy[0]  = cypro[0];
								cz[0] = czpro[0];
								_p[0] = ppro[0];
								h_sec[0] = h_secpro[0];
							}
						}
					}else{//We have just one proton candidate and one pi- candidate so no need to do MM cuts to legitimize them
						_pim = physics::Make_4Vector(ppim[0],cxpim[0],cypim[0],czpim[0],mpi);
						_pro = physics::Make_4Vector(ppro[0],cxpro[0],cypro[0],czpro[0],mp);
						ideal_top[1] = true;
						//_hists->Histogram::MM_Fill(envi,1,physics::MM_event(0,0,_elec,_pro,_pim),0,0);
						//_hists->Histogram::MM_Fill(envi,1,physics::MM_event(0,1,_elec,_pro,_pim),0,1);
						top_possible[1]=true;
						d[2] = dpim[0];
						t[2] = tpim[0];
						cx[2] = cxpim[0]; 
						cy[2]  = cypim[0];
						cz[2] = czpim[0];
						_p[2] = ppim[0];
						h_sec[2] = h_secpim[0];
						d[0] = dpro[0];
						t[0] = tpro[0];
						cx[0] = cxpro[0]; 
						cy[0]  = cypro[0];
						cz[0] = czpro[0];
						_p[0] = ppro[0];
						h_sec[0] = h_secpro[0];
					}
				}
			}
		}
		//Pim Missing Topology Look 
		//if(cuts::pim_miss(envi)){
		if(envi->was_top(2)){
			//std::cout<<std::endl <<"I am in the pim missing piece  "; 
			if(good_pro > 0 && good_pip >0){
				if( good_pro > 1){
					if(good_pip > 1){//Multiple proton and pi+
						for(int i = 0; i < good_pro; i++){
							for( int j = 0; j < good_pip; j++){
								test_pro = physics::Make_4Vector(ppro[i],cxpro[i],cypro[i],czpro[i],mp);
								test_pip = physics::Make_4Vector(ppip[j],cxpip[j],cypip[j],czpip[j],mpi);
							//	_hists->Histogram::MM_Fill(envi,2,physics::MM_event(0,0,_elec,test_pro,test_pip),0,0);
							//	_hists->Histogram::MM_Fill(envi,2,physics::MM_event(0,1,_elec,test_pro,test_pip),0,1);
								if((physics::MM_event(0,0,_elec,test_pro,test_pip) > (pim_center-pim_sig))&&(physics::MM_event(0,0,_elec,test_pro,test_pip) < (pim_center+pim_sig)) && (pro_idx[i]!=pip_idx[j])){
									_pip = test_pip;
									_pro = test_pro;
								//	top_possible[2]=true;
									d[1] = dpip[j];
									t[1] = tpip[j];
									cx[1] = cxpip[j]; 
									cy[1]  = cypip[j];
									cz[1] = czpip[j];
									_p[1] = ppip[j];
									h_sec[1] = h_secpip[j];
									d[0] = dpro[i];
									t[0] = tpro[i];
									cx[0] = cxpro[i]; 
									cy[0]  = cypro[i];
									cz[0] = czpro[i];
									_p[0] = ppro[i];
									h_sec[0] = h_secpro[i];
								}
							}
						}
					}else{//Multiple proton, but only 1 good pi+
						_pip = physics::Make_4Vector(ppip[0],cxpip[0],cypip[0],czpip[0],mpi);
						for(int i = 0; i < good_pro; i++){
							test_pro = physics::Make_4Vector(ppro[i],cxpro[i],cypro[i],czpro[i],mp);
						//	_hists->Histogram::MM_Fill(envi,2,physics::MM_event(0,0,_elec,test_pro,_pip),0,0);
						//	_hists->Histogram::MM_Fill(envi,2,physics::MM_event(0,1,_elec,test_pro,_pip),0,1);
							if((physics::MM_event(0,0,_elec,test_pro,_pip) > (pim_center-pim_sig))&&(physics::MM_event(0,0,_elec,test_pro,_pip) < (pim_center+pim_sig))&& (pro_idx[i]!=pip_idx[0])){
								_pro = test_pro;
							//	top_possible[2]=true;
								d[1] = dpip[0];
								t[1] = tpip[0];
								cx[1] = cxpip[0]; 
								cy[1]  = cypip[0];
								cz[1] = czpip[0];
								_p[1] = ppip[0];
								h_sec[1] = h_secpip[0];
								d[0] = dpro[i];
								t[0] = tpro[i];
								cx[0] = cxpro[i]; 
								cy[0]  = cypro[i];
								cz[0] = czpro[i];
								_p[0] = ppro[i];
								h_sec[0] = h_secpro[i];
							}
						}
					}
				}else{
					if(good_pip > 1){//Multiple pi+, but only 1 good proton
						_pro = physics::Make_4Vector(ppro[0],cxpro[0],cypro[0],czpro[0],mp);
						for(int i = 0; i < good_pip; i++){
							test_pip = physics::Make_4Vector(ppip[i],cxpip[i],cypip[i],czpip[i],mpi);
						//	_hists->Histogram::MM_Fill(envi,2,physics::MM_event(0,0,_elec,_pro,test_pip),0,0);
						//	_hists->Histogram::MM_Fill(envi,2,physics::MM_event(0,1,_elec,_pro,test_pip),0,1);
							if((physics::MM_event(0,0,_elec,_pro,test_pip) > (pim_center-pim_sig))&&(physics::MM_event(0,0,_elec,_pro,test_pip) < (pim_center+pim_sig))&& (pro_idx[0]!=pip_idx[i])){
								_pip = test_pip;
							//	top_possible[2]=true;
								d[1] = dpip[i];
								t[1] = tpip[i];
								cx[1] = cxpip[i]; 
								cy[1]  = cypip[i];
								cz[1] = czpip[i];
								_p[1] = ppip[i];
								h_sec[1] = h_secpip[i];
								d[0] = dpro[0];
								t[0] = tpro[0];
								cx[0] = cxpro[0]; 
								cy[0]  = cypro[0];
								cz[0] = czpro[0];
								_p[0] = ppro[0];
								h_sec[0] = h_secpro[0];
							}
						}
					}else{//We have just one proton candidate and one pi+ candidate so no need to do MM cuts to legitimize them
						if(pro_idx[0] != pip_idx[0]){
							_pip = physics::Make_4Vector(ppip[0],cxpip[0],cypip[0],czpip[0],mpi);
							_pro = physics::Make_4Vector(ppro[0],cxpro[0],cypro[0],czpro[0],mp);
							ideal_top[2] = true;
						//	_hists->Histogram::MM_Fill(envi,2,physics::MM_event(0,0,_elec,_pro,_pip),0,0);
						//	_hists->Histogram::MM_Fill(envi,2,physics::MM_event(0,1,_elec,_pro,_pip),0,1);
							top_possible[2]=true;
							d[1] = dpip[0];
							t[1] = tpip[0];
							cx[1] = cxpip[0]; 
							cy[1]  = cypip[0];
							cz[1] = czpip[0];
							_p[1] = ppip[0];
							h_sec[1] = h_secpip[0];
							d[0] = dpro[0];
							t[0] = tpro[0];
							cx[0] = cxpro[0]; 
							cy[0]  = cypro[0];
							cz[0] = czpro[0];
							_p[0] = ppro[0];
							h_sec[0] = h_secpro[0];
						}else{//They've been identified as the same particle, but we've also not measured a pi- therefore this event just cannot count for anything
							top_possible[2]=false;
						}
					}
				}
			}
		}
		//Zero Missing Topology Look 
		//if(cuts::z_miss(envi)){
		if(envi->was_top(3)){
			//std::cout<<std::endl <<"I am in the z missing piece  "; 
			if(good_pro >0 && good_pip > 0 && good_pim >0){
				if(good_pro>1){//More than 1 IDed Proton
					if(good_pip > 1){//More than 1 IDed Proton and Pi+
						if(good_pim>1){//More than 1 IDed Proton, Pi+, and Pi-
							for(int i = 0; i<p_loop; i++){
								for(int j = 0; j<pip_loop; j++){
									for(int k = 0; k<pim_loop; k++){
										test_pro = physics::Make_4Vector(ppro[i],cxpro[i],cypro[i],czpro[i],mp);
										test_pip = physics::Make_4Vector(ppip[j],cxpip[j],cypip[j],czpip[j],mpi);
										test_pim = physics::Make_4Vector(ppim[k],cxpim[k],cypim[k],czpim[k],mpi);
									//	_hists->Histogram::MM_Fill(envi,3,physics::MM_event(0,0,_elec,test_pro,test_pip,test_pim),0,0);
									//	_hists->Histogram::MM_Fill(envi,3,physics::MM_event(0,1,_elec,test_pro,test_pip,test_pim),0,1);
										if((physics::MM_event(0,1,_elec,test_pro,test_pip,test_pim) > (MM_zero_center2-MM_zero_sigma2))&&(physics::MM_event(0,1,_elec,test_pro,test_pip,test_pim) < (MM_zero_center2+MM_zero_sigma2))){
											if(pro_idx[i]!= pip_idx[j]){
												_pro = test_pro;
												_pip = test_pip;
												_pim = test_pim;
											//	top_possible[3] = true;
												d[1] = dpip[j];
												t[1] = tpip[j];
												cx[1] = cxpip[j]; 
												cy[1]  = cypip[j];
												cz[1] = czpip[j];
												_p[1] = ppip[j];
												h_sec[1] = h_secpip[j];
												d[0] = dpro[i];
												t[0] = tpro[i];
												cx[0] = cxpro[i]; 
												cy[0]  = cypro[i];
												cz[0] = czpro[i];
												_p[0] = ppro[i];
												h_sec[0] = h_secpro[i];
												d[2] = dpim[k];
												t[2] = tpim[k];
												cx[2] = cxpim[k]; 
												cy[2]  = cypim[k];
												cz[2] = czpim[k];
												_p[2] = ppim[k];
												h_sec[2] = h_secpim[k];
											}
										}
									}
								}
							}
						}else{//More than 1 IDed Proton and Pi+, but 1 good Pi-
							for(int i = 0; i<p_loop; i++){
								for(int j = 0; j<pip_loop; j++){
									test_pro = physics::Make_4Vector(ppro[i],cxpro[i],cypro[i],czpro[i],mp);
									test_pip = physics::Make_4Vector(ppip[j],cxpip[j],cypip[j],czpip[j],mpi);
									_pim = physics::Make_4Vector(ppim[0],cxpim[0],cypim[0],czpim[0],mpi);
								//	_hists->Histogram::MM_Fill(envi,3,physics::MM_event(0,0,_elec,test_pro,test_pip,_pim),0,0);
								//	_hists->Histogram::MM_Fill(envi,3,physics::MM_event(0,1,_elec,test_pro,test_pip,_pim),0,1);
									if((physics::MM_event(0,1,_elec,test_pro,test_pip,_pim) > (MM_zero_center2-MM_zero_sigma2))&&(physics::MM_event(0,1,_elec,test_pro,test_pip,_pim) < (MM_zero_center2+MM_zero_sigma2))){
										if(pro_idx[i]!= pip_idx[j]){
											_pro = test_pro;
											_pip = test_pip;
										//	top_possible[3] = true;
											d[1] = dpip[j];
											t[1] = tpip[j];
											cx[1] = cxpip[j]; 
											cy[1]  = cypip[j];
											cz[1] = czpip[j];
											_p[1] = ppip[j];
											h_sec[1] = h_secpip[j];
											d[0] = dpro[i];
											t[0] = tpro[i];
											cx[0] = cxpro[i]; 
											cy[0]  = cypro[i];
											cz[0] = czpro[i];
											_p[0] = ppro[i];
											h_sec[0] = h_secpro[i];
											d[2] = dpim[0];
											t[2] = tpim[0];
											cx[2] = cxpim[0]; 
											cy[2]  = cypim[0];
											cz[2] = czpim[0];
											_p[2] = ppim[0];
											h_sec[2] = h_secpim[0];
										}
									}
								}
							}
						}
					}else{
						if(good_pim>1){//More than 1 IDed Proton and Pi-, but 1 good Pi+
							for(int i = 0; i<p_loop; i++){
								for(int j = 0; j<pim_loop; j++){
									test_pro = physics::Make_4Vector(ppro[i],cxpro[i],cypro[i],czpro[i],mp);
									test_pim = physics::Make_4Vector(ppim[j],cxpim[j],cypim[j],czpim[j],mpi);
									_pip = physics::Make_4Vector(ppip[0],cxpip[0],cypip[0],czpip[0],mpi);
								//	_hists->Histogram::MM_Fill(envi,3,physics::MM_event(0,0,_elec,test_pro,_pip,test_pim),0,0);
								//	_hists->Histogram::MM_Fill(envi,3,physics::MM_event(0,1,_elec,test_pro,_pip,test_pim),0,1);
									if((physics::MM_event(0,1,_elec,test_pro,_pip,test_pim) > (MM_zero_center2-MM_zero_sigma2))&&(physics::MM_event(0,1,_elec,test_pro,_pip,test_pim) < (MM_zero_center2+MM_zero_sigma2))){
										if(pro_idx[i]!= pip_idx[0]){
											_pro = test_pro;
											_pim = test_pim;
										//	top_possible[3] = true;
											d[1] = dpip[0];
											t[1] = tpip[0];
											cx[1] = cxpip[0]; 
											cy[1]  = cypip[0];
											cz[1] = czpip[0];
											_p[1] = ppip[0];
											h_sec[1] = h_secpip[0];
											d[0] = dpro[i];
											t[0] = tpro[i];
											cx[0] = cxpro[i]; 
											cy[0]  = cypro[i];
											cz[0] = czpro[i];
											_p[0] = ppro[i];
											h_sec[0] = h_secpro[i];
											d[2] = dpim[j];
											t[2] = tpim[j];
											cx[2] = cxpim[j]; 
											cy[2]  = cypim[j];
											cz[2] = czpim[j];
											_p[2] = ppim[j];
											h_sec[2] = h_secpim[j];
										}
									}
								}
							}
						}else{//More than 1 IDed Proton, but 1 good pi+ and pi-
							for(int i = 0; i<p_loop; i++){
								test_pro = physics::Make_4Vector(ppro[i],cxpro[i],cypro[i],czpro[i],mp);
								_pim = physics::Make_4Vector(ppim[0],cxpim[0],cypim[0],czpim[0],mpi);
								_pip = physics::Make_4Vector(ppip[0],cxpip[0],cypip[0],czpip[0],mpi);
							//	_hists->Histogram::MM_Fill(envi,3,physics::MM_event(0,0,_elec,test_pro,_pip,_pim),0,0);
							//	_hists->Histogram::MM_Fill(envi,3,physics::MM_event(0,1,_elec,test_pro,_pip,_pim),0,1);
								if((physics::MM_event(0,1,_elec,test_pro,_pip,_pim) > (MM_zero_center2-MM_zero_sigma2))&&(physics::MM_event(0,1,_elec,test_pro,_pip,_pim) < (MM_zero_center2+MM_zero_sigma2))){
									if(pro_idx[i]!= pip_idx[0]){
										_pro = test_pro;
									//	top_possible[3] = true;
										d[1] = dpip[0];
										t[1] = tpip[0];
										cx[1] = cxpip[0]; 
										cy[1]  = cypip[0];
										cz[1] = czpip[0];
										_p[1] = ppip[0];
										h_sec[1] = h_secpip[0];
										d[0] = dpro[i];
										t[0] = tpro[i];
										cx[0] = cxpro[i]; 
										cy[0]  = cypro[i];
										cz[0] = czpro[i];
										_p[0] = ppro[i];
										h_sec[0] = h_secpro[i];
										d[2] = dpim[0];
										t[2] = tpim[0];
										cx[2] = cxpim[0]; 
										cy[2]  = cypim[0];
										cz[2] = czpim[0];
										_p[2] = ppim[0];
										h_sec[2] = h_secpim[0];
									}
								}
							}
						}
					}
				}else{//1 good proton
					if(good_pip > 1){//1 good proton, but multiple pi+
						if(good_pim>1){//1 good proton, but multiple pi+ and pi-
							for(int i = 0; i<pip_loop; i++){
								for(int j = 0; j<pim_loop; j++){
									test_pip = physics::Make_4Vector(ppip[i],cxpip[i],cypip[i],czpip[i],mpi);
									test_pim = physics::Make_4Vector(ppim[j],cxpim[j],cypim[j],czpim[j],mpi);
									_pro = physics::Make_4Vector(ppro[0],cxpro[0],cypro[0],czpro[0],mp);
								//	_hists->Histogram::MM_Fill(envi,3,physics::MM_event(0,0,_elec,_pro,test_pip,test_pim),0,0);
								//	_hists->Histogram::MM_Fill(envi,3,physics::MM_event(0,1,_elec,_pro,test_pip,test_pim),0,1);
									if((physics::MM_event(0,1,_elec,_pro,test_pip,test_pim) > (MM_zero_center2-MM_zero_sigma2))&&(physics::MM_event(0,1,_elec,_pro,test_pip,test_pim) < (MM_zero_center2+MM_zero_sigma2))){
										if(pro_idx[0]!= pip_idx[i]){
											_pip = test_pip;
											_pim = test_pim;
										//	top_possible[3] = true;
											d[1] = dpip[i];
											t[1] = tpip[i];
											cx[1] = cxpip[i]; 
											cy[1]  = cypip[i];
											cz[1] = czpip[i];
											_p[1] = ppip[i];
											h_sec[1] = h_secpip[i];
											d[0] = dpro[0];
											t[0] = tpro[0];
											cx[0] = cxpro[0]; 
											cy[0]  = cypro[0];
											cz[0] = czpro[0];
											_p[0] = ppro[0];
											h_sec[0] = h_secpro[0];
											d[2] = dpim[j];
											t[2] = tpim[j];
											cx[2] = cxpim[j]; 
											cy[2]  = cypim[j];
											cz[2] = czpim[j];
											_p[2] = ppim[j];
											h_sec[2] = h_secpim[j];
										}
									}
								}
							}
						}else{//1 good proton and pi-, but multiple pi+
							for(int i = 0; i<pip_loop; i++){
								test_pip = physics::Make_4Vector(ppip[i],cxpip[i],cypip[i],czpip[i],mpi);
								_pim = physics::Make_4Vector(ppim[0],cxpim[0],cypim[0],czpim[0],mpi);
								_pro = physics::Make_4Vector(ppro[0],cxpro[0],cypro[0],czpro[0],mp);
							//	_hists->Histogram::MM_Fill(envi,3,physics::MM_event(0,0,_elec,_pro,test_pip,_pim),0,0);
							//	_hists->Histogram::MM_Fill(envi,3,physics::MM_event(0,1,_elec,_pro,test_pip,_pim),0,1);
								if((physics::MM_event(0,1,_elec,_pro,test_pip,_pim) > (MM_zero_center2-MM_zero_sigma2))&&(physics::MM_event(0,1,_elec,_pro,test_pip,_pim) < (MM_zero_center2+MM_zero_sigma2))){
									if(pro_idx[i]!= pip_idx[0]){
										_pip = test_pip;
									//	top_possible[3] = true;
										d[1] = dpip[i];
										t[1] = tpip[i];
										cx[1] = cxpip[i]; 
										cy[1]  = cypip[i];
										cz[1] = czpip[i];
										_p[1] = ppip[i];
										h_sec[1] = h_secpip[i];
										d[0] = dpro[0];
										t[0] = tpro[0];
										cx[0] = cxpro[0]; 
										cy[0]  = cypro[0];
										cz[0] = czpro[0];
										_p[0] = ppro[0];
										h_sec[0] = h_secpro[0];
										d[2] = dpim[0];
										t[2] = tpim[0];
										cx[2] = cxpim[0]; 
										cy[2]  = cypim[0];
										cz[2] = czpim[0];
										_p[2] = ppim[0];
										h_sec[2] = h_secpim[0];
									}
								}
							}
						}
					}else{//1 good proton and pi+
						if(good_pim>1){//1 good proton and pi+, but multiple pi-
							for(int i = 0; i<pip_loop; i++){
								test_pim = physics::Make_4Vector(ppim[i],cxpim[i],cypim[i],czpim[i],mpi);
								_pip = physics::Make_4Vector(ppip[0],cxpip[0],cypip[0],czpip[0],mpi);
								_pro = physics::Make_4Vector(ppro[0],cxpro[0],cypro[0],czpro[0],mp);
							//	_hists->Histogram::MM_Fill(envi,3,physics::MM_event(0,0,_elec,_pro,_pip,test_pim),0,0);
							//	_hists->Histogram::MM_Fill(envi,3,physics::MM_event(0,1,_elec,_pro,_pip,test_pim),0,1);
								if((physics::MM_event(0,1,_elec,_pro,_pip,test_pim) > (MM_zero_center2-MM_zero_sigma2))&&(physics::MM_event(0,1,_elec,_pro,_pip,test_pim) < (MM_zero_center2+MM_zero_sigma2))){
									if(pro_idx[0]!= pip_idx[0]){
										_pim = test_pim;
										top_possible[3] = true;
										d[1] = dpip[i];
										t[1] = tpip[i];
										cx[1] = cxpip[i]; 
										cy[1]  = cypip[i];
										cz[1] = czpip[i];
										_p[1] = ppip[i];
										h_sec[1] = h_secpip[i];
										d[0] = dpro[0];
										t[0] = tpro[0];
										cx[0] = cxpro[0]; 
										cy[0]  = cypro[0];
										cz[0] = czpro[0];
										_p[0] = ppro[0];
										h_sec[0] = h_secpro[0];
										d[2] = dpim[0];
										t[2] = tpim[0];
										cx[2] = cxpim[0]; 
										cy[2]  = cypim[0];
										cz[2] = czpim[0];
										_p[2] = ppim[0];
										h_sec[2] = h_secpim[0];
									}
								}
							}
						}else{//1 good proton and pi+ and pi-
							if(pro_idx[0]!=pip_idx[0]){
								ideal_top[3] = true;
								top_possible[3]=true;
								d[1] = dpip[0];
								t[1] = tpip[0];
								cx[1] = cxpip[0]; 
								cy[1]  = cypip[0];
								cz[1] = czpip[0];
								_p[1] = ppip[0];
								h_sec[1] = h_secpip[0];
								d[0] = dpro[0];
								t[0] = tpro[0];
								cx[0] = cxpro[0]; 
								cy[0]  = cypro[0];
								cz[0] = czpro[0];
								_p[0] = ppro[0];
								h_sec[0] = h_secpro[0];
								d[2] = dpim[0];
								t[2] = tpim[0];
								cx[2] = cxpim[0]; 
								cy[2]  = cypim[0];
								cz[2] = czpim[0];
								_p[2] = ppim[0];
								h_sec[2] = h_secpim[0];
								_pip = physics::Make_4Vector(ppip[0],cxpip[0],cypip[0],czpip[0],mpi);
								_pro = physics::Make_4Vector(ppro[0],cxpro[0],cypro[0],czpro[0],mp);
								_pim = physics::Make_4Vector(ppim[0],cxpim[0],cypim[0],czpim[0],mpi);
							//	_hists->Histogram::MM_Fill(envi,3,physics::MM_event(0,0,_elec,_pro,_pip,_pim),0,0);
							//	_hists->Histogram::MM_Fill(envi,3,physics::MM_event(0,1,_elec,_pro,_pip,_pim),0,1);
							}else{
								top_possible[3]=false;
							}
						}
					}
				}
			}
		}
	}//Oh this closes it off at the moment	

	//Check to see that four vectors have been properly assigned
	for(int w = 0; w < 4; w++){
		if(top_possible[w]){
			switch(w){
				case 0: 
				_assigned_4vecs = (physics::Check_4Vec(_elec)&&physics::Check_4Vec(_pip)&&physics::Check_4Vec(_pim));
				break;
				case 1: 
				_assigned_4vecs = (physics::Check_4Vec(_elec)&&physics::Check_4Vec(_pro)&&physics::Check_4Vec(_pim));
				break;
				case 2: 
				_assigned_4vecs = (physics::Check_4Vec(_elec)&&physics::Check_4Vec(_pip)&&physics::Check_4Vec(_pro));
				break;
				case 3: 
				_assigned_4vecs = (physics::Check_4Vec(_elec)&&physics::Check_4Vec(_pip)&&physics::Check_4Vec(_pim)&&physics::Check_4Vec(_pro));
				break;
			}
			if(_assigned_4vecs == false){
				std::cout<<std::endl<<"Improperly Filled Four Vectors for " <<topologies[w+1] <<" elec: " <<physics::Check_4Vec(_elec) <<" pro: " <<physics::Check_4Vec(_pro) <<" pip: " <<physics::Check_4Vec(_pip) <<" pim: " <<physics::Check_4Vec(_pim)<<std::endl;
			}
		}
	}
	
	bool topo[4]={false,false,false,false};
	int part;

	//Event_Class Selection
	if(good_electron == 1){//Good Electron
		//P Missing
		if(top_possible[0]){//Good Pip and Pim
			//physics::Print_4Vec(_elec);
			//physics::Print_4Vec(_pip);
			//physics::Print_4Vec(_pim);
			MM_p = physics::MM_event(0,0,_elec,_pip,_pim);
			MM_p2 = physics::MM_event(0,1,_elec,_pip,_pim);
			//std::cout<<"MM squared for proton missing topology: " <<MM_p2 <<std::endl;
			_hists->Histogram::MM_Fill(envi,0,MM_p,0,0,ideal_top[0]);//Added this up in the event particle determination step
			_hists->Histogram::MM_Fill(envi,0,MM_p2,0,1,ideal_top[0]);
			
			//_hists->Histogram::MM_Fill(envi,4,MM_p,0,0);
			//_hists->Histogram::MM_Fill(envi,4,MM_p2,0,1);
			if((MM_p > (p_center-p_sig))&&(MM_p < (p_center+p_sig))){//Missing Mass Cut on Proton Mass
				_hists->Histogram::MM_Fill(envi,0,MM_p,1,0,ideal_top[0]);
				_hists->Histogram::MM_Fill(envi,0,MM_p2,1,1,ideal_top[0]);
				
				//_hists->Histogram::MM_Fill(envi,4,MM_p,1,0);
				//_hists->Histogram::MM_Fill(envi,4,MM_p2,1,1);
				topo[0]=true;
				//std::cout<<std::endl <<"PRO Missing topology passed";
				/*if(!top_possible[3]){
					_pro = _beam + _target - _elec - _pip - _pim; 
				}*/
				_hists->Histogram::WQ2_Fill(envi,1,10,_W,_Q2);
				_hists->Histogram::Fid_Fill(envi,1,physics::get_theta(data->Branches::cz(0)),physics::get_phi(data->Branches::cx(0),data->Branches::cy(0)),0,10,0,_W,data->Branches::p(0));
				_hists->Histogram::SF_Fill(envi,1,data->Branches::p(0),data->Branches::etot(0),10,0,_W,sector[0]);
				_hists->Histogram::CC_Fill(envi,1,data->Branches::cc_sect(0),data->Branches::cc_segm(0),data->Branches::nphe(0),10,0);
				for(int par = 1; par<3; par++){
					_hists->Histogram::Fid_Fill(envi,1,physics::get_theta(cz[par]),physics::get_phi(cx[par],cy[par]),par+1,6,0,_W,_p[par]);
					_hists->Histogram::DT_Fill(envi,1,par,_p[par], d[par], t[par], data->Branches::sc_r(0), data->Branches::sc_t(0),6,0,_W,physics::get_sector(physics::get_phi(cx[par],cy[par])));
				}
			}else{
				_hists->Histogram::MM_Fill(envi,0,MM_p,2,0,ideal_top[0]);
				_hists->Histogram::MM_Fill(envi,0,MM_p2,2,1,ideal_top[0]);
				
				//_hists->Histogram::MM_Fill(envi,4,MM_p,2,0);
				//_hists->Histogram::MM_Fill(envi,4,MM_p2,2,1);
				//_hists->Histogram::WQ2_Fill(envi,1,10,_W,_Q2);
				_hists->Histogram::Fid_Fill(envi,1,physics::get_theta(data->Branches::cz(0)),physics::get_phi(data->Branches::cx(0),data->Branches::cy(0)),0,10,1,_W,data->Branches::p(0));
				_hists->Histogram::SF_Fill(envi,1,data->Branches::p(0),data->Branches::etot(0),10,1,_W,sector[0]);
				_hists->Histogram::CC_Fill(envi,1,data->Branches::cc_sect(0),data->Branches::cc_segm(0),data->Branches::nphe(0),10,1);
				for(int par = 1; par<3; par++){
					_hists->Histogram::Fid_Fill(envi,1,physics::get_theta(cz[par]),physics::get_phi(cx[par],cy[par]),par+1,6,1,_W,_p[par]);
					_hists->Histogram::DT_Fill(envi,1,par,_p[par], d[par], t[par], data->Branches::sc_r(0), data->Branches::sc_t(0),6,1,_W,physics::get_sector(physics::get_phi(cx[par],cy[par])));
				}
			}
		}
		//Pip Missing
		if(top_possible[1]){//Good proton and Pim
			MM_pip = physics::MM_event(0,0,_elec,_pro,_pim);
			MM_pip2 = physics::MM_event(0,1,_elec,_pro,_pim);
			_hists->Histogram::MM_Fill(envi,1,MM_pip,0,0,ideal_top[1]);
			_hists->Histogram::MM_Fill(envi,1,MM_pip2,0,1,ideal_top[1]);
			
			//_hists->Histogram::MM_Fill(envi,4,MM_pip,0,0);
			//_hists->Histogram::MM_Fill(envi,4,MM_pip2,0,1);
			if((MM_pip > (pip_center-pip_sig))&&(MM_pip < (pip_center+pip_sig))){//Missing Mass Cut on proton Mass
				_hists->Histogram::MM_Fill(envi,1,MM_pip,1,0,ideal_top[1]);
				_hists->Histogram::MM_Fill(envi,1,MM_pip2,1,1,ideal_top[1]);
				
				//_hists->Histogram::MM_Fill(envi,4,MM_pip,1,0);
				//_hists->Histogram::MM_Fill(envi,4,MM_pip2,1,1);
				topo[1]=true;
				//std::cout<<std::endl <<"PIP Missing topology passed";
				/*if(!top_possible[3]){
					_pip = _beam + _target - _elec - _pro - _pim; 
				}*/
				//Fill Event Selection 
				//Electron
				_hists->Histogram::WQ2_Fill(envi,2,10,_W,_Q2);
				_hists->Histogram::Fid_Fill(envi,2,physics::get_theta(data->Branches::cz(0)),physics::get_phi(data->Branches::cx(0),data->Branches::cy(0)),0,10,0,_W,data->Branches::p(0));
				_hists->Histogram::SF_Fill(envi,2,data->Branches::p(0),data->Branches::etot(0),10,0,_W,sector[0]);
				_hists->Histogram::CC_Fill(envi,2,data->Branches::cc_sect(0),data->Branches::cc_segm(0),data->Branches::nphe(0),10,0);
				//Hadron
				for(int i = 0; i<2; i++){
					switch(i){
						case 0:
						part = 0; 
						break;
						case 1:
						part = 2; 
						break;
					}
					_hists->Histogram::Fid_Fill(envi,2,physics::get_theta(cz[part]),physics::get_phi(cx[part],cy[part]),part+1,6,0,_W,_p[part]);
					_hists->Histogram::DT_Fill(envi,2,part,_p[part], d[part], t[part], data->Branches::sc_r(0), data->Branches::sc_t(0),6,0,_W,physics::get_sector(physics::get_phi(cx[part],cy[part])));
				}
			}else{
				_hists->Histogram::MM_Fill(envi,1,MM_pip,2,0,ideal_top[1]);
				_hists->Histogram::MM_Fill(envi,1,MM_pip2,2,1,ideal_top[1]);
				
				//_hists->Histogram::MM_Fill(envi,4,MM_pip,2,0);
				//_hists->Histogram::MM_Fill(envi,4,MM_pip2,2,1);
				//_hists->Histogram::WQ2_Fill(envi,2,10,_W,_Q2);
				_hists->Histogram::Fid_Fill(envi,2,physics::get_theta(data->Branches::cz(0)),physics::get_phi(data->Branches::cx(0),data->Branches::cy(0)),0,10,1,_W,data->Branches::p(0));
				_hists->Histogram::SF_Fill(envi,2,data->Branches::p(0),data->Branches::etot(0),10,1,_W,sector[0]);
				_hists->Histogram::CC_Fill(envi,2,data->Branches::cc_sect(0),data->Branches::cc_segm(0),data->Branches::nphe(0),10,1);
				for(int i = 0; i<2; i++){
					switch(i){
						case 0:
						part = 0; 
						break;
						case 1:
						part = 2; 
						break;
					}
					_hists->Histogram::Fid_Fill(envi,2,physics::get_theta(cz[part]),physics::get_phi(cx[part],cy[part]),part+1,6,1,_W,_p[part]);
					_hists->Histogram::DT_Fill(envi,2,part,_p[part], d[part], t[part], data->Branches::sc_r(0), data->Branches::sc_t(0),6,1,_W,physics::get_sector(physics::get_phi(cx[part],cy[part])));
				}
			}
		}
		//Pim Missing //Need to introduce some proton/Pion discrimination here, but will leave for now
		if(top_possible[2]){//Good Proton and Pip
			MM_pim = physics::MM_event(0,0,_elec,_pro,_pip);
			MM_pim2 = physics::MM_event(0,1,_elec,_pro,_pip);
			_hists->Histogram::MM_Fill(envi,2,MM_pim,0,0,ideal_top[2]);
			_hists->Histogram::MM_Fill(envi,2,MM_pim2,0,1,ideal_top[2]);
			
			//_hists->Histogram::MM_Fill(envi,4,MM_pim,0,0);
			//_hists->Histogram::MM_Fill(envi,4,MM_pim2,0,1);
			if((MM_pim > (pim_center-pim_sig))&&(MM_pim < (pim_center+pim_sig))){//Missing Mass Cut on Proton Mass
				_hists->Histogram::MM_Fill(envi,2,MM_pim,1,0,ideal_top[2]);
				_hists->Histogram::MM_Fill(envi,2,MM_pim2,1,1,ideal_top[2]);
				
				//_hists->Histogram::MM_Fill(envi,4,MM_pim,1,0);
				//_hists->Histogram::MM_Fill(envi,4,MM_pim2,1,1);
				topo[2]=true;
			//	std::cout<<std::endl <<"PIM Missing topology passed";
				_hists->Histogram::WQ2_Fill(envi,3,10,_W,_Q2);
				_hists->Histogram::Fid_Fill(envi,3,physics::get_theta(data->Branches::cz(0)),physics::get_phi(data->Branches::cx(0),data->Branches::cy(0)),0,10,0,_W,data->Branches::p(0));
				_hists->Histogram::SF_Fill(envi,3,data->Branches::p(0),data->Branches::etot(0),10,0,_W,sector[0]);
				_hists->Histogram::CC_Fill(envi,3,data->Branches::cc_sect(0),data->Branches::cc_segm(0),data->Branches::nphe(0),10,0);
				/*if(!top_possible[3]){
					_pim = _beam + _target - _elec - _pro - _pip; 
				}*/
				for(int par = 0; par<2; par++){
					_hists->Histogram::Fid_Fill(envi,3,physics::get_theta(cz[par]),physics::get_phi(cx[par],cy[par]),par+1,6,0,_W,_p[par]);
					_hists->Histogram::DT_Fill(envi,3,par,_p[par], d[par], t[par], data->Branches::sc_r(0), data->Branches::sc_t(0),6,0,_W,physics::get_sector(physics::get_phi(cx[par],cy[par])));
				}
			}else{
				_hists->Histogram::MM_Fill(envi,2,MM_pim,2,0,ideal_top[2]);
				_hists->Histogram::MM_Fill(envi,2,MM_pim2,2,1,ideal_top[2]);
				
				//_hists->Histogram::MM_Fill(envi,4,MM_pim,2,0);
				//_hists->Histogram::MM_Fill(envi,4,MM_pim2,2,1);
				//_hists->Histogram::WQ2_Fill(envi,3,10,_W,_Q2);
				_hists->Histogram::Fid_Fill(envi,3,physics::get_theta(data->Branches::cz(0)),physics::get_phi(data->Branches::cx(0),data->Branches::cy(0)),0,10,1,_W,data->Branches::p(0));
				_hists->Histogram::SF_Fill(envi,3,data->Branches::p(0),data->Branches::etot(0),10,1,_W,sector[0]);
				_hists->Histogram::CC_Fill(envi,3,data->Branches::cc_sect(0),data->Branches::cc_segm(0),data->Branches::nphe(0),10,1);
				for(int par = 0; par<2; par++){
					_hists->Histogram::Fid_Fill(envi,3,physics::get_theta(cz[par]),physics::get_phi(cx[par],cy[par]),par+1,6,1,_W,_p[par]);
					_hists->Histogram::DT_Fill(envi,3,par,_p[par], d[par], t[par], data->Branches::sc_r(0), data->Branches::sc_t(0),6,1,_W,physics::get_sector(physics::get_phi(cx[par],cy[par])));
				}
			}
		}
		//None Missing
		if(top_possible[3]){//Good Proton Pip and Pim 
			MM_z = physics::MM_event(0,0,_elec,_pro,_pip,_pim);
			MM_z2 = physics::MM_event(0,1,_elec,_pro,_pip,_pim);
			_hists->Histogram::MM_Fill(envi,3,MM_z,0,0,ideal_top[3]);
			_hists->Histogram::MM_Fill(envi,3,MM_z2,0,1,ideal_top[3]);
			//_hists->Histogram::MM_Fill(envi,4,MM_z,0,0);
			//_hists->Histogram::MM_Fill(envi,4,MM_z2,0,1);
			if((MM_z2 > (MM_zero_center2-MM_zero_sigma2))&&(MM_z2 < (MM_zero_center2+MM_zero_sigma2))){//Missing Mass Cut on Proton Mass
				_hists->Histogram::MM_Fill(envi,3,MM_p,1,0,ideal_top[3]);//Added this up in the event particle determination step
				_hists->Histogram::MM_Fill(envi,3,MM_p2,1,1,ideal_top[3]);
				//_hists->Histogram::MM_Fill(envi,4,MM_z,1,0);
				//_hists->Histogram::MM_Fill(envi,4,MM_z2,1,1);
				topo[3]=true;
				//std::cout<<std::endl <<"ZERO Missing topology passed";
				_hists->Histogram::WQ2_Fill(envi,4,10,_W,_Q2);
				_hists->Histogram::Fid_Fill(envi,4,physics::get_theta(data->Branches::cz(0)),physics::get_phi(data->Branches::cx(0),data->Branches::cy(0)),0,10,0,_W,data->Branches::p(0));
				_hists->Histogram::SF_Fill(envi,4,data->Branches::p(0),data->Branches::etot(0),10,0,_W,sector[0]);
				_hists->Histogram::CC_Fill(envi,4,data->Branches::cc_sect(0),data->Branches::cc_segm(0),data->Branches::nphe(0),10,0);
				for(int par = 0; par<3; par++){
					_hists->Histogram::Fid_Fill(envi,4,physics::get_theta(cz[par]),physics::get_phi(cx[par],cy[par]),par+1,6,0,_W,_p[par]);
					_hists->Histogram::DT_Fill(envi,4,par,_p[par], d[par], t[par], data->Branches::sc_r(0), data->Branches::sc_t(0),6,0,_W,physics::get_sector(physics::get_phi(cx[par],cy[par])));
				}
			}else{
				_hists->Histogram::MM_Fill(envi,3,MM_p,2,0,ideal_top[3]);//Added this up in the event particle determination step
				_hists->Histogram::MM_Fill(envi,3,MM_p2,2,1,ideal_top[3]);
				
				//_hists->Histogram::MM_Fill(envi,4,MM_z,2,0);
				//_hists->Histogram::MM_Fill(envi,4,MM_z2,2,1);
				//_hists->Histogram::WQ2_Fill(envi,4,10,_W,_Q2);
				_hists->Histogram::Fid_Fill(envi,4,physics::get_theta(data->Branches::cz(0)),physics::get_phi(data->Branches::cx(0),data->Branches::cy(0)),0,10,1,_W,data->Branches::p(0));
				_hists->Histogram::SF_Fill(envi,4,data->Branches::p(0),data->Branches::etot(0),10,1,_W,sector[0]);
				_hists->Histogram::CC_Fill(envi,4,data->Branches::cc_sect(0),data->Branches::cc_segm(0),data->Branches::nphe(0),10,1);
				for(int par = 0; par<3; par++){
					_hists->Histogram::Fid_Fill(envi,4,physics::get_theta(cz[par]),physics::get_phi(cx[par],cy[par]),par+1,6,1,_W,_p[par]);
					_hists->Histogram::DT_Fill(envi,4,par,_p[par], d[par], t[par], data->Branches::sc_r(0), data->Branches::sc_t(0),6,1,_W,physics::get_sector(physics::get_phi(cx[par],cy[par])));
				}
			}
		}
	}
	//std::cout<<std::endl <<"topologies hit: " <<"pmiss " <<envi->was_top(0) <<" pipmis " <<envi->was_top(1) <<" pimmiss " <<envi->was_top(2) <<" zmiss " <<envi->was_top(3) <<std::endl;
	//Assign topology of the event and fill things for combined topologies
	if(topo[3]){
		_top = 4; 
		//std::cout<<std::endl <<"All topology is go";
	}else{
		for(int o =0;o<3;o++){
			if(topo[o]){
				_top = o+1;
				//std::cout<<std::endl <<_top <<" topology is go";
			}
		}
		switch(_top){
			case 1:
				_pro = _beam + _target - _elec - _pip - _pim;
			break;
			case 2:
				_pip = _beam + _target - _elec - _pro - _pim;
			break;
			case 3:
				_pim = _beam + _target - _elec - _pro - _pip;
			break;
		}
	}

	if(_top!=0){
		//std::cout<<" and we made it into the combined event selection?" <<std::endl;
		_valid = true;
		_hists->Histogram::WQ2_Fill(envi,5,10,_W,_Q2);
		_hists->Histogram::Fid_Fill(envi,5,physics::get_theta(data->Branches::cz(0)),physics::get_phi(data->Branches::cx(0),data->Branches::cy(0)),0,10,0,_W,data->Branches::p(0));
		_hists->Histogram::SF_Fill(envi,5,data->Branches::p(0),data->Branches::etot(0),10,0,_W,sector[0]);
		_hists->Histogram::CC_Fill(envi,5,data->Branches::cc_sect(0),data->Branches::cc_segm(0),data->Branches::nphe(0),10,0);
		for(int par = 0; par<3; par++){
			if(topo[3]){
				_hists->Histogram::Fid_Fill(envi,5,physics::get_theta(cz[par]),physics::get_phi(cx[par],cy[par]),par+1,6,0,_W,_p[par]);
				_hists->Histogram::DT_Fill(envi,5,par,_p[par], d[par], t[par], data->Branches::sc_r(0), data->Branches::sc_t(0),6,0,_W,physics::get_sector(physics::get_phi(cx[par],cy[par])));
				//std::cout<<std::endl<<"None Missing" <<std::endl;
			}else if(topo[0] && par !=0){
				_hists->Histogram::Fid_Fill(envi,5,physics::get_theta(cz[par]),physics::get_phi(cx[par],cy[par]),par+1,6,0,_W,_p[par]);
				_hists->Histogram::DT_Fill(envi,5,par,_p[par], d[par], t[par], data->Branches::sc_r(0), data->Branches::sc_t(0),6,0,_W,physics::get_sector(physics::get_phi(cx[par],cy[par])));
				//std::cout<<"p Missing" <<std::endl;
				//std::cout<<"look at vals for particle "<<species[par+1] <<": _p: " <<_p[par] <<" d: " <<d[par] <<" t: " <<t[par] <<" sec: " <<physics::get_sector(physics::get_phi(cx[par],cy[par]));
			}else if(topo[2] && par!=1){
				_hists->Histogram::Fid_Fill(envi,5,physics::get_theta(cz[par]),physics::get_phi(cx[par],cy[par]),par+1,6,0,_W,_p[par]);
				_hists->Histogram::DT_Fill(envi,5,par,_p[par], d[par], t[par], data->Branches::sc_r(0), data->Branches::sc_t(0),6,0,_W,physics::get_sector(physics::get_phi(cx[par],cy[par])));	
				//std::cout<<std::endl<<"pip Missing" <<std::endl;
				//std::cout<<"look at vals for particle "<<species[par+1] <<": _p: " <<_p[par] <<" d: " <<d[par] <<" t: " <<t[par] <<" sec: " <<physics::get_sector(physics::get_phi(cx[par],cy[par]));
			}else if(topo[3] && par!=2){
				_hists->Histogram::Fid_Fill(envi,5,physics::get_theta(cz[par]),physics::get_phi(cx[par],cy[par]),par+1,6,0,_W,_p[par]);
				_hists->Histogram::DT_Fill(envi,5,par,_p[par], d[par], t[par], data->Branches::sc_r(0), data->Branches::sc_t(0),6,0,_W,physics::get_sector(physics::get_phi(cx[par],cy[par])));
				//std::cout<<std::endl<<"pim Missing" <<std::endl;
				//std::cout<<"look at vals for particle "<<species[par+1] <<": _p: " <<_p[par] <<" d: " <<d[par] <<" t: " <<t[par] <<" sec: " <<physics::get_sector(physics::get_phi(cx[par],cy[par]));
			}
		}
	}else if(top_possible[0] || top_possible[1] || top_possible[2] || top_possible[3]){//Didn't pass event selection, but still had the relevant particles detected
		_hists->Histogram::Fid_Fill(envi,5,physics::get_theta(data->Branches::cz(0)),physics::get_phi(data->Branches::cx(0),data->Branches::cy(0)),0,10,1,_W,data->Branches::p(0));
		_hists->Histogram::SF_Fill(envi,5,data->Branches::p(0),data->Branches::etot(0),10,1,_W,sector[0]);
		_hists->Histogram::CC_Fill(envi,5,data->Branches::cc_sect(0),data->Branches::cc_segm(0),data->Branches::nphe(0),10,1);
		for(int par = 0; par<3; par++){
			if(top_possible[3]){//None Missing Topology
				_hists->Histogram::Fid_Fill(envi,5,physics::get_theta(cz[par]),physics::get_phi(cx[par],cy[par]),par+1,6,1,_W,_p[par]);
				_hists->Histogram::DT_Fill(envi,5,par,_p[par], d[par], t[par], data->Branches::sc_r(0), data->Branches::sc_t(0),6,1,_W,physics::get_sector(physics::get_phi(cx[par],cy[par])));
			}else if(top_possible[0] && par!=0){//Pro Missing
				_hists->Histogram::Fid_Fill(envi,5,physics::get_theta(cz[par]),physics::get_phi(cx[par],cy[par]),par+1,6,1,_W,_p[par]);
				_hists->Histogram::DT_Fill(envi,5,par,_p[par], d[par], t[par], data->Branches::sc_r(0), data->Branches::sc_t(0),6,1,_W,physics::get_sector(physics::get_phi(cx[par],cy[par])));
			}else if(top_possible[1] && par!=1){//Pip Missing
				_hists->Histogram::Fid_Fill(envi,5,physics::get_theta(cz[par]),physics::get_phi(cx[par],cy[par]),par+1,6,1,_W,_p[par]);
				_hists->Histogram::DT_Fill(envi,5,par,_p[par], d[par], t[par], data->Branches::sc_r(0), data->Branches::sc_t(0),6,1,_W,physics::get_sector(physics::get_phi(cx[par],cy[par])));	
			}else if(top_possible[2] && par!=2){//Pim Missing
				_hists->Histogram::Fid_Fill(envi,5,physics::get_theta(cz[par]),physics::get_phi(cx[par],cy[par]),par+1,6,1,_W,_p[par]);
				_hists->Histogram::DT_Fill(envi,5,par,_p[par], d[par], t[par], data->Branches::sc_r(0), data->Branches::sc_t(0),6,1,_W,physics::get_sector(physics::get_phi(cx[par],cy[par])));
			}
		}
	}
}

float Event_Class::Get_px(int i){
	float _px = -99; 
	switch(i){
		case 0:
			_px = _beam[0];
		break;
		case 1:
			_px = _target[0];
		break;
		case 2:
			_px = _elec[0];
		break;
		case 3:
			_px = _pro[0];
		break;
		case 4:
			_px = _pip[0];
		break;
		case 5:
			_px = _pim[0];
		break;
	}
	return _px;
}

float Event_Class::Get_py(int i){
	float _py = -99; 
	switch(i){
		case 0:
			_py = _beam[1];
		break;
		case 1:
			_py = _target[1];
		break;
		case 2:
			_py = _elec[1];
		break;
		case 3:
			_py = _pro[1];
		break;
		case 4:
			_py = _pip[1];
		break;
		case 5:
			_py = _pim[1];
		break;
	}
	return _py;
}

float Event_Class::Get_pz(int i)
{
	float _pz = -99; 
	switch(i){
		case 0:
			_pz = _beam[2];
		break;
		case 1:
			_pz = _target[2];
		break;
		case 2:
			_pz = _elec[2];
		break;
		case 3:
			_pz = _pro[2];
		break;
		case 4:
			_pz = _pip[2];
		break;
		case 5:
			_pz = _pim[2];
		break;
	}
	return _pz;
}

float Event_Class::Get_p0(int i)
{
	float _p0 = -99; 
	switch(i){
		case 0:
			_p0 = _beam[3];
		break;
		case 1:
			_p0 = _target[3];
		break;
		case 2:
			_p0 = _elec[3];
		break;
		case 3:
			_p0 = _pro[3];
		break;
		case 4:
			_p0 = _pip[3];
		break;
		case 5:
			_p0 = _pim[3];
		break;
	}
	return _p0;
}

float Event_Class::Get_hel(){
	return _helicity; 
}

float Event_Class::Get_top(){
	return _top; 
}

float Event_Class::Get_pid(int i){
	int _pid = 0;
	switch(i){
		case 0:
			_pid = ELECTRON;
		break;
		case 1:
			_pid = PROTON;
		break;
		case 2:
			_pid = ELECTRON;
		break;
		case 3:
			_pid = PROTON;
		break;
		case 4:
			_pid = PION;
		break;
		case 5:
			_pid = -PION;
		break;
		default:
			_pid = 0; 
		break;
	}
	return _pid; 
}

bool Event_Class::is_valid(){
	return _valid; 
}

int Event_Class::Get_run_type(){
	return _run_type;
}

int Event_Class::Get_ppip(int idx){
	return ppip[idx];
}

/*
void Fill_Tree(forest tree, int Event_Class_n){
	tree.forest::fill_evnt(Event_Class_n);
	tree.forest::fill_apart(4);//Four particles
	for(int i = 0; i<4; i++){
		switch(i){
		case 0:
			tree.forest::fill_px(_elec[0], i);
			tree.forest::fill_py(_elec[1], i);
			tree.forest::fill_pz(_elec[2], i);
			tree.forest::fill_p0(_elec[3], i);
			tree.forest::fill_pid(ELECTRON, i);
		break;
		case 1:
			tree.forest::fill_px(_pro[0], i);
			tree.forest::fill_py(_pro[1], i);
			tree.forest::fill_pz(_pro[2], i);
			tree.forest::fill_p0(_pro[3], i);
			tree.forest::fill_pid(PROTON, i);
		break;
		case 2:
			tree.forest::fill_px(_pip[0], i);
			tree.forest::fill_py(_pip[1], i);
			tree.forest::fill_pz(_pip[2], i);
			tree.forest::fill_p0(_pip[3], i);
			tree.forest::fill_pid(PION, i);
		break;
		case 3:
			tree.forest::fill_px(_pim[0], i);
			tree.forest::fill_py(_pim[1], i);
			tree.forest::fill_pz(_pim[2], i);
			tree.forest::fill_p0(_pim[3], i);
			tree.forest::fill_pid(-PION, i);
		break;
		}	
	}
	tree.forest::fill_hel(_hel);
	tree.forest::fill_top(_top);
}*/
