#include "event_class.hpp"


Event_Class::Event_Class(std::shared_ptr<Branches> data, std::shared_ptr<Histogram> _hists, int run_type){ 
		bool fid_e_pass = false;
		bool sf_e_pass = false;
		bool cc_e_pass = false;
		int num_parts = data->Branches::gpart();
		float theta[num_parts];
		float phi[num_parts];
		float sector[num_parts];
		float _p[num_parts];
		//std::cout<<std::endl <<"gpart: " <<num_parts <<std::endl;
		for(int i=0; i < num_parts; i++){
			theta[i] = physics::get_theta(data->Branches::cz(i));
			phi[i] = physics::get_phi(data->Branches::cx(i),data->Branches::cy(i));
			sector[i] = physics::get_sector(phi[i]);
			_p[i] = data->Branches::p(i);
			//std::cout<<"theta " <<i <<": " <<theta[i] <<std::endl;
			//std::cout<<"phi " <<i <<": " <<phi[i]<<std::endl;
			//std::cout<<"sector " <<i <<": " <<sector[i]<<std::endl;
			//std::cout<<"momentum " <<i <<": " <<_p[i]<<std::endl;
		}

		_W = physics::WP(0,data);
		_Q2 = physics::Qsquared(0,data);
		//Pre ID Filling
			//Event-Wide
			//Electrons
		_hists->Histogram::WQ2_Fill(0,0,_W,_Q2);
		_hists->Histogram::Fid_Fill(0,theta[0],phi[0],0,0,0,_W,_p[0]);
		_hists->Histogram::SF_Fill(0,data->Branches::p(0),data->Branches::etot(0),0,0,_W,sector[0]);
		_hists->Histogram::CC_Fill(0,data->Branches::cc_sect(0),data->Branches::cc_segm(0),data->Branches::nphe(0),0,0);
		
		
		//electron ID
		if(data->Branches::q(0)==-1 && data->Branches::cc(0)!=0 && data->Branches::dc(0)!=0 && data->Branches::sc(0)!=0 && data->Branches::ec(0)!=0){//Sanity are q = -1, cc, dc, ec, and sc hits
			_hists->Histogram::WQ2_Fill(0,1,_W,_Q2);
			_hists->Histogram::Fid_Fill(0,theta[0],phi[0],0,1,0,_W,data->Branches::p(0));
			_hists->Histogram::SF_Fill(0,data->Branches::p(0),data->Branches::etot(0),1,0,_W,sector[0]);
			_hists->Histogram::CC_Fill(0,data->Branches::cc_sect(0),data->Branches::cc_segm(0),data->Branches::nphe(0),1,0);
			if(cuts::fid_cut(0,data->Branches::p(0),data->Branches::cx(0),data->Branches::cy(0),data->Branches::cz(0))){
				fid_e_pass = true;
				_hists->Histogram::WQ2_Fill(0,2,_W,_Q2);
				_hists->Histogram::Fid_Fill(0,theta[0],phi[0],0,2,0,_W,data->Branches::p(0));
				_hists->Histogram::SF_Fill(0,data->Branches::p(0),data->Branches::etot(0),2,0,_W,sector[0]);
				_hists->Histogram::CC_Fill(0,data->Branches::cc_sect(0),data->Branches::cc_segm(0),data->Branches::nphe(0),2,0);
			}else{
				_hists->Histogram::Fid_Fill(0,theta[0],phi[0],0,2,1,_W,data->Branches::p(0));
				_hists->Histogram::SF_Fill(0,data->Branches::p(0),data->Branches::etot(0),2,1,_W,sector[0]);
				_hists->Histogram::CC_Fill(0,data->Branches::cc_sect(0),data->Branches::cc_segm(0),data->Branches::nphe(0),2,1);
			}
			if(cuts::sf_cut(data->Branches::p(0),data->Branches::etot(0),data->Branches::cx(0),data->Branches::cy(0))){
				sf_e_pass = true;
				_hists->Histogram::WQ2_Fill(0,3,_W,_Q2);
				_hists->Histogram::Fid_Fill(0,theta[0],phi[0],0,3,0,_W,data->Branches::p(0));
				_hists->Histogram::SF_Fill(0,data->Branches::p(0),data->Branches::etot(0),3,0,_W,sector[0]);
				_hists->Histogram::CC_Fill(0,data->Branches::cc_sect(0),data->Branches::cc_segm(0),data->Branches::nphe(0),3,0);
			}else{
				_hists->Histogram::Fid_Fill(0,theta[0],phi[0],0,3,1,_W,data->Branches::p(0));
				_hists->Histogram::SF_Fill(0,data->Branches::p(0),data->Branches::etot(0),3,1,_W,sector[0]);
				_hists->Histogram::CC_Fill(0,data->Branches::cc_sect(0),data->Branches::cc_segm(0),data->Branches::nphe(0),3,1);
			}
			if(cuts::min_cc(data->Branches::cc_segm(0),data->Branches::cc_sect(0),data->Branches::nphe(0))){
				cc_e_pass = true;
				_hists->Histogram::WQ2_Fill(0,4,_W,_Q2);
				_hists->Histogram::Fid_Fill(0,theta[0],phi[0],0,4,0,_W,data->Branches::p(0));
				_hists->Histogram::SF_Fill(0,data->Branches::p(0),data->Branches::etot(0),4,0,_W,sector[0]);
				_hists->Histogram::CC_Fill(0,data->Branches::cc_sect(0),data->Branches::cc_segm(0),data->Branches::nphe(0),4,0);
			}else{
				_hists->Histogram::Fid_Fill(0,theta[0],phi[0],0,4,1,_W,data->Branches::p(0));
				_hists->Histogram::SF_Fill(0,data->Branches::p(0),data->Branches::etot(0),4,1,_W,sector[0]);
				_hists->Histogram::CC_Fill(0,data->Branches::cc_sect(0),data->Branches::cc_segm(0),data->Branches::nphe(0),4,1);
			}
			if(fid_e_pass && sf_e_pass){
				_hists->Histogram::WQ2_Fill(0,5,_W,_Q2);
				_hists->Histogram::Fid_Fill(0,theta[0],phi[0],0,5,0,_W,data->Branches::p(0));
				_hists->Histogram::SF_Fill(0,data->Branches::p(0),data->Branches::etot(0),5,0,_W,sector[0]);
				_hists->Histogram::CC_Fill(0,data->Branches::cc_sect(0),data->Branches::cc_segm(0),data->Branches::nphe(0),5,0);
			}else{
				_hists->Histogram::Fid_Fill(0,theta[0],phi[0],0,5,1,_W,data->Branches::p(0));
				_hists->Histogram::SF_Fill(0,data->Branches::p(0),data->Branches::etot(0),5,1,_W,sector[0]);
				_hists->Histogram::CC_Fill(0,data->Branches::cc_sect(0),data->Branches::cc_segm(0),data->Branches::nphe(0),5,1);
			}
			if(fid_e_pass && cc_e_pass){
				_hists->Histogram::WQ2_Fill(0,6,_W,_Q2);
				_hists->Histogram::Fid_Fill(0,theta[0],phi[0],0,6,0,_W,data->Branches::p(0));
				_hists->Histogram::SF_Fill(0,data->Branches::p(0),data->Branches::etot(0),6,0,_W,sector[0]);
				_hists->Histogram::CC_Fill(0,data->Branches::cc_sect(0),data->Branches::cc_segm(0),data->Branches::nphe(0),6,0);
			}else{
				_hists->Histogram::Fid_Fill(0,theta[0],phi[0],0,6,1,_W,data->Branches::p(0));
				_hists->Histogram::SF_Fill(0,data->Branches::p(0),data->Branches::etot(0),6,1,_W,sector[0]);
				_hists->Histogram::CC_Fill(0,data->Branches::cc_sect(0),data->Branches::cc_segm(0),data->Branches::nphe(0),6,1);
			}
			if(sf_e_pass && cc_e_pass){
				_hists->Histogram::WQ2_Fill(0,7,_W,_Q2);
				_hists->Histogram::Fid_Fill(0,theta[0],phi[0],0,7,0,_W,data->Branches::p(0));
				_hists->Histogram::SF_Fill(0,data->Branches::p(0),data->Branches::etot(0),7,0,_W,sector[0]);
				_hists->Histogram::CC_Fill(0,data->Branches::cc_sect(0),data->Branches::cc_segm(0),data->Branches::nphe(0),7,0);
			}else{
				_hists->Histogram::Fid_Fill(0,theta[0],phi[0],0,7,1,_W,data->Branches::p(0));
				_hists->Histogram::SF_Fill(0,data->Branches::p(0),data->Branches::etot(0),7,1,_W,sector[0]);
				_hists->Histogram::CC_Fill(0,data->Branches::cc_sect(0),data->Branches::cc_segm(0),data->Branches::nphe(0),7,1);
			}
			if(fid_e_pass && cc_e_pass && sf_e_pass){
				_hists->Histogram::WQ2_Fill(0,8,_W,_Q2);
				_hists->Histogram::Fid_Fill(0,theta[0],phi[0],0,8,0,_W,data->Branches::p(0));
				_hists->Histogram::SF_Fill(0,data->Branches::p(0),data->Branches::etot(0),8,0,_W,sector[0]);
				_hists->Histogram::CC_Fill(0,data->Branches::cc_sect(0),data->Branches::cc_segm(0),data->Branches::nphe(0),8,0);
			}else{
				_hists->Histogram::Fid_Fill(0,theta[0],phi[0],0,8,1,_W,data->Branches::p(0));
				_hists->Histogram::SF_Fill(0,data->Branches::p(0),data->Branches::etot(0),8,1,_W,sector[0]);
				_hists->Histogram::CC_Fill(0,data->Branches::cc_sect(0),data->Branches::cc_segm(0),data->Branches::nphe(0),8,1);
			}
			if((data->Branches::id(0))==ELECTRON){
				_hists->Histogram::WQ2_Fill(0,9,_W,_Q2);
				_hists->Histogram::Fid_Fill(0,theta[0],phi[0],0,9,0,_W,data->Branches::p(0));
				_hists->Histogram::SF_Fill(0,data->Branches::p(0),data->Branches::etot(0),9,0,_W,sector[0]);
				_hists->Histogram::CC_Fill(0,data->Branches::cc_sect(0),data->Branches::cc_segm(0),data->Branches::nphe(0),9,0);
			}
			else{
				_hists->Histogram::Fid_Fill(0,theta[0],phi[0],0,9,1,_W,data->Branches::p(0));
				_hists->Histogram::SF_Fill(0,data->Branches::p(0),data->Branches::etot(0),9,1,_W,sector[0]);
				_hists->Histogram::CC_Fill(0,data->Branches::cc_sect(0),data->Branches::cc_segm(0),data->Branches::nphe(0),9,1);
			}
		}else{
			_hists->Histogram::Fid_Fill(0,theta[0],phi[0],0,1,1,_W,data->Branches::p(0));
			_hists->Histogram::SF_Fill(0,data->Branches::p(0),data->Branches::etot(0),1,1,_W,sector[0]);
			_hists->Histogram::CC_Fill(0,data->Branches::cc_sect(0),data->Branches::cc_segm(0),data->Branches::nphe(0),1,1);
		}
		
		//Hadron Loop
		int q_h = 0; //desired charge for hadron
		bool fid_h = false; 
		bool dt_h = false;
		int hid_num = -99; 

		for(int h = 1; h < num_parts; h++){
			for(int part = 0; part < 3; part++){//particle loop
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
				_hists->Histogram::Fid_Fill(0,theta[h],phi[h],part+1,0,0,_W,data->Branches::p(h));
				//_hists->Histogram::DT_Fill(0,part,data->Branches::p(h), data->Branches::sc_r(h), data->Branches::sc_t(h), data->Branches::sc_r(0), data->Branches::sc_t(0),0,0,_W,sector[h]);
				//Sanity Cut
				if(data->Branches::q(h)==q_h && data->Branches::dc(h)!=0 && data->Branches::sc(h)!=0){
					_hists->Histogram::Fid_Fill(0,theta[h],phi[h],part+1,1,0,_W,data->Branches::p(h));
					//_hists->Histogram::DT_Fill(0,part,data->Branches::p(h), data->Branches::sc_r(h), data->Branches::sc_t(h), data->Branches::sc_r(0), data->Branches::sc_t(0),1,0,_W,sector[h]);
					if(cuts::fid_cut(part+1,data->Branches::p(h),data->Branches::cx(h),data->Branches::cy(h),data->Branches::cz(h))){
						fid_h = true; 
						_hists->Histogram::Fid_Fill(0,theta[h],phi[h],part+1,2,0,_W,data->Branches::p(h));
						//_hists->Histogram::DT_Fill(0,part,data->Branches::p(h), data->Branches::sc_r(h), data->Branches::sc_t(h), data->Branches::sc_r(0), data->Branches::sc_t(0),2,0,_W,sector[h]);
					}else{
						_hists->Histogram::Fid_Fill(0,theta[h],phi[h],part+1,2,1,_W,data->Branches::p(h));
						//_hists->Histogram::DT_Fill(0,part,data->Branches::p(h), data->Branches::sc_r(h), data->Branches::sc_t(h), data->Branches::sc_r(0), data->Branches::sc_t(0),2,1,_W,sector[h]);
					}	
					if(cuts::delta_t(part, data->Branches::p(h), data->Branches::sc_r(0), data->Branches::sc_r(h), data->Branches::sc_t(0), data->Branches::sc_t(h))){
						dt_h = true; 
						_hists->Histogram::Fid_Fill(0,theta[h],phi[h],part+1,3,0,_W,data->Branches::p(h));
						//_hists->Histogram::DT_Fill(0,part,data->Branches::p(h), data->Branches::sc_r(h), data->Branches::sc_t(h), data->Branches::sc_r(0), data->Branches::sc_t(0),3,0,_W,sector[h]);
					}else{
						_hists->Histogram::Fid_Fill(0,theta[h],phi[h],part+1,3,1,_W,data->Branches::p(h));
						//_hists->Histogram::DT_Fill(0,part,data->Branches::p(h), data->Branches::sc_r(h), data->Branches::sc_t(h), data->Branches::sc_r(0), data->Branches::sc_t(0),3,1,_W,sector[h]);
					}
					if(dt_h && fid_h){
						_hists->Histogram::Fid_Fill(0,theta[h],phi[h],part+1,4,0,_W,data->Branches::p(h));
						//_hists->Histogram::DT_Fill(0,part,data->Branches::p(h), data->Branches::sc_r(h), data->Branches::sc_t(h), data->Branches::sc_r(0), data->Branches::sc_t(0),4,0,_W,sector[h]);
					}else{
						_hists->Histogram::Fid_Fill(0,theta[h],phi[h],part+1,4,1,_W,data->Branches::p(h));
						//_hists->Histogram::DT_Fill(0,part,data->Branches::p(h), data->Branches::sc_r(h), data->Branches::sc_t(h), data->Branches::sc_r(0), data->Branches::sc_t(0),4,1,_W,sector[h]);
					}
					if(data->Branches::id(h) == hid_num){
						_hists->Histogram::Fid_Fill(0,theta[h],phi[h],part+1,5,0,_W,data->Branches::p(h));
						//_hists->Histogram::DT_Fill(0,part,data->Branches::p(h), data->Branches::sc_r(h), data->Branches::sc_t(h), data->Branches::sc_r(0), data->Branches::sc_t(0),5,0,_W,sector[h]);
					}else{
						_hists->Histogram::Fid_Fill(0,theta[h],phi[h],part+1,5,1,_W,data->Branches::p(h));
						//_hists->Histogram::DT_Fill(0,part,data->Branches::p(h), data->Branches::sc_r(h), data->Branches::sc_t(h), data->Branches::sc_r(0), data->Branches::sc_t(0),5,1,_W,sector[h]);
					}
				}else{
					_hists->Histogram::Fid_Fill(0,theta[h],phi[h],part+1,1,1,_W,data->Branches::p(h));
					//_hists->Histogram::DT_Fill(0,part,data->Branches::p(h), data->Branches::sc_r(h), data->Branches::sc_t(h), data->Branches::sc_r(0), data->Branches::sc_t(0),1,1,_W,sector[h]);
				}
			}
		}

		

		//Event_Class Selection


		//Variable Determination


}

float Event_Class::Get_px(int i){
	float _px = -99; 
	switch(i){
		case 0:
			_px = _elec[0];
		break;
		case 1:
			_px = _prot[0];
		break;
		case 2:
			_px = _pip[0];
		break;
		case 3:
			_px = _pim[0];
		break;
	}
	return _px;
}

float Event_Class::Get_py(int i){
	float _py = -99; 
	switch(i){
		case 0:
			_py = _elec[1];
		break;
		case 1:
			_py = _prot[1];
		break;
		case 2:
			_py = _pip[1];
		break;
		case 3:
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
			_pz = _elec[2];
		break;
		case 1:
			_pz = _prot[2];
		break;
		case 2:
			_pz = _pip[2];
		break;
		case 3:
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
			_p0 = _elec[3];
		break;
		case 1:
			_p0 = _prot[3];
		break;
		case 2:
			_p0 = _pip[3];
		break;
		case 3:
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
			_pid = PION;
		break;
		case 3:
			_pid = -PION;
		break;
		default:
			_pid = 10; 
		break;
	}
	return _pid; 
}

bool Event_Class::is_valid(){
	return _valid; 
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
			tree.forest::fill_px(_prot[0], i);
			tree.forest::fill_py(_prot[1], i);
			tree.forest::fill_pz(_prot[2], i);
			tree.forest::fill_p0(_prot[3], i);
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
