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

		_beam = k_mu_e16;
		_target = p_mu;
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
				_elec = physics::Make_4Vector(data->Branches::p(0),data->Branches::p(0),data->Branches::cy(0),data->Branches::cz(0),me);
				good_electron++;
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
		int p_loop = 0;
		int pip_loop = 0;
		int pim_loop = 0; 

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
				_hists->Histogram::Fid_Fill(0,theta[h],phi[h],part+1,0,0,_W,data->Branches::p(h));
				_hists->Histogram::DT_Fill(0,part,data->Branches::p(h), data->Branches::sc_r(h), data->Branches::sc_t(h), data->Branches::sc_r(0), data->Branches::sc_t(0),0,0,_W,sector[h]);
				//Sanity Cut
				if(data->Branches::q(h)==q_h && data->Branches::dc(h)!=0 && data->Branches::sc(h)!=0){
					_hists->Histogram::Fid_Fill(0,theta[h],phi[h],part+1,1,0,_W,data->Branches::p(h));
					_hists->Histogram::DT_Fill(0,part,data->Branches::p(h), data->Branches::sc_r(h), data->Branches::sc_t(h), data->Branches::sc_r(0), data->Branches::sc_t(0),1,0,_W,sector[h]);
					if(cuts::fid_cut(part+1,data->Branches::p(h),data->Branches::cx(h),data->Branches::cy(h),data->Branches::cz(h))){
						fid_h = true; 
						_hists->Histogram::Fid_Fill(0,theta[h],phi[h],part+1,2,0,_W,data->Branches::p(h));
						_hists->Histogram::DT_Fill(0,part,data->Branches::p(h), data->Branches::sc_r(h), data->Branches::sc_t(h), data->Branches::sc_r(0), data->Branches::sc_t(0),2,0,_W,sector[h]);
					}else{
						_hists->Histogram::Fid_Fill(0,theta[h],phi[h],part+1,2,1,_W,data->Branches::p(h));
						_hists->Histogram::DT_Fill(0,part,data->Branches::p(h), data->Branches::sc_r(h), data->Branches::sc_t(h), data->Branches::sc_r(0), data->Branches::sc_t(0),2,1,_W,sector[h]);
					}	
					if(cuts::delta_t_cut(part, data->Branches::p(h), data->Branches::sc_r(0), data->Branches::sc_r(h), data->Branches::sc_t(0), data->Branches::sc_t(h))){
						dt_h = true; 
						_hists->Histogram::Fid_Fill(0,theta[h],phi[h],part+1,3,0,_W,data->Branches::p(h));
						_hists->Histogram::DT_Fill(0,part,data->Branches::p(h), data->Branches::sc_r(h), data->Branches::sc_t(h), data->Branches::sc_r(0), data->Branches::sc_t(0),3,0,_W,sector[h]);
					}else{
						_hists->Histogram::Fid_Fill(0,theta[h],phi[h],part+1,3,1,_W,data->Branches::p(h));
						_hists->Histogram::DT_Fill(0,part,data->Branches::p(h), data->Branches::sc_r(h), data->Branches::sc_t(h), data->Branches::sc_r(0), data->Branches::sc_t(0),3,1,_W,sector[h]);
					}
					if(dt_h && fid_h){
						_hists->Histogram::Fid_Fill(0,theta[h],phi[h],part+1,4,0,_W,data->Branches::p(h));
						_hists->Histogram::DT_Fill(0,part,data->Branches::p(h), data->Branches::sc_r(h), data->Branches::sc_t(h), data->Branches::sc_r(0), data->Branches::sc_t(0),4,0,_W,sector[h]);
						switch(part){
							case 0:
							_prot = physics::Make_4Vector(data->Branches::p(h),data->Branches::p(h),data->Branches::cy(h),data->Branches::cz(h),mp);
							good_proton++;
							pro_idx[p_loop]=h;
							p_loop++;
							//std::cout<<" assign pro idx:" <<h <<std::endl;
							d[0] = data->Branches::sc_r(h);
							t[0] = data->Branches::sc_t(h);
							cx[0] = data->Branches::cx(h);
							cy[0] = data->Branches::cy(h);
							cz[0] = data->Branches::cz(h);
							_p[0] = data->Branches::p(h);
							h_sec[0] = sector[h];
							check_idx[0] = h; 
							break;
							case 1:
							_pip = physics::Make_4Vector(data->Branches::p(h),data->Branches::p(h),data->Branches::cy(h),data->Branches::cz(h),mpi);
							good_pip++;
							pip_idx[pip_loop]=h;
							pip_loop++;
							//std::cout<<" assign pip idx:" <<h <<std::endl;
							d[1] = data->Branches::sc_r(h);
							t[1] = data->Branches::sc_t(h);
							cx[1] = data->Branches::cx(h);
							cy[1] = data->Branches::cy(h);
							cz[1] = data->Branches::cz(h);
							_p[1] = data->Branches::p(h);
							h_sec[1] = sector[h];
							check_idx[1] = h; 
							break;
							case 2:
							_pim = physics::Make_4Vector(data->Branches::p(h),data->Branches::p(h),data->Branches::cy(h),data->Branches::cz(h),mpi);
							good_pim=h;
							pim_idx[pim_loop]=h;
							pim_loop++;
							check_idx[2] = h; 
							//std::cout<<" assign pim idx:" <<h <<std::endl;
							d[2] = data->Branches::sc_r(h);
							t[2] = data->Branches::sc_t(h);
							cx[2] = data->Branches::cx(h);
							cy[2] = data->Branches::cy(h);
							cz[2] = data->Branches::cz(h);
							_p[2] = data->Branches::p(h);
							h_sec[2] = sector[h];
							break;
						}
					}else{
						_hists->Histogram::Fid_Fill(0,theta[h],phi[h],part+1,4,1,_W,data->Branches::p(h));
						_hists->Histogram::DT_Fill(0,part,data->Branches::p(h), data->Branches::sc_r(h), data->Branches::sc_t(h), data->Branches::sc_r(0), data->Branches::sc_t(0),4,1,_W,sector[h]);
					}
					if(data->Branches::id(h) == hid_num){
						_hists->Histogram::Fid_Fill(0,theta[h],phi[h],part+1,5,0,_W,data->Branches::p(h));
						_hists->Histogram::DT_Fill(0,part,data->Branches::p(h), data->Branches::sc_r(h), data->Branches::sc_t(h), data->Branches::sc_r(0), data->Branches::sc_t(0),5,0,_W,sector[h]);
					}else{
						_hists->Histogram::Fid_Fill(0,theta[h],phi[h],part+1,5,1,_W,data->Branches::p(h));
						_hists->Histogram::DT_Fill(0,part,data->Branches::p(h), data->Branches::sc_r(h), data->Branches::sc_t(h), data->Branches::sc_r(0), data->Branches::sc_t(0),5,1,_W,sector[h]);
					}
				}else{
					_hists->Histogram::Fid_Fill(0,theta[h],phi[h],part+1,1,1,_W,data->Branches::p(h));
					_hists->Histogram::DT_Fill(0,part,data->Branches::p(h), data->Branches::sc_r(h), data->Branches::sc_t(h), data->Branches::sc_r(0), data->Branches::sc_t(0),1,1,_W,sector[h]);
				}
			}
		}
		
		bool double_id = false; 
		int double_idx = -99;
		int num_double_id = 0;

		bool have_prot = false;
		bool have_pip = false;
		bool have_pim = false;

		if(check_idx[2]==good_pim){
			have_pim = true; 
		}
		if(check_idx[1]==good_pip){
			have_pip = true; 
		}
		if(check_idx[0]==good_proton){
			have_prot = true; 
		}

		//std::cout<<"pro_idx:" <<pro_idx[0];
		//std::cout<<" pip_idx:" <<pip_idx[0];
		//std::cout<<" pim_idx:" <<pim_idx[0] <<std::endl;

		//std::cout<<std::endl <<"p_loop: " <<p_loop << "   pip_loop: " <<pip_loop <<std::endl;
		//std::cout<<"good_proton: " <<good_proton << "   good_pip: " <<good_pip <<std::endl;
		

		/*TLorentzVector test_p,test_pip;
		if(good_proton >1){//Detected more than 1 proton?
			if(good_pip > 1){//Detected more than 1 pip?
				for(int i = 0 ; i< p_loop ; i++){//Look at all indexes for proton
					for( int j = 0; j < pip_loop; j++){//Look at all indexes for proton
						if( i == j ){//Same particle IDed as both proton and pip
							double_id= true; 
							double_idx = i; 
							num_double_id++; 
						}
					}
				}
				if(num_double_id < 1){//IF there were no double IDs then default to the first proton or pip measured
					_pip = physics::Make_4Vector(data->Branches::p(pip_idx[0]),data->Branches::p(pip_idx[0]),data->Branches::cy(pip_idx[0]),data->Branches::cz(pip_idx[0]),mpi);
					_prot = physics::Make_4Vector(data->Branches::p(pro_idx[0]),data->Branches::p(pro_idx[0]),data->Branches::cy(pro_idx[0]),data->Branches::cz(pro_idx[0]),mp);
				}else{
					//This shouldn't be too big of an issue
					if(num_double_id > 1){
						ppip++;//WIll see how big an issue this is
					}
					if(num_double_id == 1){
						test_pip = physics::Make_4Vector(data->Branches::p(pip_idx[double_idx]),data->Branches::p(pip_idx[double_idx]),data->Branches::cy(pip_idx[double_idx]),data->Branches::cz(pip_idx[double_idx]),mpi);
						test_p = physics::Make_4Vector(data->Branches::p(pro_idx[double_idx]),data->Branches::p(pro_idx[double_idx]),data->Branches::cy(pro_idx[double_idx]),data->Branches::cz(pro_idx[double_idx]),mp);
					}
				}
			}else if(good_pip==1){	
				_pip = physics::Make_4Vector(data->Branches::p(pip_idx[0]),data->Branches::p(pip_idx[0]),data->Branches::cy(pip_idx[0]),data->Branches::cz(pip_idx[0]),mpi);
			}
		}else if(good_proton==1){
			_prot = physics::Make_4Vector(data->Branches::p(pro_idx[0]),data->Branches::p(pro_idx[0]),data->Branches::cy(pro_idx[0]),data->Branches::cz(pro_idx[0]),mp);
		}*/
		
		bool topo[4]={false,false,false,false};
		int idx=-99; 
		int par = -99; 

		//Event_Class Selection
		if(good_electron == 1){//Good Electron
			//P Missing
			if((have_pim) && have_pip){//Good Pip and Pim
				MM_p = physics::MM_event(0,0,_elec,_pip,_pim);
				MM_p2 = physics::MM_event(0,1,_elec,_pip,_pim);
				_hists->Histogram::MM_Fill(0,MM_p,0,0);
				_hists->Histogram::MM_Fill(0,MM_p2,0,1);
				_hists->Histogram::MM_Fill(4,MM_p,0,0);
				_hists->Histogram::MM_Fill(4,MM_p2,0,1);
				if((MM_p > (p_center-p_sig))&&(MM_p < (p_center+p_sig))){//Missing Mass Cut on Proton Mass
					_hists->Histogram::MM_Fill(0,MM_p,1,0);
					_hists->Histogram::MM_Fill(0,MM_p2,1,1);
					_hists->Histogram::MM_Fill(4,MM_p,1,0);
					_hists->Histogram::MM_Fill(4,MM_p2,1,1);
					topo[0]=true;
					_prot = _beam + _target - _elec - _pip - _pim; 
					_hists->Histogram::WQ2_Fill(1,10,_W,_Q2);
					_hists->Histogram::Fid_Fill(1,theta[0],phi[0],0,10,0,_W,data->Branches::p(0));
					_hists->Histogram::SF_Fill(1,data->Branches::p(0),data->Branches::etot(0),10,0,_W,sector[0]);
					_hists->Histogram::CC_Fill(1,data->Branches::cc_sect(0),data->Branches::cc_segm(0),data->Branches::nphe(0),10,0);
					for(int i = 0; i<2; i++){
						switch(i){
							case 0:
							idx = pip_idx[0];
							par = 1; 
							break;
							case 1:
							idx = pim_idx[0];
							par = 2; 
							break;
						}
						_hists->Histogram::Fid_Fill(1,physics::get_theta(cz[par]),physics::get_phi(cx[par],cy[par]),par+1,6,0,_W,data->Branches::p(idx));
						_hists->Histogram::DT_Fill(1,par,_p[par], d[par], t[par], data->Branches::sc_r(0), data->Branches::sc_t(0),6,0,_W,h_sec[par]);
					}
				}else{
					_hists->Histogram::MM_Fill(0,MM_p,2,0);
					_hists->Histogram::MM_Fill(0,MM_p2,2,1);
					_hists->Histogram::MM_Fill(4,MM_p,2,0);
					_hists->Histogram::MM_Fill(4,MM_p2,2,1);
					_hists->Histogram::WQ2_Fill(1,10,_W,_Q2);
					_hists->Histogram::Fid_Fill(1,theta[0],phi[0],0,10,1,_W,data->Branches::p(0));
					_hists->Histogram::SF_Fill(1,data->Branches::p(0),data->Branches::etot(0),10,1,_W,sector[0]);
					_hists->Histogram::CC_Fill(1,data->Branches::cc_sect(0),data->Branches::cc_segm(0),data->Branches::nphe(0),10,1);
					for(int i = 0; i<2; i++){
						switch(i){
							case 0:
							idx = pip_idx[0];
							par = 1; 
							break;
							case 1:
							idx = pim_idx[0];
							par = 2; 
							break;
						}
						_hists->Histogram::Fid_Fill(1,physics::get_theta(cz[par]),physics::get_phi(cx[par],cy[par]),par+1,6,1,_W,data->Branches::p(idx));
						_hists->Histogram::DT_Fill(1,par,_p[par], d[par], t[par], data->Branches::sc_r(0), data->Branches::sc_t(0),6,1,_W,h_sec[par]);
					}
				}
			}
			//Pip Missing
			if((have_prot) && (have_pim)){//Good proton and Pim
				MM_pip = physics::MM_event(0,0,_elec,_prot,_pim);
				MM_pip2 = physics::MM_event(0,1,_elec,_prot,_pim);
				_hists->Histogram::MM_Fill(1,MM_pip,0,0);
				_hists->Histogram::MM_Fill(1,MM_pip2,0,1);
				_hists->Histogram::MM_Fill(4,MM_pip,0,0);
				_hists->Histogram::MM_Fill(4,MM_pip2,0,1);
				if((MM_pip > (pip_center-pip_sig))&&(MM_pip < (pip_center+pip_sig))){//Missing Mass Cut on proton Mass
					_hists->Histogram::MM_Fill(1,MM_pip,1,0);
					_hists->Histogram::MM_Fill(1,MM_pip2,1,1);
					_hists->Histogram::MM_Fill(4,MM_pip,1,0);
					_hists->Histogram::MM_Fill(4,MM_pip2,1,1);
					topo[1]=true;
					_pip = _beam + _target - _elec - _prot - _pim; 
					//Fill Event Selection 
					//Electron
					_hists->Histogram::WQ2_Fill(2,10,_W,_Q2);
					_hists->Histogram::Fid_Fill(2,theta[0],phi[0],0,10,0,_W,data->Branches::p(0));
					_hists->Histogram::SF_Fill(2,data->Branches::p(0),data->Branches::etot(0),10,0,_W,sector[0]);
					_hists->Histogram::CC_Fill(2,data->Branches::cc_sect(0),data->Branches::cc_segm(0),data->Branches::nphe(0),10,0);
					//Hadron
					for(int i = 0; i<2; i++){
						switch(i){
							case 0:
							idx = pro_idx[0];
							par = 0; 
							break;
							case 1:
							idx = pim_idx[0];
							par = 2; 
							break;
						}
						_hists->Histogram::Fid_Fill(2,physics::get_theta(cz[par]),physics::get_phi(cx[par],cy[par]),par+1,6,0,_W,data->Branches::p(idx));
						_hists->Histogram::DT_Fill(2,par,_p[par], d[par], t[par], data->Branches::sc_r(0), data->Branches::sc_t(0),6,0,_W,h_sec[par]);
					}
				}else{
					_hists->Histogram::MM_Fill(1,MM_pip,2,0);
					_hists->Histogram::MM_Fill(1,MM_pip2,2,1);
					_hists->Histogram::MM_Fill(4,MM_pip,2,0);
					_hists->Histogram::MM_Fill(4,MM_pip2,2,1);
					_hists->Histogram::WQ2_Fill(2,10,_W,_Q2);
					_hists->Histogram::Fid_Fill(2,theta[0],phi[0],0,10,1,_W,data->Branches::p(0));
					_hists->Histogram::SF_Fill(2,data->Branches::p(0),data->Branches::etot(0),10,1,_W,sector[0]);
					_hists->Histogram::CC_Fill(2,data->Branches::cc_sect(0),data->Branches::cc_segm(0),data->Branches::nphe(0),10,1);
					for(int i = 0; i<2; i++){
						switch(i){
							case 0:
							idx = pro_idx[0];
							par = 0; 
							break;
							case 1:
							idx = pim_idx[0];
							par = 2; 
							break;
						}
						_hists->Histogram::Fid_Fill(2,physics::get_theta(cz[par]),physics::get_phi(cx[par],cy[par]),par+1,6,1,_W,data->Branches::p(idx));
						_hists->Histogram::DT_Fill(2,par,_p[par], d[par], t[par], data->Branches::sc_r(0), data->Branches::sc_t(0),6,1,_W,h_sec[par]);
					}
				}
			}
			//Pim Missing //Need to introduce some proton/Pion discrimination here, but will leave for now
			if((have_prot) && (have_pip) && !double_id){//Good Proton and Pip
				MM_pim = physics::MM_event(0,0,_elec,_prot,_pip);
				MM_pim2 = physics::MM_event(0,1,_elec,_prot,_pip);
				_hists->Histogram::MM_Fill(2,MM_pim,0,0);
				_hists->Histogram::MM_Fill(2,MM_pim2,0,1);
				_hists->Histogram::MM_Fill(4,MM_pim,0,0);
				_hists->Histogram::MM_Fill(4,MM_pim2,0,1);
				if((MM_pim > (pim_center-pim_sig))&&(MM_pim < (pim_center+pim_sig))){//Missing Mass Cut on Proton Mass
					_hists->Histogram::MM_Fill(2,MM_pim,1,0);
					_hists->Histogram::MM_Fill(2,MM_pim2,1,1);
					_hists->Histogram::MM_Fill(4,MM_pim,1,0);
					_hists->Histogram::MM_Fill(4,MM_pim2,1,1);
					topo[2]=true;
					_hists->Histogram::WQ2_Fill(3,10,_W,_Q2);
					_hists->Histogram::Fid_Fill(3,theta[0],phi[0],0,10,0,_W,data->Branches::p(0));
					_hists->Histogram::SF_Fill(3,data->Branches::p(0),data->Branches::etot(0),10,0,_W,sector[0]);
					_hists->Histogram::CC_Fill(3,data->Branches::cc_sect(0),data->Branches::cc_segm(0),data->Branches::nphe(0),10,0);
					_pim = _beam + _target - _elec - _prot - _pip; 
					for(int i = 0; i<2; i++){
						switch(i){
							case 0:
							idx = pro_idx[0];
							par = 0; 
							break;
							case 1:
							idx = pip_idx[0];
							par = 1; 
							break;
						}
						_hists->Histogram::Fid_Fill(3,physics::get_theta(cz[par]),physics::get_phi(cx[par],cy[par]),par+1,6,0,_W,data->Branches::p(idx));
						_hists->Histogram::DT_Fill(3,par,_p[par], d[par], t[par], data->Branches::sc_r(0), data->Branches::sc_t(0),6,0,_W,h_sec[par]);
					}
				}else{
					_hists->Histogram::MM_Fill(2,MM_pim,2,0);
					_hists->Histogram::MM_Fill(2,MM_pim2,2,1);
					_hists->Histogram::MM_Fill(4,MM_pim,2,0);
					_hists->Histogram::MM_Fill(4,MM_pim2,2,1);
					_hists->Histogram::WQ2_Fill(3,10,_W,_Q2);
					_hists->Histogram::Fid_Fill(3,theta[0],phi[0],0,10,1,_W,data->Branches::p(0));
					_hists->Histogram::SF_Fill(3,data->Branches::p(0),data->Branches::etot(0),10,1,_W,sector[0]);
					_hists->Histogram::CC_Fill(3,data->Branches::cc_sect(0),data->Branches::cc_segm(0),data->Branches::nphe(0),10,1);
					for(int i = 0; i<2; i++){
						switch(i){
							case 0:
							idx = pro_idx[0];
							par = 0; 
							break;
							case 1:
							idx = pip_idx[0];
							par = 1; 
							break;
						}
						_hists->Histogram::Fid_Fill(3,physics::get_theta(cz[par]),physics::get_phi(cx[par],cy[par]),par+1,6,1,_W,data->Branches::p(idx));
						_hists->Histogram::DT_Fill(3,par,_p[par], d[par], t[par], data->Branches::sc_r(0), data->Branches::sc_t(0),6,1,_W,h_sec[par]);
					}
				}
			}
			//None Missing
			if((have_prot)&&(have_pim)&&(have_pip) && !double_id){//Good Proton Pip and Pim 
				MM_z = physics::MM_event(0,0,_elec,_prot,_pip,_pim);
				MM_z2 = physics::MM_event(0,1,_elec,_prot,_pip,_pim);
				_hists->Histogram::MM_Fill(3,MM_z,0,0);
				_hists->Histogram::MM_Fill(3,MM_z2,0,1);
				_hists->Histogram::MM_Fill(4,MM_z,0,0);
				_hists->Histogram::MM_Fill(4,MM_z2,0,1);
				if((MM_z2 > (MM_zero_center2-MM_zero_sigma2))&&(MM_z2 < (MM_zero_center2+MM_zero_sigma2))){//Missing Mass Cut on Proton Mass
					_hists->Histogram::MM_Fill(3,MM_z,1,0);
					_hists->Histogram::MM_Fill(3,MM_z2,1,1);
					_hists->Histogram::MM_Fill(4,MM_z,1,0);
					_hists->Histogram::MM_Fill(4,MM_z2,1,1);
					topo[3]=true;
					_hists->Histogram::WQ2_Fill(4,10,_W,_Q2);
					_hists->Histogram::Fid_Fill(4,theta[0],phi[0],0,10,0,_W,data->Branches::p(0));
					_hists->Histogram::SF_Fill(4,data->Branches::p(0),data->Branches::etot(0),10,0,_W,sector[0]);
					_hists->Histogram::CC_Fill(4,data->Branches::cc_sect(0),data->Branches::cc_segm(0),data->Branches::nphe(0),10,0);
					for(int i = 0; i<3; i++){
						switch(i){
							case 0:
							idx = pro_idx[0];
							par = 0; 
							break;
							case 1:
							idx = pip_idx[0];
							par = 1; 
							break;
							case 2:
							idx = pim_idx[0];
							par = 2;
							break;
						}
						_hists->Histogram::Fid_Fill(4,physics::get_theta(cz[par]),physics::get_phi(cx[par],cy[par]),par+1,6,0,_W,data->Branches::p(idx));
						_hists->Histogram::DT_Fill(4,par,_p[par], d[par], t[par], data->Branches::sc_r(0), data->Branches::sc_t(0),6,0,_W,h_sec[par]);
					}
				}else{
					_hists->Histogram::MM_Fill(3,MM_z,2,0);
					_hists->Histogram::MM_Fill(3,MM_z2,2,1);
					_hists->Histogram::MM_Fill(4,MM_z,2,0);
					_hists->Histogram::MM_Fill(4,MM_z2,2,1);
					_hists->Histogram::WQ2_Fill(4,10,_W,_Q2);
					_hists->Histogram::Fid_Fill(4,theta[0],phi[0],0,10,1,_W,data->Branches::p(0));
					_hists->Histogram::SF_Fill(4,data->Branches::p(0),data->Branches::etot(0),10,1,_W,sector[0]);
					_hists->Histogram::CC_Fill(4,data->Branches::cc_sect(0),data->Branches::cc_segm(0),data->Branches::nphe(0),10,1);
					for(int i = 0; i<3; i++){
						switch(i){
							case 0:
							idx = pro_idx[0];
							par = 0; 
							break;
							case 1:
							idx = pip_idx[0];
							par = 1; 
							break;
							case 2:
							idx = pim_idx[0];
							par = 2;
							break;
						}
						_hists->Histogram::Fid_Fill(4,physics::get_theta(cz[par]),physics::get_phi(cx[par],cy[par]),par+1,6,1,_W,data->Branches::p(idx));
						_hists->Histogram::DT_Fill(4,par,_p[par], d[par], t[par], data->Branches::sc_r(0), data->Branches::sc_t(0),6,1,_W,h_sec[par]);
					}
				}
			}
		}
		//Assign topology of the event
		if(topo[3]){
			_top = 4; 
		}else{
			for(int o =0;o<3;o++){
				if(topo[o]){
					_top = topo[o]+1;
				}
			}
		}
		if(_top!=0){
			_valid = true;
			_hists->Histogram::WQ2_Fill(5,10,_W,_Q2);
			_hists->Histogram::Fid_Fill(5,theta[0],phi[0],0,10,0,_W,data->Branches::p(0));
			_hists->Histogram::SF_Fill(5,data->Branches::p(0),data->Branches::etot(0),10,0,_W,sector[0]);
			_hists->Histogram::CC_Fill(5,data->Branches::cc_sect(0),data->Branches::cc_segm(0),data->Branches::nphe(0),10,0);
			for(int i = 0; i<3; i++){
				switch(i){
					case 0:
					idx = pro_idx[0];
					par = 0; 
					break;
					case 1:
					idx = pip_idx[0];
					par = 1; 
					break;
					case 2:
					idx = pim_idx[0];
					par = 2;
					break;
				}
				if(topo[3]){
					_hists->Histogram::Fid_Fill(5,physics::get_theta(cz[par]),physics::get_phi(cx[par],cy[par]),par+1,6,0,_W,data->Branches::p(idx));
					_hists->Histogram::DT_Fill(5,par,_p[par], d[par], t[par], data->Branches::sc_r(0), data->Branches::sc_t(0),6,0,_W,h_sec[par]);
				}else if(topo[0] && i !=0){
					_hists->Histogram::Fid_Fill(5,physics::get_theta(cz[par]),physics::get_phi(cx[par],cy[par]),par+1,6,0,_W,data->Branches::p(idx));
					_hists->Histogram::DT_Fill(5,par,_p[par], d[par], t[par], data->Branches::sc_r(0), data->Branches::sc_t(0),6,0,_W,h_sec[par]);
				}else if(topo[2] && i!=1){
					_hists->Histogram::Fid_Fill(5,physics::get_theta(cz[par]),physics::get_phi(cx[par],cy[par]),par+1,6,0,_W,data->Branches::p(idx));
					_hists->Histogram::DT_Fill(5,par,_p[par], d[par], t[par], data->Branches::sc_r(0), data->Branches::sc_t(0),6,0,_W,h_sec[par]);	
				}else if(topo[3] && i!=2){
					_hists->Histogram::Fid_Fill(5,physics::get_theta(cz[par]),physics::get_phi(cx[par],cy[par]),par+1,6,0,_W,data->Branches::p(idx));
					_hists->Histogram::DT_Fill(5,par,_p[par], d[par], t[par], data->Branches::sc_r(0), data->Branches::sc_t(0),6,0,_W,h_sec[par]);
				}
			}
		}else{
			_hists->Histogram::Fid_Fill(5,theta[0],phi[0],0,10,1,_W,data->Branches::p(0));
			_hists->Histogram::SF_Fill(5,data->Branches::p(0),data->Branches::etot(0),10,1,_W,sector[0]);
			_hists->Histogram::CC_Fill(5,data->Branches::cc_sect(0),data->Branches::cc_segm(0),data->Branches::nphe(0),10,1);
			for(int i = 0; i<3; i++){
				switch(i){
					case 0:
					idx = pro_idx[0];
					par = 0; 
					break;
					case 1:
					idx = pip_idx[0];
					par = 1; 
					break;
					case 2:
					idx = pim_idx[0];
					par = 2;
					break;
				}
				if(topo[3]){
					_hists->Histogram::Fid_Fill(5,physics::get_theta(cz[par]),physics::get_phi(cx[par],cy[par]),par+1,6,1,_W,data->Branches::p(idx));
					_hists->Histogram::DT_Fill(5,par,_p[par], d[par], t[par], data->Branches::sc_r(0), data->Branches::sc_t(0),6,1,_W,h_sec[par]);
				}else if(topo[0] && i !=0){
					_hists->Histogram::Fid_Fill(5,physics::get_theta(cz[par]),physics::get_phi(cx[par],cy[par]),par+1,6,1,_W,data->Branches::p(idx));
					_hists->Histogram::DT_Fill(5,par,_p[par], d[par], t[par], data->Branches::sc_r(0), data->Branches::sc_t(0),6,1,_W,h_sec[par]);
				}else if(topo[2] && i!=1){
					_hists->Histogram::Fid_Fill(5,physics::get_theta(cz[par]),physics::get_phi(cx[par],cy[par]),par+1,6,1,_W,data->Branches::p(idx));
					_hists->Histogram::DT_Fill(5,par,_p[par], d[par], t[par], data->Branches::sc_r(0), data->Branches::sc_t(0),6,1,_W,h_sec[par]);	
				}else if(topo[3] && i!=2){
					_hists->Histogram::Fid_Fill(5,physics::get_theta(cz[par]),physics::get_phi(cx[par],cy[par]),par+1,6,1,_W,data->Branches::p(idx));
					_hists->Histogram::DT_Fill(5,par,_p[par], d[par], t[par], data->Branches::sc_r(0), data->Branches::sc_t(0),6,1,_W,h_sec[par]);
				}
			}
		}


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

int Event_Class::Get_ppip(){
	return ppip;
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
