#include "event.hpp"


Event::Event(std::shared_ptr<Branches> data, std::shared_ptr<Histogram> _hists, int run_type, int data_set){ 
	/*if(run_type %2 != 0){
		_beam = k_mu_e16;
	}else if(run_type%2 == 0){
		_beam = k_mu_e1f;
	}*/
	_beam = data->Branches::Par_4Vec(0)
	_elec= data->Branches::Par_4Vec(2);
	_gamma = _beam - _elec;
	_target = data->Branches::Par_4Vec(1);
	_pro= data->Branches::Par_4Vec(3);
	_pip= data->Branches::Par_4Vec(4);
	_pim= data->Branches::Par_4Vec(5);

	_fc_curr = data->Branches::fc_tot(); 

	//Get all relevant four vectors into the center of mass frame
	physics::COM_gp(_beam,_elec,_pro,_pip,_pim,_target,_gamma);


	_alpha1 = physics::alpha(0,_pro,_pip,_pim,_target);	//{}
	_alpha2 = physics::alpha(1,_pro,_pip,_pim,_target);	//{}
	_alpha3 = physics::alpha(2,_pro,_pip,_pim,_target);	//{}
	_theta1 = physics::get_theta(_pro); 	//proton
	_theta2 = physics::get_theta(_pip);	//pip
	_theta3 = physics::get_theta(_pim);	//pim
	_MMt1 = (_pro + _pip).Mag(); //MM proton/pip
	_MMt2 = (_pro + _pim).Mag();	//MM proton/pim
	_MMt3 = (_pim + _pip).Mag();	//MM Pip/pim

	_W = NAN;
	_Q2 = NAN; 

	_helicity = data->Branches::hel(); 
}

float Event::Get_px(int i){
	float _px = -99; 
	switch(i){
		case 0:
			_px = _elec[0];
		break;
		case 1:
			_px = _pro[0];
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

float Event::Get_py(int i){
	float _py = -99; 
	switch(i){
		case 0:
			_py = _elec[1];
		break;
		case 1:
			_py = _pro[1];
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

float Event::Get_pz(int i)
{
	float _pz = -99; 
	switch(i){
		case 0:
			_pz = _elec[2];
		break;
		case 1:
			_pz = _pro[2];
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

float Event::Get_p0(int i)
{
	float _p0 = -99; 
	switch(i){
		case 0:
			_p0 = _elec[3];
		break;
		case 1:
			_p0 = _pro[3];
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

float Event::Get_hel(){
	return _helicity; 
}

float Event::Get_top(){
	return _top; 
}

float Event::Get_pid(int i){
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

bool Event::is_valid(){
	return _valid; 
}

int Event::Get_ppip(int idx){
	return ppip[idx];
}

/*
void Fill_Tree(forest tree, int Event_n){
	tree.forest::fill_evnt(Event_n);
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
