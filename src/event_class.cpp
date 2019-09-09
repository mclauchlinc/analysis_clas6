#include "event_class.hpp"


Event_Class::Event_Class(std::shared_ptr<Branches> data, std::shared_ptr<Histogram> _hists, int run_type){ 
		//Pre ID Filling
		

		//electron ID
		

		//Hadron Loop


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
