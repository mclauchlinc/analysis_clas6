#include "event_class.hpp"

Event::Event(std::shared_ptr<Branches> data, std::shared_ptr<Histogram> _hists, int run_type){ 
		//Pre ID Filling
		_W = 10; 

		//electron ID
		

		//Hadron Loop


		//Event Selection


		//Variable Determination


}

float Event::Get_px(int i){
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

float Event::Get_py(int i){
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

float Event::Get_pz(int i)
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

float Event::Get_p0(int i)
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
	}
	return _pid; 
}
