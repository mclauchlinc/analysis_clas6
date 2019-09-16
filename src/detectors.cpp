#include "detectors.hpp"


int detect::cc_segment(int cc_segm){
	int seg = -1;

	if(cc_segm > 0 && cc_segm <200){
		seg = ((cc_segm)/10)-1;
		//std::cout<<"cc_seg: " <<cc_segm <<" -> " <<seg <<"             " <<std::endl;
	}
	if(cc_segm > 1000 && cc_segm <1200){
		seg = ((cc_segm-1000)/10)-1;
		//std::cout<<"cc_seg: " <<cc_segm <<" -> " <<seg<<"             " <<std::endl;
	}
	if(cc_segm > 2000 && cc_segm <2200){
		seg = ((cc_segm-2000)/10)-1;
		//std::cout<<"cc_seg: " <<cc_segm <<" -> " <<seg<<"             " <<std::endl;
	}
	return seg;
}

//Defining left, right, or coincidence hits in the CC
int detect::cc_lrc(int cc_segm){
	int po = -1;
	if(cc_segm > 0 && cc_segm <200){
		po = 0;//left
	}
	if(cc_segm > 1000 && cc_segm <1200){
		po = 1;//coincident
	}
	if(cc_segm > 2000 && cc_segm <2200){
		po = 2;//right
	}
	
	return po;
}