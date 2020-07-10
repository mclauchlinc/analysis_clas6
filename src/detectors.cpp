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

TVector3 detect::cc_sector_center(TVector3 v_){
	float _phi_ = physics::get_phi(v_[0],v_[1]);
	int _sector_ = physics::get_sector(_phi_);
	float _x_ = NAN;
	float _y_ = NAN;
	float _z_ = v_[3];
	_x_ = physics::CCX_Rotate(v_[0],v_[1],_sector_);
	_y_ = physics::CCY_Rotate(v_[0],v_[1],_sector_);
	TVector3 _v_ = {_x_,_y_,_z_}; 
	return _v_;
}
TVector3 detect::cc_sector_center(float x_, float y_, float z_){
	float _phi_ = physics::get_phi(x_,y_);
	int _sector_ = physics::get_sector(_phi_);
	float _x_ = NAN;
	float _y_ = NAN;
	float _z_ = z_;
	_x_ = physics::CCX_Rotate(x_,y_,_sector_);
	_y_ = physics::CCY_Rotate(x_,y_,_sector_);
	TVector3 _v_ = {_x_,_y_,_z_}; 
	return _v_;
}

float detect::cc_theta(float cx_sc_, float cy_sc_, float cz_sc_, float x_sc_, float y_sc_, float z_sc_){
	TVector3 n(cx_sc_, cy_sc_, cz_sc_);
	TVector3 p0(x_sc_, y_sc_, z_sc_);
	TVector3 s(Acc,Bcc,Ccc);
	//TVector3 s(physics::X_Rotate(Acc,Bcc,physics::get_sector(physics::get_phi(cx_sc_,cy_sc_))),physics::Y_Rotate(Acc,Bcc,physics::get_sector(physics::get_phi(cx_sc_,cy_sc_))),Ccc);//Get Projective plane for particular secgtor we're working in
	float tmag = TMath::Abs((physics::Dot_Product(s,p0)+Dcc)/physics::Dot_Product(s,n)); 
	TVector3 t = tmag*n;
	TVector3 _p_ = p0 + t; 
	return TMath::ACos(_p_[2]/physics::Vec3_Mag(_p_))*(180.0/TMath::Pi());
}

float detect::cc_theta(std::shared_ptr<Branches> data_, int idx_){
	return detect::cc_theta(data_->dc_cxsc(idx_), data_->dc_cysc(idx_), data_->dc_czsc(idx_), data_->dc_xsc(idx_), data_->dc_ysc(idx_), data_->dc_zsc(idx_));
}

float detect::cc_phi(float cx_sc_, float cy_sc_, float cz_sc_, float x_sc_, float y_sc_, float z_sc_){
	TVector3 n(cx_sc_, cy_sc_, cz_sc_);
	TVector3 p0(x_sc_, y_sc_, z_sc_);
	TVector3 s(Acc,Bcc,Ccc);
	//TVector3 s(physics::X_Rotate(Acc,Bcc,physics::get_sector(physics::get_phi(cx_sc_,cy_sc_))),physics::Y_Rotate(Acc,Bcc,physics::get_sector(physics::get_phi(cx_sc_,cy_sc_))),Ccc);//Get Projective plane for particular secgtor we're working in
	float tmag = TMath::Abs((physics::Dot_Product(s,p0)+Dcc)/physics::Dot_Product(s,n)); 
	TVector3 t = tmag*n;
	TVector3 _p_ = p0 + t; 
	return TMath::ATan2(_p_[1],_p_[0])*(180.0/TMath::Pi());
}

float detect::cc_phi(std::shared_ptr<Branches> data_, int idx_){
	return detect::cc_phi(data_->dc_cxsc(idx_), data_->dc_cysc(idx_), data_->dc_czsc(idx_), data_->dc_xsc(idx_), data_->dc_ysc(idx_), data_->dc_zsc(idx_));
}

float detect::cc_x(float cx_sc_, float cy_sc_, float cz_sc_, float x_sc_, float y_sc_, float z_sc_, int sec_){
	TVector3 n(cx_sc_, cy_sc_, cz_sc_);
	TVector3 p0(x_sc_, y_sc_, z_sc_);
	TVector3 s(Acc,Bcc,Ccc);
	//TVector3 s(physics::X_Rotate(Acc,Bcc,physics::get_sector(physics::get_phi(cx_sc_,cy_sc_))),physics::Y_Rotate(Acc,Bcc,physics::get_sector(physics::get_phi(cx_sc_,cy_sc_))),Ccc);//Get Projective plane for particular secgtor we're working in
	float tmag = TMath::Abs((physics::Dot_Product(s,p0)+Dcc)/physics::Dot_Product(s,n)); 
	TVector3 t = tmag*n;
	TVector3 _p_ = p0 + t; 
	return physics::X_Rotate(_p_[0],_p_[1],sec_);
}

float detect::cc_y(float cx_sc_, float cy_sc_, float cz_sc_, float x_sc_, float y_sc_, float z_sc_, int sec_){
	TVector3 n(cx_sc_, cy_sc_, cz_sc_);
	TVector3 p0(x_sc_, y_sc_, z_sc_);
	TVector3 s(Acc,Bcc,Ccc);
	//TVector3 s(physics::X_Rotate(Acc,Bcc,physics::get_sector(physics::get_phi(cx_sc_,cy_sc_))),physics::Y_Rotate(Acc,Bcc,physics::get_sector(physics::get_phi(cx_sc_,cy_sc_))),Ccc);//Get Projective plane for particular secgtor we're working in
	float tmag = TMath::Abs((physics::Dot_Product(s,p0)+Dcc)/physics::Dot_Product(s,n)); 
	TVector3 t = tmag*n;
	TVector3 _p_ = p0 + t; 
	return physics::Y_Rotate(_p_[0],_p_[1],sec_);
}

float detect::cc_x(std::shared_ptr<Branches> data_, int idx_){
	return detect::cc_x(data_->dc_cxsc(idx_), data_->dc_cysc(idx_), data_->dc_czsc(idx_), data_->dc_xsc(idx_), data_->dc_ysc(idx_), data_->dc_zsc(idx_), data_->sc_sect(idx_));
}

float detect::cc_y(std::shared_ptr<Branches> data_, int idx_){
	return detect::cc_y(data_->dc_cxsc(idx_), data_->dc_cysc(idx_), data_->dc_czsc(idx_), data_->dc_xsc(idx_), data_->dc_ysc(idx_), data_->dc_zsc(idx_), data_->sc_sect(idx_)); 
}



