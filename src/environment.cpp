#include "environment.hpp"

Environment::Environment(){

}
void Environment::env_data_set(int data_){
	_data_set = data_; 
}
void Environment::env_dc_hit(bool dc_hit_){//Was a hit on the DC required?
	_dc_hit = dc_hit_;
}
void Environment::env_cc_hit(bool cc_hit_) //Was a hit on the CC required for electrons?
{
	_cc_hit = cc_hit_;
}
void Environment::env_ec_hit(bool ec_hit_)//Was a hit on the EC required for electrons?
{
	_ec_hit = ec_hit_;
}
void Environment::env_sc_hit(bool sc_hit_)//Was a hit on the SC required?
{
	_sc_hit = sc_hit_;
}
void Environment::env_stat_cut(bool stat_cut_)//Was a nonzero stat required?
{
	_stat_cut = stat_cut_;
}
void Environment::env_dc_stat_cut(bool dc_stat_cut_) //Was a nonzero dc_stat required?
{
	_dc_stat_cut = dc_stat_cut_;
}
void Environment::env_sim(bool sim_) //Was this simulation data?
{
	_sim = sim_;
}
void Environment::env_top(int top_, bool yes_top_)//Which topologies are included {pmiss,pipmiss,pimmiss,zeromiss,combined}
{
	_top[top_] = yes_top_;
}
void Environment::env_p_corr(bool p_corr_)//Were momentum corrections performed?
{
	_p_corr = p_corr_;
}
/*void Environment::env_fid(int species_, bool fid_)//Were fiducial cuts performed {e,p,pip,pim}
{
	_fid[species_] = fid_;
}
void Environment::env_dt(int had_, bool dt_)//Were delta t cuts performed {p,pip,pim}
{
	_dt[had_] = dt_;
}*/
void Environment::env_pvpip(bool pvpip_)//Were procedures allowed for dual identified positive hadrons?
{
	_pvpip = pvpip_;
}
void Environment::env_npart(int npar_)//Number of particles allowed in an event 
{
	_npart = npar_;
}
void Environment::env_Qmin(float Qmin_)//Lower range of Q squared allowed
{
	_Qmin = Qmin_;
}
void Environment::env_Qmax(float Qmax_)//Upper range of Q squared allowed
{
	_Qmax = Qmax_; 
}
void Environment::env_Wmin(float Wmin_)//Lower limit on W allowed
{
	_Wmin = Wmin_;
}
void Environment::env_Wmax(float Wmax_)//Upper limit on W allowed 
{
	_Wmax = Wmax_;
}
void Environment::env_golden(bool gold_)//Was this a golden run?
{
	_golden = gold_;
}
void Environment::env_frac(float frac_)//Fraction of the run analyzed
{
	_frac = frac_; 
}
void Environment::env_num_file(int num_file_)//Number of files analyzed
{
	_num_file = _num_file + num_file_;
}
void Environment::env_skim(bool skim_)//Were the files initially skimmed
{
	_skim = skim_; 
}
void Environment::env_skim_type(int skim_type_)//What type of skim was performed? See documentation 
{
	_skim_type = skim_type_; 
}
void Environment::env_efficiency(bool eff_)//Efficiencies of detector used
{
	_efficiency = eff_; 
}
void Environment::env_empty(bool empty_)//Is this an empty target run?
{
	_empty = empty_;
}
void Environment::env_COM(bool COM_)//Are the four vectors already put voido the center of mass?
{
	_COM = COM_;
}
//For cut versions look voido the notes on the analysis
void Environment::env_fid_ver(int species_, int fid_ver_)//Version of the fiducial cut. 
{
	_fid_ver[species_] = fid_ver_; 
}
void Environment::env_dt_ver(int species_, int dt_ver_)//Version of delta t cut. 
{
	_dt_ver[species_] = dt_ver_;
}
void Environment::env_MM_ver(int species_, int MM_ver_)//Version of Missing Mass cut. 
{
	_MM_ver[species_] = MM_ver_;
}
void Environment::env_cc_ver(int ver_)//Version of Min-CC cut. 
{
	_cc_ver = ver_;
}
void Environment::env_ec_ver(int ver_)//Version of Min EC cut. 
{
	_ec_ver = ver_; 
}
void Environment::env_sf_ver(int ver_)//Version of Sampling Fraction cut. 
{
	_sf_ver = ver_; 
}
//Eid Cuts applied
void Environment::env_eid_fid(bool e_fid_)//Did eid include a fiducial cut?
{
	_eid_fid = e_fid_;
}
void Environment::env_eid_ec(bool e_ec_)//Did eid include a min ec cut?
{
	_eid_ec = e_ec_;
}
void Environment::env_eid_cc(bool e_cc_)//Did eid include a min cc cut?
{
	_eid_cc = e_cc_; 
}
void Environment::env_eid_sf(bool e_sf_)//Did eid include a sampling fraction cut?
{
	_eid_sf = e_sf_; 
}
void Environment::env_eid_dt(bool e_dt_){
	_eid_dt = e_dt_; 
}
//Hadron cuts applied {p,pip,pim}
void Environment::env_hid_fid(int had_, bool fid_)//Did hid include a fiducial cut?
{
	_hid_fid[had_] = fid_;
}
void Environment::env_hid_dt(int had_, bool dt_)//Did hid include a delta t cut?
{
	_hid_dt[had_] = dt_;
}
void Environment::env_hid_beta(int had_, bool beta_)//Did hid include a beta cut?
{
	_hid_beta[had_] = beta_;
}
void Environment::env_hid_e(bool hid_e_)//Did hid cut out a band of electrons for pim ID?
{
	_hid_e = hid_e_;
}
//Plot info
void Environment::env_fid_plot(int species_, bool fid_plot_)//Were fiducial cuts performed {e,p,pip,pim}
{
	_fid_plot[species_] = fid_plot_;
}
void Environment::env_dt_plot(int had_, bool dt_plot_) //Were delta t cuts performed {p,pip,pim}
{
	_dt_plot[had_] = dt_plot_;
}
void Environment::env_cc_plot(bool cc_plot_) //Were Min CC plots made? 
{
	_cc_plot = cc_plot_;
}
void Environment::env_ec_plot(bool ec_plot_) //were Min EC plots made?
{
	_ec_plot = ec_plot_;
}
void Environment::env_sf_plot(bool sf_plot_) //Were Sampling Fraction plots made?
{
	_sf_plot = sf_plot_;
}
void Environment::env_MM_plot(int species_, bool MM_plot_) //Were missing mass plots made? {pmiss,pipmiss,pimmiss,zmiss}
{
	_MM_plot[species_] = MM_plot_;
}
void Environment::env_WQ2_plot(bool WQ2_plot_) //Were W Qsquared plots made?
{
	_WQ2_plot = WQ2_plot_;
}
void Environment::env_p_dep_plot(bool p_dep_plot_) //Is there momentum dependence in the relevant plots?
{
	_p_dep_plot = p_dep_plot_;
}
void Environment::env_W_dep_plot(bool W_dep_plot_) //Is there W dependence in the relevant plots? 
{
	_W_dep_plot = W_dep_plot_;
}
void Environment::env_Friend_plot(bool friend_plot){
	_friend_plot = friend_plot;
}

void Environment::env_fitting(bool fitting_, int fit_){
	_fitting = fitting_;
	_fit = fit_;
}
//Switch to the Was 

int Environment::was_data_set(){
	return _data_set; 
}

bool Environment::was_dc_hit()//Was a hit on the DC required?
{
	return _dc_hit; 
}
bool Environment::was_cc_hit() //Was a hit on the CC required for electrons?
{
	return _cc_hit; 
}
bool Environment::was_ec_hit()//Was a hit on the EC required for electrons?
{
	return _ec_hit; 
}
bool Environment::was_sc_hit()//Was a hit on the SC required?
{
	return _sc_hit; 
}
bool Environment::was_stat_cut()//Was a nonzero stat required?
{
	return _stat_cut; 
}
bool Environment::was_dc_stat_cut() //Was a nonzero dc_stat required?
{
	return _dc_stat_cut; 
}
bool Environment::was_sim() //Was this simulation data?
{
	return _sim; 
}
bool Environment::was_top(int species_)//Which topologies are included {pmiss,pipmiss,pimmiss,zeromiss,combined}
{
	return _top[species_];
}
bool Environment::was_p_corr()//Were momentum corrections performed?
{
	return _p_corr; 
}
/*bool Environment::was_fid(int species_)//Were fiducial cuts performed {e,p,pip,pim}
{
	return _fid[species_];
}
bool Environment::was_dt(int species_)//Were delta t cuts performed {p,pip,pim}
{
	return _dt[species_];
}*/
bool Environment::was_pvpip()//Were procedures allowed for dual identified positive hadrons?
{
	return _pvpip; 
}
int Environment::was_npart()//Number of particles allowed in an event 
{
	return _npart; 
}
float Environment::was_Qmin()//Lower range of Q squared allowed
{
	return _Qmin; 
}
float Environment::was_Qmax()//Upper range of Q squared allowed
{
	return _Qmax; 
}
float Environment::was_Wmin()//Lower limit on W allowed
{
	return _Wmin; 
}
float Environment::was_Wmax()//Upper limit on W allowed 
{
	return _Wmax; 
}
bool Environment::was_golden()//Was this a golden run?
{
	return _golden; 
}
float Environment::was_frac()//Fraction of the run analyzed
{
	return _frac; 
}
int Environment::was_num_file()//Number of files analyzed
{
	return _num_file; 
}
bool Environment::was_skim()//Were the files initially skimmed
{
	return _skim; 
}
int Environment::was_skim_type()//What type of skim was performed? See documentation 
{
	return _skim_type; 
}
bool Environment::was_efficiency()//Efficiencies of detector used
{	
	return _efficiency; 
}
bool Environment::was_empty()//Is this an empty target run?
{
	return _empty; 
}
bool Environment::was_COM()//Are the four vectors already put into the center of mass?
{
	return _COM; 
}
//For cut versions look into the notes on the analysis
int Environment::was_fid_ver(int species_)//Version of the fiducial cut. 
{
	return _fid_ver[species_];
}
int Environment::was_dt_ver(int species_)//Version of delta t cut. 
{
	return _dt_ver[species_];
}
int Environment::was_MM_ver(int species_)//Version of Missing Mass cut. 
{
	return _MM_ver[species_];
}
int Environment::was_cc_ver()//Version of Min-CC cut. 
{
	return _cc_ver;
}
int Environment::was_ec_ver()//Version of Min EC cut. 
{
	return _ec_ver; 
}
int Environment::was_sf_ver()//Version of Sampling Fraction cut. 
{
	return _sf_ver; 
}
//Eid Cuts applied
bool Environment::was_eid_fid()//Did eid include a fiducial cut?
{
	return _eid_fid; 
}
bool Environment::was_eid_ec()//Did eid include a min ec cut?
{
	return _eid_ec; 
}
bool Environment::was_eid_cc()//Did eid include a min cc cut?
{
	return _eid_cc; 
}
bool Environment::was_eid_sf()//Did eid include a sampling fraction cut?
{
	return _eid_sf; 
}

bool Environment::was_eid_dt(){
	return _eid_dt;
}
//Hadron cuts applied {p,pip,pim}
bool Environment::was_hid_fid(int species_)//Did hid include a fiducial cut?
{
	return _hid_fid; 
}
bool Environment::was_hid_dt(int species_)//Did hid include a delta t cut?
{
	return _hid_dt; 
}
bool Environment::was_hid_beta(int species_)//Did hid include a beta cut?
{
	return _hid_beta; 
}
bool Environment::was_hid_e()//Did hid cut out a band of electrons for pim ID?
{
	return _hid_e; 
}

bool Environment::was_fid_plot(int species_)//Were fiducial cuts performed {e,p,pip,pim}
{
	return _fid_plot[species_];
}
bool Environment::was_dt_plot(int had_) //Were delta t cuts performed {p,pip,pim}
{
	return _dt_plot[had_];
}
bool Environment::was_cc_plot() //Were Min CC plots made? 
{
	return _cc_plot;
}
bool Environment::was_ec_plot() //were Min EC plots made?
{
	return _ec_plot; 
}
bool Environment::was_sf_plot() //Were Sampling Fraction plots made?
{
	return _sf_plot; 
}
bool Environment::was_MM_plot(int species_) //Were missing mass plots made? {pmiss,pipmiss,pimmiss,zmiss}
{
	return _MM_plot[species_];
}
bool Environment::was_WQ2_plot() //Were W Qsquared plots made?
{
	return _WQ2_plot; 
}
bool Environment::was_p_dep_plot() //Is there momentum dependence in the relevant plots?
{
	return _p_dep_plot;
}
bool Environment::was_W_dep_plot() //Is there W dependence in the relevant plots? 
{
	return _W_dep_plot;
}
bool Environment::was_Friend_plot()
{
	return _friend_plot;
}

bool Environment::was_fitting(){
	return _fitting;
}

int Environment::was_fit_type(){
	return _fit;
}
/*
void copy_envi(Environment *cop_envi){
	_data_set = cop_envi->Environment::was_data_set(); //Which data set? {1,2}->{e16,e1f}
  	_dc_hit= cop_envi->Environment::was_dc_hit(); //Was a hit on the DC required?
  	_cc_hit= cop_envi->Environment::was_cc_hit();  //Was a hit on the CC required for electrons?
  	_ec_hit = cop_envi->Environment::was_ec_hit(); //Was a hit on the EC required for electrons?
  	_sc_hit = cop_envi->Environment::was_sc_hit(); //Was a hit on the SC required?
  	_stat_cut = cop_envi->Environment::was_stat_cut(); //Was a nonzero stat required?
  	_dc_stat_cut = cop_envi->Environment::was_dc_stat_cut();  //Was a nonzero dc_stat required?
  	_sim = cop_envi->Environment::was_sim();  //Was this simulation data?
  	for(int i = 0; i< 5; i++){
  		_top[i] = cop_envi->Environment::was_top(i);
  	}
  	for(int j = 0; j< 4; j++){
  		_fid[j] = cop_envi->Environment::was_fid(i);
  		_fid_ver[j] = cop_envi->Environment::was_fid_ver(j);
  		_MM_ver[j] = cop_envi->Environment::was_MM_ver(j);
  	}
  	for(int k = 0; k<3; k++){
  		_dt[k] = cop_envi->Environment::was_dt(k);
  		_dt_ver[k] = cop_envi->Environment::was_dt_ver(k);
  		_hid_fid[k] = cop_envi->Environment::was_hid_fid(k);
  		_hid_dt[k] = cop_envi->Environment::was_hid_dt(k);
  		_hid_beta[k] = cop_envi->Environment::was_hid_beta(k);
  	}
  	_p_corr = cop_envi->Environment::was_p_corr(); //Were momentum corrections performed?
  	_pvpip = cop_envi->Environment::was_pvpip(); //Were procedures allowed for dual identified positive hadrons?
  	_npart = cop_envi->Environment::was_npart(); //Number of particles allowed in an event 
  	_Qmin = cop_envi->Environment::was_Qmin(); //Lower range of Q squared allowed
  	_Qmax = cop_envi->Environment::was_Qmax(); //Upper range of Q squared allowed
  	_Wmin = cop_envi->Environment::was_Wmin(); //Lower limit on W allowed
  	_Wmax = cop_envi->Environment::was_Wmax(); //Upper limit on W allowed 
  	_golden = cop_envi->Environment::was_golden();//Was this a golden run?
  	_frac = cop_envi->Environment::was_frac(); //Fraction of the run analyzed
  	_num_file = cop_envi->Environment::was_num_file(); //Number of files analyzed
  	_skim = cop_envi->Environment::was_skim();//Were the files initially skimmed
  	_skim_type = cop_envi->Environment::was_skim_type(); //What type of skim was performed? See documentation 
  	_efficiency = cop_envi->Environment::was_efficiency(); //Efficiencies of detector used
  	_empty = cop_envi->Environment::was_empty();//Is this an empty target run?
  	_COM = cop_envi->Environment::was_COM(); //Are the four vectors already put into the center of mass?
  //For cut versions look into the notes on the analysis 
  	_cc_ver = cop_envi->Environment::was_cc_ver(); //Version of Min-CC cut. 
  	_ec_ver = cop_envi->Environment::was_cc_ver(); //Version of Min EC cut. 
  	_sf_ver = cop_envi->Environment::was_cc_ver(); //Version of Sampling Fraction cut. 
  //Eid Cuts applied
  	_eid_fid = cop_envi->Environment::was_eid_fid(); //Did eid include a fiducial cut?
  	_eid_ec = cop_envi->Environment::was_eid_ec(); //Did eid include a min ec cut?
  	_eid_cc = cop_envi->Environment::was_eid_cc(); //Did eid include a min cc cut?
  	_eid_sf = cop_envi->Environment::was_eid_sf(); //Did eid include a sampling fraction cut?
  //Hadron cuts applied {p,pip,pim}
  	_hid_e = cop_envi->Environment::was_hid_e(); //Did hid cut out a band of electrons for pim ID?

}
*/