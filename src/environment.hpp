#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include "constants.hpp"


//Create a data chain for a given data set
//Optional number of files to load

class Environment {
 private:
  int _data_set = 0; //Which data set? {1,2}->{e16,e1f}
  bool _dc_hit= false; //Was a hit on the DC required?
  bool _cc_hit= false;  //Was a hit on the CC required for electrons?
  bool _ec_hit = false; //Was a hit on the EC required for electrons?
  bool _sc_hit = false; //Was a hit on the SC required?
  bool _stat_cut = false; //Was a nonzero stat required?
  bool _dc_stat_cut = false;  //Was a nonzero dc_stat required?
  bool _sim = false;  //Was this simulation data?
  bool _top[5] = {false,false,false,false,false}; //Which topologies are included {pmiss,pipmiss,pimmiss,zeromiss,combined}
  bool _p_corr = false; //Were momentum corrections performed?
  bool _pvpip = false; //Were procedures allowed for dual identified positive hadrons?
  int _npart = 0; //Number of particles allowed in an event 
  float _Qmin = NAN; //Lower range of Q squared allowed
  float _Qmax = NAN; //Upper range of Q squared allowed
  float _Wmin = NAN; //Lower limit on W allowed
  float _Wmax = NAN; //Upper limit on W allowed 
  bool _golden = false;//Was this a golden run?
  float _frac = NAN; //Fraction of the run analyzed
  int _num_file = 0; //Number of files analyzed
  bool _skim = false;//Were the files initially skimmed
  int _skim_type = -1; //What type of skim was performed? See documentation 
  bool _efficiency = false; //Efficiencies of detector used
  bool _empty = true;//Is this an empty target run?
  bool _COM = false; //Are the four vectors already put into the center of mass?
  //For cut versions look into the notes on the analysis
  int _fid_ver[4] = {-1,-1,-1,-1};//Version of the fiducial cut. 
  int _dt_ver[3] = {-1,-1,-1};//Version of delta t cut. 
  int _MM_ver[4] = {-1,-1,-1,-1};//Version of Missing Mass cut. 
  int _cc_ver = -1; //Version of Min-CC cut. 
  int _ec_ver = -1; //Version of Min EC cut. 
  int _sf_ver = -1; //Version of Sampling Fraction cut. 
  //Eid Cuts applied
  bool _eid_fid = false; //Did eid include a fiducial cut?
  bool _eid_ec = false; //Did eid include a min ec cut?
  bool _eid_cc = false; //Did eid include a min cc cut?
  bool _eid_sf = false; //Did eid include a sampling fraction cut?
  //Hadron cuts applied {p,pip,pim}
  bool _hid_fid[3] = {false,false,false};//Did hid include a fiducial cut?
  bool _hid_dt[3] = {false,false,false};//Did hid include a delta t cut?
  bool _hid_beta[3] = {false,false,false};//Did hid include a beta cut?
  bool _hid_e = false; //Did hid cut out a band of electrons for pim ID?
  //Plots being made
  bool _fid_plot[4] = {false,false,false,false};//Were fiducial cuts performed {e,p,pip,pim}
  bool _dt_plot[3] = {false,false,false}; //Were delta t cuts performed {p,pip,pim}
  bool _cc_plot = false; //Were Min CC plots made? 
  bool _ec_plot = false; //were Min EC plots made?
  bool _sf_plot = false; //Were Sampling Fraction plots made?
  bool _MM_plot[4] = {false, false, false, false}; //Were missing mass plots made? {pmiss,pipmiss,pimmiss,zmiss}
  bool _WQ2_plot = false; //Were W Qsquared plots made?
  bool _p_dep_plot = false; //Is there momentum dependence in the relevant plots?
  bool _W_dep_plot = false; //Is there W dependence in the relevant plots? 
  

 public:
  Environment();
  ~Environment(){};
  void env_data_set(int data_); //Which data set? {1,2}->{e16,e1f}
  void env_dc_hit(bool dc_hit_); //Was a hit on the DC required?
  void env_cc_hit(bool cc_hit_);  //Was a hit on the CC required for electrons?
  void env_ec_hit(bool ec_hit_); //Was a hit on the EC required for electrons?
  void env_sc_hit(bool sc_hit_); //Was a hit on the SC required?
  void env_stat_cut(bool stat_cut_); //Was a nonzero stat required?
  void env_dc_stat_cut(bool dc_stat_cut_);  //Was a nonzero dc_stat required?
  void env_sim(bool sim_);  //Was this simulation data?
  void env_top(int top_, bool yes_top_); //Which topologies are included {pmiss,pipmiss,pimmiss,zeromiss,combined}
  void env_p_corr(bool p_corr_); //Were momentum corrections performed?
  void env_pvpip(bool pvpip_); //Were procedures allowed for dual identified positive hadrons?
  void env_npart(int npar_); //Number of particles allowed in an event 
  void env_Qmin(float Qmin_); //Lower range of Q squared allowed
  void env_Qmax(float Qmax_); //Upper range of Q squared allowed
  void env_Wmin(float Wmin_); //Lower limit on W allowed
  void env_Wmax(float Wmax_); //Upper limit on W allowed 
  void env_golden(bool gold_);//Was this a golden run?
  void env_frac(float frac_); //Fraction of the run analyzed
  void env_num_file(int num_file_); //Number of files analyzed
  void env_skim(bool skim_);//Were the files initially skimmed
  void env_skim_type(int skim_type_); //What type of skim was performed? See documentation 
  void env_efficiency(bool eff_); //Efficiencies of detector used
  void env_empty(bool empty_);//Is this an empty target run?
  void env_COM(bool COM_); //Are the four vectors already put voido the center of mass?
  //For cut versions look voido the notes on the analysis
  void env_fid_ver(int species_, int fid_ver_);//Version of the fiducial cut. 
  void env_dt_ver(int species_, int dt_ver_);//Version of delta t cut. 
  void env_MM_ver(int species_, int MM_ver_);//Version of Missing Mass cut. 
  void env_cc_ver(int ver_); //Version of Min-CC cut. 
  void env_ec_ver(int ver_); //Version of Min EC cut. 
  void env_sf_ver(int ver_); //Version of Sampling Fraction cut. 
  //Eid Cuts applied
  void env_eid_fid(bool e_fid_); //Did eid include a fiducial cut?
  void env_eid_ec(bool e_ec_); //Did eid include a min ec cut?
  void env_eid_cc(bool e_cc_); //Did eid include a min cc cut?
  void env_eid_sf(bool e_sf_); //Did eid include a sampling fraction cut?
  //Hadron cuts applied {p,pip,pim}
  void env_hid_fid(int had_, bool fid_);//Did hid include a fiducial cut?
  void env_hid_dt(int had_, bool dt_);//Did hid include a delta t cut?
  void env_hid_beta(int had_, bool beta_);//Did hid include a beta cut?
  void env_hid_e(bool hid_e_); //Did hid cut out a band of electrons for pim ID?
  //Plot info
  void env_fid_plot(int species_, bool fid_plot_);//Were fiducial cuts performed {e,p,pip,pim}
  void env_dt_plot(int had_, bool dt_plot_); //Were delta t cuts performed {p,pip,pim}
  void env_cc_plot(bool cc_plot_); //Were Min CC plots made? 
  void env_ec_plot(bool ec_plot_); //were Min EC plots made?
  void env_sf_plot(bool sf_plot_); //Were Sampling Fraction plots made?
  void env_MM_plot(int species_, bool MM_plot_); //Were missing mass plots made? {pmiss,pipmiss,pimmiss,zmiss}
  void env_WQ2_plot(bool WQ2_plot_); //Were W Qsquared plots made?
  void env_p_dep_plot(bool p_dep_plot_); //Is there momentum dependence in the relevant plots?
  void env_W_dep_plot(bool W_dep_plot_); //Is there W dependence in the relevant plots? 

  int was_data_set();
  bool was_dc_hit(); //Was a hit on the DC required?
  bool was_cc_hit();  //Was a hit on the CC required for electrons?
  bool was_ec_hit(); //Was a hit on the EC required for electrons?
  bool was_sc_hit(); //Was a hit on the SC required?
  bool was_stat_cut(); //Was a nonzero stat required?
  bool was_dc_stat_cut();  //Was a nonzero dc_stat required?
  bool was_sim();  //Was this simulation data?
  bool was_top(int species_); //Which topologies are included {pmiss,pipmiss,pimmiss,zeromiss,combined}
  bool was_p_corr(); //Were momentum corrections performed?
  bool was_pvpip(); //Were procedures allowed for dual identified positive hadrons?
  int was_npart(); //Number of particles allowed in an event 
  float was_Qmin(); //Lower range of Q squared allowed
  float was_Qmax(); //Upper range of Q squared allowed
  float was_Wmin(); //Lower limit on W allowed
  float was_Wmax(); //Upper limit on W allowed 
  bool was_golden();//Was this a golden run?
  float was_frac(); //Fraction of the run analyzed
  int was_num_file(); //Number of files analyzed
  bool was_skim();//Were the files initially skimmed
  int was_skim_type(); //What type of skim was performed? See documentation 
  bool was_efficiency(); //Efficiencies of detector used
  bool was_empty();//Is this an empty target run?
  bool was_COM(); //Are the four vectors already put into the center of mass?
  //For cut versions look into the notes on the analysis
  int was_fid_ver(int species_);//Version of the fiducial cut. 
  int was_dt_ver(int species_);//Version of delta t cut. 
  int was_MM_ver(int species_);//Version of Missing Mass cut. 
  int was_cc_ver(); //Version of Min-CC cut. 
  int was_ec_ver(); //Version of Min EC cut. 
  int was_sf_ver(); //Version of Sampling Fraction cut. 
  //Eid Cuts applied
  bool was_eid_fid(); //Did eid include a fiducial cut?
  bool was_eid_ec(); //Did eid include a min ec cut?
  bool was_eid_cc(); //Did eid include a min cc cut?
  bool was_eid_sf(); //Did eid include a sampling fraction cut?
  //Hadron cuts applied {p,pip,pim}
  bool was_hid_fid(int species_);//Did hid include a fiducial cut?
  bool was_hid_dt(int species_);//Did hid include a delta t cut?
  bool was_hid_beta(int species_);//Did hid include a beta cut?
  bool was_hid_e(); //Did hid cut out a band of electrons for pim ID?
  //Plot creation and filling
  bool was_fid_plot(int species_);//Were fiducial cuts performed {e,p,pip,pim}
  bool was_dt_plot(int had_); //Were delta t cuts performed {p,pip,pim}
  bool was_cc_plot(); //Were Min CC plots made? 
  bool was_ec_plot(); //were Min EC plots made?
  bool was_sf_plot(); //Were Sampling Fraction plots made?
  bool was_MM_plot(int species_); //Were missing mass plots made? {pmiss,pipmiss,pimmiss,zmiss}
  bool was_WQ2_plot(); //Were W Qsquared plots made?
  bool was_p_dep_plot(); //Is there momentum dependence in the relevant plots?
  bool was_W_dep_plot(); //Is there W dependence in the relevant plots? 

  //void copy_envi(Environment cop_envi);
};

#endif
