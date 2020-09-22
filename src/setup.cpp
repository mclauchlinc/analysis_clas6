#include "setup.hpp"

void Setup::set_envi(std::shared_ptr<Environment> setup, int run_type, int fit){
	if(fit >= 0){
		setup->Environment::env_fitting(true,fit);
	}else{
		setup->Environment::env_fitting(false,fit);
	}
	switch(run_type%2){
		case 1:
			setup->Environment::env_data_set(1); //Which data set? {1,2}->{e16,e1f}
		break;
		case 0:
			setup->Environment::env_data_set(2); //Which data set? {1,2}->{e16,e1f}
		break;
		default:
			setup->Environment::env_data_set(1); //Which data set? {1,2}->{e16,e1f}
		break;
	}
	if(run_type >= 3){
		setup->Environment::env_sim(true);
	}
	if(run_type < 5 || run_type == 7){//Everything running on normal banks
		setup->Environment::env_dc_hit(true); //Was a hit on the DC required?
		if(run_type < 3){//Only for experiment do we work with the CC
			setup->Environment::env_cc_hit(true);  //Was a hit on the CC required for electrons?
			setup->Environment::env_cc_ver(1); //Version of Min-CC cut. 
			setup->Environment::env_eid_cc(true); //Did eid include a min cc cut?
			setup->Environment::env_cc_plot(true); //Were Min CC plots made?
			setup->Environment::env_sim(false); 
		}
		setup->Environment::env_ec_hit(true); //Was a hit on the EC required for electrons?
		setup->Environment::env_sc_hit(true); //Was a hit on the SC required?
		setup->Environment::env_stat_cut(false); //Was a nonzero stat required?
		setup->Environment::env_dc_stat_cut(false);  //Was a nonzero dc_stat required?
		//setup->Environment::env_sim(false);  //Was this simulation data?
		//Which topologies are included {pmiss,pipmiss,pimmiss,zeromiss,combined}
		setup->Environment::env_top(0, true);//Proton Missing
		setup->Environment::env_top(1, true);//Pip Missing
		setup->Environment::env_top(2, true);//Pim Missing
		setup->Environment::env_top(3, true); //Zero MissingWhich topologies are included {pmiss,pipmiss,pimmiss,zeromiss,combined}
		setup->Environment::env_p_corr(false); //Were momentum corrections performed?
		
		setup->Environment::env_pvpip(true); //Were procedures allowed for dual identified positive hadrons?
		setup->Environment::env_npart(-1); //Number of particles allowed in an event (-1 means no limits)
		setup->Environment::env_Qmin(Q2minAna); //Lower range of Q squared allowed
		setup->Environment::env_Qmax(Q2maxAna); //Upper range of Q squared allowed
		setup->Environment::env_Wmin(WminAna); //Lower limit on W allowed
		setup->Environment::env_Wmax(WmaxAna); //Upper limit on W allowed 
		setup->Environment::env_golden(false);//Was this a golden run?
		setup->Environment::env_frac(1.0); //Fraction of the run analyzed
		setup->Environment::env_num_file(-1); //Number of files analyzed. This is assigned in the program from the get go. -1 means all of them
		setup->Environment::env_skim(true);//Were the files initially skimmed
		setup->Environment::env_skim_type(1); //What type of skim was performed? See documentation 
		setup->Environment::env_efficiency(false); //Efficiencies of detector used
		setup->Environment::env_empty(false);//Is this an empty target run?
		setup->Environment::env_COM(false); //Are the four vectors already put into the center of mass?
		//For cut versions look setup->Environment:: the notes on the analysis
		setup->Environment::env_fid_ver(0, 1); //Which electron fiducial cut version was used? See documentation on fiducial cut history
		setup->Environment::env_fid_ver(1, 1);//Which proton fiducial cut version was used? See documentation on fiducial cut history
		setup->Environment::env_fid_ver(2, 1);//Which pip fiducial cut version was used? See documentation on fiducial cut history
		setup->Environment::env_fid_ver(3, 1);//Which pim fiducial cut version was used? See documentation on fiducial cut history
		setup->Environment::env_dt_ver(0, 0);//which proton delta t cut version was used? See documentation on delta t cuts
		setup->Environment::env_dt_ver(1, 0);//which pip delta t cut version was used? See documentation on delta t cuts
		setup->Environment::env_dt_ver(2, 0);//which pim delta t cut version was used? See documentation on delta t cuts
		setup->Environment::env_MM_ver(0, 1);//Which P missing cut version was used? See documentation on missing mass cuts
		setup->Environment::env_MM_ver(1, 1);//Which Pip missing cut version was used? See documentation on missing mass cuts
		setup->Environment::env_MM_ver(2, 1);//Which Pim missing cut version was used? See documentation on missing mass cuts
		setup->Environment::env_MM_ver(3, 1);//Which Z missing cut version was used? See documentation on missing mass cuts 
		
		setup->Environment::env_ec_ver(0); //Version of Min EC cut. 
		setup->Environment::env_sf_ver(0); //Version of Sampling Fraction cut. 
		//Eid Cuts applied
		setup->Environment::env_eid_fid(true); //Did eid include a fiducial cut?
		setup->Environment::env_eid_ec(true); //Did eid include a min ec cut?
		
		setup->Environment::env_eid_sf(true); //Did eid include a sampling fraction cut?
		//Hadron cuts applied {p,pip,pim}
		setup->Environment::env_hid_fid(0, true);//Did proton ID include a fiducial cut?
		setup->Environment::env_hid_fid(1, true);//Did pip ID include a fiducial cut?
		setup->Environment::env_hid_fid(2, true);//Did pim ID include a fiducial cut?
		setup->Environment::env_hid_dt(0, true);//Did proton ID include a delta t cut?
		setup->Environment::env_hid_dt(1, true);//Did pip ID include a delta t cut?
		setup->Environment::env_hid_dt(2, true);//Did pim ID include a delta t cut?
		setup->Environment::env_hid_beta(0, true);//Did proton ID include a beta cut?
		setup->Environment::env_hid_beta(1, true);//Did pip ID include a beta cut?
		setup->Environment::env_hid_beta(2, true);//Did pim ID include a beta cut?
		setup->Environment::env_hid_e(true); //Did hid cut out a band of electrons for pim ID?

		//Were fiducial cuts plotted? {e,p,pip,pim}
		setup->Environment::env_fid_plot(0, true); //Are electron fiducial plots being created?
		setup->Environment::env_fid_plot(1, true); //Are proton fiducial plots being created?
		setup->Environment::env_fid_plot(2, true);	// Are Pip fiducial plots being created?
		setup->Environment::env_fid_plot(3, true);//Are PIM fiducial plots being created?
		//Were delta t plots made performed {p,pip,pim}
		setup->Environment::env_dt_plot(0, true);//Are elec delta t plots being created?
		setup->Environment::env_dt_plot(1, true);//Are pro delta t plots being created?
		setup->Environment::env_dt_plot(2, true); //Are PIp delta t plots being created?
		setup->Environment::env_dt_plot(3, true); //Are PIM delta t plots being created?
		
	  	setup->Environment::env_ec_plot(true); //were Min EC plots made?
	  	setup->Environment::env_sf_plot(true); //Were Sampling Fraction plots made?
	  	setup->Environment::env_MM_plot(0,true); //Were missing mass plots made? {pmiss,pipmiss,pimmiss,zmiss}
	  	setup->Environment::env_MM_plot(1,true);
	  	setup->Environment::env_MM_plot(2,true);
	  	setup->Environment::env_MM_plot(3,true);
	  	setup->Environment::env_WQ2_plot(true); //Were W Qsquared plots made?
	  	setup->Environment::env_p_dep_plot(true); //Is there momentum dependence in the relevant plots?
	  	setup->Environment::env_W_dep_plot(true); //Is there W dependence in the relevant plots? 
	  	setup->Environment::env_Friend_plot(true);//Construct the multi dimensional histogram and fill it
  	}else if(run_type >= 5 && run_type!=7){
  		setup->Environment::env_npart(4); //Number of particles allowed in an event (-1 means no limits)
  		setup->Environment::env_Qmin(Q2minAna); //Lower range of Q squared allowed
		setup->Environment::env_Qmax(Q2maxAna); //Upper range of Q squared allowed
		setup->Environment::env_Wmin(WminAna); //Lower limit on W allowed
		setup->Environment::env_Wmax(WmaxAna); //Upper limit on W allowed 
  		setup->Environment::env_WQ2_plot(true); //Were W Qsquared plots made?
  		setup->Environment::env_MM_plot(0,false); //Were missing mass plots made? {pmiss,pipmiss,pimmiss,zmiss}
	  	setup->Environment::env_MM_plot(1,false);
	  	setup->Environment::env_MM_plot(2,false);
	  	setup->Environment::env_MM_plot(3,false);

  	}
}

void Setup::make_envi_file(const std::string& output_name, std::shared_ptr<Environment> envi){
	std::string curr_dir = fun::get_current_dir();
	std::string envi_name = "$curr/$name/$name_envi.txt";
	fun::replace(envi_name, "$curr", curr_dir);
	fun::replace(envi_name, "$name", output_name);
	fun::replace(envi_name, "$name", output_name);
	std::cout<<"	envi file name: " <<envi_name <<std::endl;
	std::ofstream envi_f;
	envi_f.open(envi_name);//main.hpp
	time_t now = time(0);
	char* dt = ctime(&now);
	if(!envi_f){
		std::cout<<"The file is not actually being made" <<std::endl;
	}else{
		envi_f <<"The environment for the " <<output_name <<" run" <<std::endl
		<<"Date Analysis was done at " <<dt <<std::endl  
		
		<<std::endl <<"--Data that was being looked at--" <<std::endl
		<<"Data set: " <<envi->Environment::was_data_set() <<std::endl //Which data set? {1,2}->{e16,e1f} //Not Setting
		<<"Sim: " <<envi->Environment::was_sim() <<std::endl  //Was this simulation data?
		<<"Empty Run: " <<envi->Environment::was_empty() <<std::endl//Is this an empty target run?
		<<"Golden Runs: " <<envi->Environment::was_golden() <<std::endl//Was this a golden run?
		<<"Fraction of Run: " <<envi->Environment::was_frac() <<std::endl //Fraction of the run analyzed
		<<"Number of files: " <<envi->Environment::was_num_file() <<std::endl //Number of files analyzed
		<<"Skim: " <<envi->Environment::was_skim() <<std::endl//Were the files initially skimmed
		<<"Type of Skim " <<envi->Environment::was_skim_type() <<std::endl //What type of skim was performed? See documentation 
		
		<<std::endl <<"--Kinematic Cuts--" <<std::endl
		<<"npart cut: " <<envi->Environment::was_npart() <<std::endl //Number of particles allowed in an event 
		<<"Q2min: " <<envi->Environment::was_Qmin() <<std::endl //Lower range of Q squared allowed
		<<"Q2max " <<envi->Environment::was_Qmax() <<std::endl //Upper range of Q squared allowed
		<<"Wmin " <<envi->Environment::was_Wmin() <<std::endl //Lower limit on W allowed
		<<"Wmax " <<envi->Environment::was_Wmax() <<std::endl //Upper limit on W allowed 
		
		<<std::endl <<"--Sanity Cut Parameters--" <<std::endl
		<<"DC hit: " <<envi->Environment::was_dc_hit() <<std::endl //Was a hit on the DC required?
		<<"CC hit:  " <<envi->Environment::was_cc_hit() <<std::endl  //Was a hit on the CC required for electrons?
		<<"EC hit: " <<envi->Environment::was_ec_hit() <<std::endl //Was a hit on the EC required for electrons?
		<<"SC hit: " <<envi->Environment::was_sc_hit() <<std::endl //Was a hit on the SC required?
		<<"Stat cut: " <<envi->Environment::was_stat_cut() <<std::endl //Was a nonzero stat required?
		<<"DC_Stat cut: " <<envi->Environment::was_dc_stat_cut() <<std::endl  //Was a nonzero dc_stat required?
		
		<<std::endl <<"--Topologies Included--" <<std::endl
		<<"P miss: " <<envi->Environment::was_top(0) <<std::endl
		<<"	P MM Version: " <<envi->Environment::was_MM_ver(0) <<std::endl
		<<"Pip miss: " <<envi->Environment::was_top(1) <<std::endl
		<<"	PIP MM Version: " <<envi->Environment::was_MM_ver(1) <<std::endl
		<<"Pim miss: " <<envi->Environment::was_top(2) <<std::endl
		<<"	PIM MM Version: " <<envi->Environment::was_MM_ver(2) <<std::endl
		<<"Zero miss: " <<envi->Environment::was_top(3) <<std::endl //Which topologies are included {pmiss,pipmiss,pimmiss,zeromiss,combined}
		<<"	Z MM Version: " <<envi->Environment::was_MM_ver(3) <<std::endl//Version of Missing Mass cut. 
		
		<<std::endl <<"--Corrections performed-- " <<std::endl
		<<"P corr: " <<envi->Environment::was_p_corr() <<std::endl //Were momentum corrections performed?
		<<"Efficiency Cuts: " <<envi->Environment::was_efficiency() <<std::endl //Efficiencies of detector used
		<<"COM Conversion: " <<envi->Environment::was_COM() <<std::endl //Are the four vectors already put envi->Environment:: the center of mass?

		<<std::endl <<"--Electron ID Cuts-- " <<std::endl
		<<"E fid: " <<envi->Environment::was_fid_plot(0) <<std::endl
		<<"EID Fid: " <<envi->Environment::was_eid_fid() <<std::endl //Did eid include a fiducial cut?
		<<"E Fid Version: " <<envi->Environment::was_fid_ver(0) <<std::endl
		<<"EID EC: " <<envi->Environment::was_eid_ec() <<std::endl //Did eid include a min ec cut?
		<<"E EC Version: " <<envi->Environment::was_ec_ver() <<std::endl //Version of Min EC cut. 
		<<"EID CC: " <<envi->Environment::was_eid_cc() <<std::endl //Did eid include a min cc cut?
		<<"E CC Version: " <<envi->Environment::was_cc_ver() <<std::endl //Version of Min-CC cut. 
		<<"EID SF: " <<envi->Environment::was_eid_sf() <<std::endl //Did eid include a sampling fraction cut?
		<<"E SF Version: " <<envi->Environment::was_sf_ver() <<std::endl //Version of Sampling Fraction cut.

		<<std::endl <<"--Proton ID Cuts-- " <<std::endl
		//<<"P Fid: " <<envi->Environment::was_fid(1) <<std::endl
		<<"P Fid Version: " <<envi->Environment::was_fid_ver(1) <<std::endl
		//<<"P dt: " <<envi->Environment::was_dt(0) <<std::endl
		<<"P DT Version: " <<envi->Environment::was_dt_ver(0) <<std::endl
		<<"Proton/Pip resolved: " <<envi->Environment::was_pvpip() <<std::endl //Were procedures allowed for dual identified positive hadrons?
		<<"P Fid: " <<envi->Environment::was_hid_fid(0) <<std::endl
		<<"P DT: " <<envi->Environment::was_hid_dt(0) <<std::endl
		<<"P Beta: " <<envi->Environment::was_hid_beta(0) <<std::endl

		<<std::endl <<"--PIP ID Cuts-- " <<std::endl 
		//<<"Pip Fid: " <<envi->Environment::was_fid(2) <<std::endl
		<<"PIP Fid Version: " <<envi->Environment::was_fid_ver(2) <<std::endl
		//<<"Pip dt: " <<envi->Environment::was_dt(1) <<std::endl
		<<"PIP DT Version: " <<envi->Environment::was_dt_ver(1) <<std::endl
		<<"PIP Fid: " <<envi->Environment::was_hid_fid(1) <<std::endl
		<<"PIP DT: " <<envi->Environment::was_hid_dt(1) <<std::endl
		<<"PIP Beta: " <<envi->Environment::was_hid_beta(1) <<std::endl

		<<std::endl <<"--PIM ID Cuts-- " <<std::endl
		//<<"Pim fid: " <<envi->Environment::was_hid_fid(2) <<std::endl//Were fiducial cuts performed {e,p,pip,pim}
		<<"Pim Fid Version: " <<envi->Environment::was_fid_ver(3) <<std::endl//Version of the fiducial cut. 
		//<<"Pim dt: " <<envi->Environment::was_hiddt(2) <<std::endl //Were delta t cuts performed {p,pip,pim}
		<<"PIM DT Version: " <<envi->Environment::was_dt_ver(2) <<std::endl//Version of delta t cut. 
		<<"PIM Fid: " <<envi->Environment::was_hid_fid(2) <<std::endl//Did hid include a fiducial cut?
		<<"PIM DT: " <<envi->Environment::was_hid_dt(2) <<std::endl//Did hid include a delta t cut?
		<<"PIM Beta: " <<envi->Environment::was_hid_beta(2) <<std::endl//Did hid include a beta cut?
		<<"PIM direct Electron cut " <<envi->Environment::was_hid_e() <<std::endl //Did hid cut out a band of electrons for pim ID?
		
		<<std::endl <<"--Plots that were made-- " <<std::endl
		<<"	Event plots: " <<std::endl
		<<"W vs Qsquared plots: " <<envi->Environment::was_WQ2_plot() <<std::endl
		<<"	Electron Plots: " <<std::endl
		<<"Fiducial plots: " <<envi->Environment::was_fid_plot(0) <<std::endl
		<<"CC plots: " <<envi->Environment::was_cc_plot() <<std::endl
		<<"EC plots: " <<envi->Environment::was_ec_plot() <<std::endl
		<<"Sampling Fraction plots: " <<envi->Environment::was_sf_plot() <<std::endl
		<<"	Proton Plots: " <<std::endl
		<<"Fiducial plots: " <<envi->Environment::was_fid_plot(1) <<std::endl
		<<"Delta T plots: " <<envi->Environment::was_dt_plot(0) <<std::endl
		<<"	Pip Plots: " <<std::endl
		<<"Fiducial plots: " <<envi->Environment::was_fid_plot(2) <<std::endl
		<<"Delta T plots: " <<envi->Environment::was_dt_plot(1) <<std::endl
		<<"	Pim Plots: " <<std::endl
		<<"Fiducial plots: " <<envi->Environment::was_fid_plot(3) <<std::endl
		<<"Delta T plots: " <<envi->Environment::was_dt_plot(2) <<std::endl
		<<"	Topologies plotted: " <<std::endl
		<<"PMiss Plot: " <<envi->Environment::was_MM_plot(0) <<std::endl
		<<"PipMiss Plot: " <<envi->Environment::was_MM_plot(1) <<std::endl
		<<"PimMiss Plot: " <<envi->Environment::was_MM_plot(2) <<std::endl
		<<"ZMiss Plot: " <<envi->Environment::was_MM_plot(3) <<std::endl;
		//<<"		Corrections Plotted: " <<std::endl
		//<<"P Corr: " <<envi->Environment::was_pcorr_plot() <<std::endl
		//<<"Efficiencies plotted: " <<envi->Environment::was_eff_plot() <<std::endl

		
		envi_f.close();
	}
}