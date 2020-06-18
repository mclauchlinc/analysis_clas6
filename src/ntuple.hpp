#ifndef NTUPLE_H_GUARD
#define	NTUPLE_H_GUARD


#include "TFile.h"
#include "TTree.h"
#include "TList.h"
#include "TDirectory.h"
#include "event_class.hpp"
#include "functions.hpp"
#include "constants.hpp"
#include <iostream>
#include "event.hpp"


//Output Event files which will contain event selected information 
class forest{
private:
	Int_t _evnt = 0; //The #event for the given file
	Int_t _apart = 6;//The particle in each event  
	Int_t _bpart = 3; 
	Float_t _px[6] = {NAN,NAN,NAN,NAN,NAN,NAN};//{beam,target,e',p',pip,pim}
	Float_t _py[6]= {NAN,NAN,NAN,NAN,NAN,NAN}; //{beam,target,e',p',pip,pim}
	Float_t _pz[6]= {NAN,NAN,NAN,NAN,NAN,NAN};//{beam,target,e',p',pip,pim}
	Float_t _p0[6]= {NAN,NAN,NAN,NAN,NAN,NAN};////{beam,target,e',p',pip,pim}Energy of the particle 
	Int_t _pid[6]= {0,0,0,0,0,0};//particle ID 
	Int_t _hel= 0; //helicity
	Int_t _top= 0; //Topology {pmiss,pipmiss,pimmiss,zero} -> {1,2,3,4}
	Int_t _fc_tot = 0;
	Float_t _MM[3]= {NAN,NAN,NAN};//{p/pip,p/pim,pip/pim}
	Float_t _theta[3]{NAN,NAN,NAN};//{p,pip,pim}
	Float_t _alpha[3]{NAN,NAN,NAN};//[{pim,p},{pp,pip}],[{p,pp},{pip,pim}],[{pip,p},{pp,pim}]
	Int_t _run_type= 0;//{1,2,3,4,5,6,7,8}->{e16,e1f,e16sim,e1f sim, e16sim_thrown, e1f sim_thrown,e16empty,e1fempty}
	Float_t _weight = NAN; 
	Int_t _COM = -1; 

	//Variables for the individual threads
	Int_t _evntb[NUM_THREADS]; //The #event for the given file
	Int_t _apartb[NUM_THREADS];//The particle in each event
	Int_t _bpartb[NUM_THREADS];  
	Float_t _pxb[NUM_THREADS][6];
	Float_t _pyb[NUM_THREADS][6]; 
	Float_t _pzb[NUM_THREADS][6];
	Float_t _p0b[NUM_THREADS][6];//Energy of the particle 
	Int_t _pidb[NUM_THREADS][6];//particle ID 
	Int_t _helb[NUM_THREADS]; //helicity
	Int_t _topb[NUM_THREADS]; //Topology {pmiss,pipmiss,pimmiss,zero} -> {1,2,3,4}
	Int_t _fc_totb[NUM_THREADS];
	Float_t _MMb[NUM_THREADS][3];//{p/pip,p/pim,pip/pim}
	Float_t _thetab[NUM_THREADS][3];//{p,pip,pim}
	Float_t _alphab[NUM_THREADS][3];//[{pim,p},{pp,pip}],[{p,pp},{pip,pim}],[{pip,p},{pp,pim}]
	Int_t _run_typeb[NUM_THREADS];//{1,2,3,4}->{e16,e1f,e16sim,e1f sim} 
	Float_t _weightb[NUM_THREADS]; 
	Int_t _COMb[NUM_THREADS];
	/*
	Int_t _evntc[NUM_THREADS]; //The #event for the given file
	Int_t _apartc[NUM_THREADS];//The particle in each event  
	Float_t _pxc[NUM_THREADS][4]; 
	Float_t _pyc[NUM_THREADS][4]; 
	Float_t _pzc[NUM_THREADS][4];
	Float_t _p0c[NUM_THREADS][4];//Energy of the particle 
	Int_t _pidc[NUM_THREADS][4];//particle ID 
	Int_t _helc[NUM_THREADS]; //helicity
	Int_t _topc[NUM_THREADS]; //Topology {pmiss,pipmiss,pimmiss,zero} -> {1,2,3,4}
	*/

	std::shared_ptr<TFile> _tree_file; //The root file that will be made
	//std::shared_ptr<TTree> 
	//Event* the_event; 
	//Event* the_eventb[NUM_THREADS];
	TTree* _the_tree; //Final tree for the root file
	TTree* _a_tree[NUM_THREADS];//Trees for the various threads 
	TList* _the_forest;
	//TTree* _b_tree[NUM_THREADS]; 
	//std::shared_ptr<TDirectory> _almanac; 
	int alive = 0; 

public:
	forest(int is_alive);//-1 thread_id gives the whole tree, but the others give pieces of trees
	//~forest();
	void Write();

	//void mktree(int thread_id);
	void mkfile(std::string tree_file_name);
	void Fill_Thread_Tree(std::shared_ptr<Event_Class> event_friend, int event_n, int thread_id);//For multithreading there need to be separate trees for each thread
	void Fill_Thread_Tree(Event event_, int event_n, int thread_id);//For multithreading there need to be separate trees for each thread
	void scan_thread_tree(int thread_id);
	void fill_evnt(int event_, int thread_id);
	void fill_apart(int apart_, int thread_id);
	void fill_px(float px_n, int i, int thread_id);
	void fill_py(float py_n, int i, int thread_id);
	void fill_pz(float pz_n, int i, int thread_id);
	void fill_p0(float p0_n, int i, int thread_id);
	void fill_pid(int pid_n, int i, int thread_id);
	void fill_hel(int hel_n, int thread_id);
	void fill_top(int top_n, int thread_id);

	void Grow_Forest();
	//void mkforest();
};



#endif