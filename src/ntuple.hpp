#ifndef NTUPLE_H_GUARD
#define	NTUPLE_H_GUARD


#include "TFile.h"
#include "TTree.h"
#include "TList.h"
#include "TDirectory.h"
#include "event_class.hpp"
#include "functions.hpp"
#include "constants.hpp"

//Output Event files which will contain event selected information 
class forest{
private:
	Int_t _evnt; //The #event for the given file
	Int_t _apart;//The particle in each event  
	Float_t _px[4]; 
	Float_t _py[4]; 
	Float_t _pz[4];
	Float_t _p0[4];//Energy of the particle 
	Int_t _pid[4];//particle ID 
	Int_t _hel; //helicity
	Int_t _top; //Topology {pmiss,pipmiss,pimmiss,zero} -> {1,2,3,4}

	std::shared_ptr<TFile> _tree_file; 
	std::shared_ptr<TTree> _the_tree;
	TTree _a_tree; 
	TList* _the_forest;
	//std::shared_ptr<TDirectory> _almanac; 

public:
	forest();
	~forest();

	void mktree(int thread_id);
	void mkfile(std::string tree_file_name);
	void Fill_Tree(std::shared_ptr<Event> event_friend, int event_n, int thread_id);//For multithreading there need to be separate trees for each thread
	void fill_evnt(int event_);
	void fill_apart(int apart_);
	void fill_px(float px_n, int i);
	void fill_py(float py_n, int i);
	void fill_pz(float pz_n, int i);
	void fill_p0(float p0_n, int i);
	void fill_pid(int pid_n, int i);
	void fill_hel(int hel_n);
	void fill_top(int top_n);

	//void Grow_Forest(std::shared_ptr<TTree> a_tree);
	//void mkforest();
};



#endif