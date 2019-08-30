#ifndef NTUPLE_H_GUARD
#define	NTUPLE_H_GUARD


#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "event_class.hpp"
#include "functions.hpp"
#include "constants.hpp"

//Making these TTrees
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

public:
	forest(std::string tree_file_name);
	~forest();

	void mktree(int thread_id);
	void Fill_Tree(std::shared_ptr<Event> event_friend, int event_n, int thread_id);

};



#endif