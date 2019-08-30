#include "ntuple.hpp"


forest::forest(){
};

forest::~forest(){
	_tree_file->Write();
}

void forest::mktree(int thread_id){
	char* tree_name;
	char* std_tree_name = "t"; 
	tree_name = std_tree_name + thread_id; 
	_the_tree = std::make_shared<TTree>(tree_name,"Tree to hold Event Fourvectors");//
	_the_tree->Branch("evnt",&_evnt,"evnt/I");
	_the_tree->Branch("apart",&_apart,"apart/I");
	_the_tree->Branch("px",&_evnt,"px[apart]/F");
	_the_tree->Branch("py",&_evnt,"py[apart]/F");
	_the_tree->Branch("pz",&_evnt,"pz[apart]/F");
	_the_tree->Branch("p0",&_evnt,"p0[apart]/F");
	_the_tree->Branch("pid",&_evnt,"pid[apart]/F");
	_the_tree->Branch("hel",&_hel,"hel/I");
	_the_tree->Branch("top",&_top,"top/I");	
}

void forest::mkfile(std::string tree_file_name){
	_tree_file = fun::Name_File(tree_file_name);//The output rootfile containng the tree
	TDirectory* _almanac = _tree_file->mkdir("physics_phorest");
	_almanac->cd();
}

void forest::Fill_Tree(std::shared_ptr<Event> event_friend, int event_n, int thread_id){
	_evnt = event_n; 
	for(int i; i<4; i++){
		_apart = i; 
		_px[i] = event_friend->Get_px(i);
		_py[i] = event_friend->Get_py(i);
		_pz[i] = event_friend->Get_pz(i);
		_p0[i] = event_friend->Get_p0(i);
		_pid[i] = event_friend->Get_pid(i);
	}
	_hel = event_friend->Get_hel();
	_top = event_friend->Get_top();
	_the_tree->TTree::Fill();	
}
/*
void forest::Grow_Forest(std::shared_ptr<TTree> a_tree){
	TTree* temp_tree = a_tree; 
	_the_forest->TList::Add(temp_tree);
}
	
void forest::mkforest(){
  _the_tree = TTree::MergeTrees(_the_forest);
  _the_tree->TTree::Write();
}*/
