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
		_px[i] = event_friend->Event::Get_px(i);
		_py[i] = event_friend->Event::Get_py(i);
		_pz[i] = event_friend->Event::Get_pz(i);
		_p0[i] = event_friend->Event::Get_p0(i);
		_pid[i] = event_friend->Event::Get_pid(i);
	}
	_hel = event_friend->Event::Get_hel();
	_top = event_friend->Event::Get_top();
	_the_tree->TTree::Fill();	
}

void forest::fill_evnt(int event_n){
	_evnt = event_n;
	_the_tree->TTree::Fill();
}

void forest::fill_apart(int apart_n){
	_apart = apart_n;
	_the_tree->TTree::Fill();
}
void forest::fill_px(float px_n, int i){
	_px[i] = px_n; 
	_the_tree->TTree::Fill();
}
void forest::fill_py(float py_n, int i){
	_py[i] = py_n;
	_the_tree->TTree::Fill();
}

void forest::fill_pz(float pz_n, int i){
	_pz[i] = pz_n;
	_the_tree->TTree::Fill();
}
void forest::fill_p0(float p0_n, int i){
	_py[i] = p0_n;
	_the_tree->TTree::Fill();
}
void forest::fill_pid(int pid_n, int i){
	_pid[i] = pid_n;
	_the_tree->TTree::Fill();
}
void forest::fill_hel(int hel_n){
	_hel = hel_n;
	_the_tree->TTree::Fill();
}
void forest::fill_top(int top_n){
	_top = top_n; 
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
