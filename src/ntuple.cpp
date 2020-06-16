#include "ntuple.hpp"


forest::forest(int is_alive){
	//The main Tree that will be made into a root file
	char tree_name[100];

	_the_tree = new TTree("k10","Tree to hold Event Fourvectors");//
	//_the_tree->Branch("EventBranch",&the_event,"evnt/I:apart/I:px/F:py/F:pz/F:p0/F:pid/I:hel/I:top/I");
	_the_tree->Branch("evnt",&_evnt,"evnt/I");
	_the_tree->Branch("apart",&_apart,"apart/I");
	_the_tree->Branch("bpart",&_bpart,"bpart/I");
	_the_tree->Branch("px",&_px,"px[apart]/F");
	_the_tree->Branch("py",&_py,"py[apart]/F");
	_the_tree->Branch("pz",&_pz,"pz[apart]/F");
	_the_tree->Branch("p0",&_p0,"p0[apart]/F");
	_the_tree->Branch("pid",&_pid,"pid[apart]/I");
	_the_tree->Branch("hel",&_hel,"hel/I");
	_the_tree->Branch("top",&_top,"top/I");
	_the_tree->Branch("fc_tot",&_fc_tot,"fc_tot/I");
	_the_tree->Branch("MM",&_MM,"MM[bpart]/F");
	_the_tree->Branch("theta",&_theta,"theta[bpart]/F");
	_the_tree->Branch("alpha",&_alpha,"alpha[bpart]/F");
	_the_tree->Branch("run_type",&_run_type,"run_type/I");
	_the_tree->Branch("weight",&_weight,"weight/I");
	Float_t _MM_sp[3]= {NAN,NAN,NAN};//{p/pip,p/pim,pip/pim}
	Float_t _theta_sp[3]{NAN,NAN,NAN};//{p,pip,pim}
	Float_t _alpha_sp[3]{NAN,NAN,NAN};//[{pim,p},{pp,pip}],[{p,pp},{pip,pim}],[{pip,p},{pp,pim}]
	Int_t _run_type= 0;//{1,2,3,4}->{e16,e1f,e16sim,e1f sim} 
	//The individual thread trees that will be merged
	for(int thread_id = 0; thread_id < NUM_THREADS; thread_id++){
		sprintf(tree_name,"t%d",thread_id);
		_a_tree[thread_id] = new TTree(tree_name,"Tree to hold Thread Event Fourvectors");//
		//_a_tree[thread_id]->Branch("EventBranch",&the_eventb,"evntb/I:apartb/I:pxb/F:pyb/F:pzb/F:p0b/F:pidb/I:helb/I:topb/I");
		_a_tree[thread_id]->Branch("evnt",&(_evntb[thread_id]),"evnt/I");
		_a_tree[thread_id]->Branch("apart",&(_apartb[thread_id]),"apart/I");
		_a_tree[thread_id]->Branch("bpart",&(_bpartb[thread_id]),"bpart/I");
		_a_tree[thread_id]->Branch("px",&(_pxb[thread_id]),"px[apart]/F");
		_a_tree[thread_id]->Branch("py",&(_pyb[thread_id]),"py[apart]/F");
		_a_tree[thread_id]->Branch("pz",&(_pzb[thread_id]),"pz[apart]/F");
		_a_tree[thread_id]->Branch("p0",&(_p0b[thread_id]),"p0[apart]/F");
		_a_tree[thread_id]->Branch("pid",&(_pidb[thread_id]),"pid[apart]/I");
		_a_tree[thread_id]->Branch("hel",&(_helb[thread_id]),"hel/I");
		_a_tree[thread_id]->Branch("top",&(_topb[thread_id]),"top/I");
		_a_tree[thread_id]->Branch("fc_tot",&(_fc_totb[thread_id]),"fc_tot/I");
		_a_tree[thread_id]->Branch("MM",&_MMb[thread_id],"MM[bpart]/F");
		_a_tree[thread_id]->Branch("theta",&_thetab[thread_id],"theta[bpart]/F");
		_a_tree[thread_id]->Branch("alpha",&_alphab[thread_id],"alpha[bpart]/F");
		_a_tree[thread_id]->Branch("run_type",&_run_typeb[thread_id],"run_type/I");
		_a_tree[thread_id]->Branch("weight",&_weightb[thread_id],"weight/I");
		//The B Trees
		/*
		_b_tree[thread_id] = std::make_shared<TTree>("TREEsb","Tree to hold Thread Event Fourvectors");//
		_b_tree[thread_id]->Branch("evnt",&_evntc[thread_id],"evnt/I");
		_b_tree[thread_id]->Branch("apart",&_apartc[thread_id],"apart/I");
		_b_tree[thread_id]->Branch("px",&_pxc[thread_id],"px[apart]/F");
		_b_tree[thread_id]->Branch("py",&_pyc[thread_id],"py[apart]/F");
		_b_tree[thread_id]->Branch("pz",&_pzc[thread_id],"pz[apart]/F");
		_b_tree[thread_id]->Branch("p0",&_p0c[thread_id],"p0[apart]/F");
		_b_tree[thread_id]->Branch("pid",&_pidc[thread_id],"pid[apart]/F");
		_b_tree[thread_id]->Branch("hel",&_helc[thread_id],"hel/I");
		_b_tree[thread_id]->Branch("top",&_topc[thread_id],"top/I");*/
	}
	alive = is_alive; 
}

//forest::~forest(){ this->Write();}//Including this and the forest::Write() are necessary for properly writing the TTree to the TFile 


void forest::Write(){
	std::cout<<"Writing Event Tree: ";
	_tree_file->cd();
	TDirectory* _almanac = _tree_file->mkdir("k10");
	_almanac->cd();
	_the_tree->Write();
	_tree_file->Close();
	std::cout<<" Written and Closed" <<std::endl;
}

/*
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
}*/

void forest::mkfile(std::string tree_file_name){
	_tree_file = fun::Name_Tree_File(tree_file_name);//The output rootfile containng the tree
}


void forest::Fill_Thread_Tree(std::shared_ptr<Event_Class> event_friend, int event_n, int thread_id){
	//the_eventb[thread_id] = new Event();
	//std::cout<<std::endl <<"Tree being filled in thread " <<thread_id; 
	//if(event_friend->Event_Class::is_valid()){
		//std::cout<<std::endl<<"Filling Tree trying for event: " <<event_n;
		TLorentzVector k[6]; 
		_evntb[thread_id] = event_n; 
		_apartb[thread_id] = 6;
		_bpartb[thread_id] = 3;  
		for(int i = 0; i<6; i++){
			_pxb[thread_id][i] = event_friend->Event_Class::Get_px(i);
			//std::cout<<"Here is pxb " <<_pxb[thread_id][i] <<" for thread " <<thread_id <<" and index " <<i <<std::endl;
			_pyb[thread_id][i] = event_friend->Event_Class::Get_py(i);
			_pzb[thread_id][i] = event_friend->Event_Class::Get_pz(i);
			_p0b[thread_id][i] = event_friend->Event_Class::Get_p0(i);
			_pidb[thread_id][i] = event_friend->Event_Class::Get_pid(i);
			k[i] = {_pxb[thread_id][i],_pyb[thread_id][i],_pzb[thread_id][i],_p0b[thread_id][i]};
			//std::cout<<std::endl <<"pid? " <<_pidb[thread_id][i];
		}
		_helb[thread_id] = event_friend->Event_Class::Get_hel();
		_topb[thread_id] = event_friend->Event_Class::Get_top();
		_run_typeb[thread_id] = event_friend->Event_Class::Get_run_type();
		_MMb[thread_id][0] = (k[3]+k[4]).Mag();
		_MMb[thread_id][1] = (k[3]+k[5]).Mag();
		_MMb[thread_id][2] = (k[5]+k[4]).Mag();
		_a_tree[thread_id]->TTree::Fill();
		//std::cout<<"   successfully filled?"; 
	//}
}

void forest::scan_thread_tree(int thread_id){
	_a_tree[thread_id]->TTree::Scan();
}

void forest::fill_evnt(int event_n, int thread_id){
	_evnt = event_n;
	_a_tree[thread_id]->TTree::Fill();
	//_a_tree[thread_id]->TTree::Scan();
}

void forest::fill_apart(int apart_n, int thread_id){
	_apart = apart_n;
	_a_tree[thread_id]->TTree::Fill();
}
void forest::fill_px(float px_n, int i, int thread_id){
	_px[i] = px_n; 
	_a_tree[thread_id]->TTree::Fill();
}
void forest::fill_py(float py_n, int i, int thread_id){
	_py[i] = py_n;
	_a_tree[thread_id]->TTree::Fill();
}

void forest::fill_pz(float pz_n, int i, int thread_id){
	_pz[i] = pz_n;
	_a_tree[thread_id]->TTree::Fill();
}
void forest::fill_p0(float p0_n, int i, int thread_id){
	_py[i] = p0_n;
	_a_tree[thread_id]->TTree::Fill();
}
void forest::fill_pid(int pid_n, int i, int thread_id){
	_pid[i] = pid_n;
	_a_tree[thread_id]->TTree::Fill();
}
void forest::fill_hel(int hel_n, int thread_id){
	_hel = hel_n;
	_a_tree[thread_id]->TTree::Fill();
}
void forest::fill_top(int top_n, int thread_id){
	_top = top_n; 
	_a_tree[thread_id]->TTree::Fill();
}
/*
void forest::move_forest_indoors(){
	for(int thread_id = 0; thread_id < NUM_THREADS; i++){
		_evntc[thread_id] = _evntb[thread_id]; 
		for(int i; i<4; i++){
			_apartc[thread_id] = apartb[thread_id]; 
			_pxc[thread_id][i] = _pxb[thread_id][i];
			_pyc[thread_id][i] = _pyb[thread_id][i];
			_pzc[thread_id][i] = _pzb[thread_id][i];
			_p0c[thread_id][i] = _p0b[thread_id][i];
			_pidc[thread_id][i] = _pidb[thread_id][i];
		}
		_helc[thread_id] = _helb[thread_id];
		_topc[thread_id] = _topb[thread_id];
		_b_tree[thread_id]->TTree::Fill();	
	}
}*/

void forest::Grow_Forest(){
	TList the_forest; 
	int num_ev[NUM_THREADS];
	for(int i = 0; i<NUM_THREADS ; i++){
		the_forest.TList::Add(_a_tree[i]);
	}
	_the_tree = TTree::MergeTrees(&the_forest,"");
}
	/*
void forest::mkforest(){
  _the_tree = TTree::MergeTrees(_the_forest);
  _the_tree->TTree::Write();
}*/
