#include "ntuple.hpp"


forest::forest(int is_alive, bool sim_){
	//The main Tree that will be made into a root file
	char tree_name[100];

	_the_tree = new TTree("r10","Tree to hold Event Fourvectors");//
	//_the_tree->Branch("EventBranch",&the_event,"evnt/I:apart/I:px/F:py/F:pz/F:p0/F:pid/I:hel/I:top/I");
	_the_tree->Branch("evnt",&_evnt[0],"evnt/I");
	_the_tree->Branch("apart",&_apart[0],"apart/I");
	_the_tree->Branch("bpart",&_bpart[0],"bpart/I");
	_the_tree->Branch("px",&_px[0],"px[apart]/F");
	_the_tree->Branch("py",&_py[0],"py[apart]/F");
	_the_tree->Branch("pz",&_pz[0],"pz[apart]/F");
	_the_tree->Branch("p0",&_p0[0],"p0[apart]/F");
	_the_tree->Branch("pid",&_pid[0],"pid[apart]/I");
	_the_tree->Branch("hel",&_hel[0],"hel/I");
	_the_tree->Branch("top",&_top[0],"top/I");
	_the_tree->Branch("fc_tot",&_fc_tot[0],"fc_tot/I");
	_the_tree->Branch("MM",&_MM[0],"MM[bpart]/F");
	_the_tree->Branch("theta",&_theta[0],"theta[bpart]/F");
	_the_tree->Branch("alpha",&_alpha[0],"alpha[bpart]/F");
	_the_tree->Branch("run_type",&_run_type[0],"run_type/I");
	_the_tree->Branch("weight",&_weight[0],"weight/F");
	_the_tree->Branch("COM",&_COM[0],"COM/I");
	//Float_t _MM_sp[3]= {NAN,NAN,NAN};//{p/pip,p/pim,pip/pim}
	//Float_t _theta_sp[3]{NAN,NAN,NAN};//{p,pip,pim}
	//Float_t _alpha_sp[3]{NAN,NAN,NAN};//[{pim,p},{pp,pip}],[{p,pp},{pip,pim}],[{pip,p},{pp,pim}]
	//Int_t _run_type= 0;//{1,2,3,4}->{e16,e1f,e16sim,e1f sim} 
	//The individual thread trees that will be merged
	for(int thread_id = 0; thread_id < NUM_THREADS; thread_id++){
		sprintf(tree_name,"k%d",thread_id);
		_a_tree[thread_id] = new TTree(tree_name,"Tree to hold Thread Event Fourvectors");//
		//_a_tree[thread_id]->Branch("EventBranch",&the_eventb,"evntb/I:apartb/I:pxb/F:pyb/F:pzb/F:p0b/F:pidb/I:helb/I:topb/I");
		_a_tree[thread_id]->Branch("evnt",&(_evntb[thread_id][0]),"evnt/I");
		_a_tree[thread_id]->Branch("apart",&(_apartb[thread_id][0]),"apart/I");
		_a_tree[thread_id]->Branch("bpart",&(_bpartb[thread_id][0]),"bpart/I");
		_a_tree[thread_id]->Branch("px",&(_pxb[thread_id][0]),"px[apart]/F");
		_a_tree[thread_id]->Branch("py",&(_pyb[thread_id][0]),"py[apart]/F");
		_a_tree[thread_id]->Branch("pz",&(_pzb[thread_id][0]),"pz[apart]/F");
		_a_tree[thread_id]->Branch("p0",&(_p0b[thread_id][0]),"p0[apart]/F");
		_a_tree[thread_id]->Branch("pid",&(_pidb[thread_id][0]),"pid[apart]/I");
		_a_tree[thread_id]->Branch("hel",&(_helb[thread_id][0]),"hel/I");
		_a_tree[thread_id]->Branch("top",&(_topb[thread_id][0]),"top/I");
		_a_tree[thread_id]->Branch("fc_tot",&(_fc_totb[thread_id][0]),"fc_tot/I");
		_a_tree[thread_id]->Branch("MM",&_MMb[thread_id][0],"MM[bpart]/F");
		_a_tree[thread_id]->Branch("theta",&_thetab[thread_id][0],"theta[bpart]/F");
		_a_tree[thread_id]->Branch("alpha",&_alphab[thread_id][0],"alpha[bpart]/F");
		_a_tree[thread_id]->Branch("run_type",&_run_typeb[thread_id][0],"run_type/I");
		_a_tree[thread_id]->Branch("weight",&_weightb[thread_id][0],"weight/F");
		_a_tree[thread_id]->Branch("COM",&_COMb[thread_id][0],"COM/I");
	}
	if(sim_){
		_thrown_tree = new TTree("t10","Tree to hold Event Fourvectors");//
		//_thrown_tree->Branch("EventBranch",&the_event,"evnt/I:apart/I:px/F:py/F:pz/F:p0/F:pid/I:hel/I:top/I");
		_thrown_tree->Branch("evnt",&_evnt[1],"evnt/I");
		_thrown_tree->Branch("apart",&_apart[1],"apart/I");
		_thrown_tree->Branch("bpart",&_bpart[1],"bpart/I");
		_thrown_tree->Branch("px",&_px[1],"px[apart]/F");
		_thrown_tree->Branch("py",&_py[1],"py[apart]/F");
		_thrown_tree->Branch("pz",&_pz[1],"pz[apart]/F");
		_thrown_tree->Branch("p0",&_p0[1],"p0[apart]/F");
		_thrown_tree->Branch("pid",&_pid[1],"pid[apart]/I");
		_thrown_tree->Branch("hel",&_hel[1],"hel/I");
		_thrown_tree->Branch("top",&_top[1],"top/I");
		_thrown_tree->Branch("fc_tot",&_fc_tot[1],"fc_tot/I");
		_thrown_tree->Branch("MM",&_MM[1],"MM[bpart]/F");
		_thrown_tree->Branch("theta",&_theta[1],"theta[bpart]/F");
		_thrown_tree->Branch("alpha",&_alpha[1],"alpha[bpart]/F");
		_thrown_tree->Branch("run_type",&_run_type[1],"run_type/I");
		_thrown_tree->Branch("weight",&_weight[1],"weight/F");
		_thrown_tree->Branch("COM",&_COM[1],"COM/I");
		//Float_t _MM_sp[3]= {NAN,NAN,NAN};//{p/pip,p/pim,pip/pim}
		//Float_t _theta_sp[3]{NAN,NAN,NAN};//{p,pip,pim}
		//Float_t _alpha_sp[3]{NAN,NAN,NAN};//[{pim,p},{pp,pip}],[{p,pp},{pip,pim}],[{pip,p},{pp,pim}]
		//Int_t _run_type= 0;//{1,2,3,4}->{e16,e1f,e16sim,e1f sim} 
		//The individual thread trees that will be merged
		for(int thread_id = 0; thread_id < NUM_THREADS; thread_id++){
			sprintf(tree_name,"u%d",thread_id);
			_athrown_tree[thread_id] = new TTree(tree_name,"Tree to hold Thread Event Fourvectors");//
			//_athrown_tree[thread_id]->Branch("EventBranch",&the_eventb,"evntb/I:apartb/I:pxb/F:pyb/F:pzb/F:p0b/F:pidb/I:helb/I:topb/I");
			_athrown_tree[thread_id]->Branch("evnt",&(_evntb[thread_id][1]),"evnt/I");
			_athrown_tree[thread_id]->Branch("apart",&(_apartb[thread_id][1]),"apart/I");
			_athrown_tree[thread_id]->Branch("bpart",&(_bpartb[thread_id][1]),"bpart/I");
			_athrown_tree[thread_id]->Branch("px",&(_pxb[thread_id][1]),"px[apart]/F");
			_athrown_tree[thread_id]->Branch("py",&(_pyb[thread_id][1]),"py[apart]/F");
			_athrown_tree[thread_id]->Branch("pz",&(_pzb[thread_id][1]),"pz[apart]/F");
			_athrown_tree[thread_id]->Branch("p0",&(_p0b[thread_id][1]),"p0[apart]/F");
			_athrown_tree[thread_id]->Branch("pid",&(_pidb[thread_id][1]),"pid[apart]/I");
			_athrown_tree[thread_id]->Branch("hel",&(_helb[thread_id][1]),"hel/I");
			_athrown_tree[thread_id]->Branch("top",&(_topb[thread_id][1]),"top/I");
			_athrown_tree[thread_id]->Branch("fc_tot",&(_fc_totb[thread_id][1]),"fc_tot/I");
			_athrown_tree[thread_id]->Branch("MM",&_MMb[thread_id][1],"MM[bpart]/F");
			_athrown_tree[thread_id]->Branch("theta",&_thetab[thread_id][1],"theta[bpart]/F");
			_athrown_tree[thread_id]->Branch("alpha",&_alphab[thread_id][1],"alpha[bpart]/F");
			_athrown_tree[thread_id]->Branch("run_type",&_run_typeb[thread_id][1],"run_type/I");
			_athrown_tree[thread_id]->Branch("weight",&_weightb[thread_id][1],"weight/F");
			_athrown_tree[thread_id]->Branch("COM",&_COMb[thread_id][1],"COM/I");
		}
	}
	
	alive = is_alive; 
}

//forest::~forest(){ this->Write();}//Including this and the forest::Write() are necessary for properly writing the TTree to the TFile 


void forest::Write(bool sim_){
	std::cout<<"Writing Event Tree\n";
	_tree_file->cd();
	TDirectory* _almanac1 = _tree_file->mkdir("r10");
	_almanac1->cd();
	_the_tree->Write();
	std::cout<<"Reconstructed Tree Written\n";
	if(sim_){
		TDirectory* _almanac2 = _tree_file->mkdir("t10");
		_almanac2->cd();
		_thrown_tree->Write();
		std::cout<<"Thrown Tree Written\n";

	}
	_tree_file->Close();//It seems that if there are multiple forests that closing this one closes other ones too? But if I don't close then it all gets very weird 
	std::cout<<"Event Tree Closed\n";
	
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

void forest::mkfile(std::string tree_file_name, bool thrown_){
	if(thrown_){
		_thr_tree_file = fun::Name_Tree_File(tree_file_name, thrown_);
	}else{
		_tree_file = fun::Name_Tree_File(tree_file_name, thrown_);
	}
	//The output rootfile containng the tree
}


/*void forest::Fill_Thread_Tree(std::shared_ptr<Event_Class> event_friend, int event_n, int thread_id){
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
}*/

void forest::Fill_Thread_Tree(Event event_, int event_n, int thread_id, bool thrown_){
		//std::cout<<"Filling tree: " <<event_n <<std::endl;
	int thr = 0; 
	if(thrown_){
		thr = 1;
	}	
		//TLorentzVector k[6]; 
		_evntb[thread_id][thr] = event_n; 
		_apartb[thread_id][thr] = 6;
		_bpartb[thread_id][thr] = 3; 
	
	

		bool want_COM = true;//Do you want things in Center of Mass Frame?

		for(int i = 0; i<6; i++){
			if(i > 1){
				_pxb[thread_id][thr][i] = event_.Event::Get_Px(i-2,want_COM);
				_pyb[thread_id][thr][i] = event_.Event::Get_Py(i-2,want_COM);
				_pzb[thread_id][thr][i] = event_.Event::Get_Pz(i-2,want_COM);
				_p0b[thread_id][thr][i] = event_.Event::Get_P0(i-2,want_COM);
				_pidb[thread_id][thr][i] = event_.Event::Get_PID(i-2);
				//k[i] = {_pxb[thread_id][thr][i],_pyb[thread_id][thr][i],_pzb[thread_id][thr][i],_p0b[thread_id][thr][i]};
			}else{
				if(i == 0){//beam
					_pxb[thread_id][thr][i] = event_.Event::Get_Beam_Comp(0,true);
					_pyb[thread_id][thr][i] = event_.Event::Get_Beam_Comp(1,true);
					_pzb[thread_id][thr][i] = event_.Event::Get_Beam_Comp(2,true);
					_p0b[thread_id][thr][i] = event_.Event::Get_Beam_Comp(3,true);
					_pidb[thread_id][thr][i] = ELECTRON;
				}else if(i ==1){ //Target
					_pxb[thread_id][thr][i] = event_.Event::Get_Target_Comp(0,true);
					_pyb[thread_id][thr][i] = event_.Event::Get_Target_Comp(1,true);
					_pzb[thread_id][thr][i] = event_.Event::Get_Target_Comp(2,true);
					_p0b[thread_id][thr][i] = event_.Event::Get_Target_Comp(3,true);
					_pidb[thread_id][thr][i] = PROTON;
				}
			}
		}
		_helb[thread_id][thr] = event_.Event::Get_Hel();
		_topb[thread_id][thr] = event_.Event::Get_Top();
		_run_typeb[thread_id][thr] = event_.Event::Get_Set();//{1,0}.{e16,e1f}
		for(int j = 0; j< 3; j++){
			_MMb[thread_id][thr][j] = event_.Event::Get_MMb(j);
			_thetab[thread_id][thr][j] = event_.Event::Get_Thetab(j);
			_alphab[thread_id][thr][j] = event_.Event::Get_Alphab(j);
		}
		_weightb[thread_id][thr] = event_.Event::Get_Weight();
		//std::cout<<"In Tree thread " <<thread_id <<" the weight is: " <<_weightb[thread_id] <<std::endl;
		if(event_.Event::Get_COM()){
			_COMb[thread_id][thr] = 1;
		}else{
			_COMb[thread_id][thr] = 0; 
		}
		if(thrown_){
			_athrown_tree[thread_id]->TTree::Fill();
		}else{
			_a_tree[thread_id]->TTree::Fill();
		}
		
		//std::cout<<"   successfully filled?"; 
	//}
}

void forest::scan_thread_tree(int thread_id, bool thrown_ ){
	if(thrown_){
		_athrown_tree[thread_id]->TTree::Scan();
	}else{
		_a_tree[thread_id]->TTree::Scan();
	}
}

void forest::fill_evnt(int event_n, int thread_id, bool thrown_){
	if(thrown_){
		_evnt[1] = event_n;
		_athrown_tree[thread_id]->TTree::Fill();
	}else{
		_evnt[0] = event_n;
		_a_tree[thread_id]->TTree::Fill();
	}
	//_a_tree[thread_id]->TTree::Scan();
}

void forest::fill_apart(int apart_n, int thread_id, bool thrown_){
	if(thrown_){
		_apart[1] = apart_n;
		_athrown_tree[thread_id]->TTree::Fill();
	}else{
		_apart[0] = apart_n;
		_a_tree[thread_id]->TTree::Fill();
	}
}
void forest::fill_px(float px_n, int i, int thread_id, bool thrown_){
	if(thrown_){
		_px[1][i] = px_n; 
		_athrown_tree[thread_id]->TTree::Fill();
	}else{
		_px[0][i] = px_n; 
		_a_tree[thread_id]->TTree::Fill();
	}
}
void forest::fill_py(float py_n, int i, int thread_id, bool thrown_){
	if(thrown_){
		_py[1][i] = py_n;
		_athrown_tree[thread_id]->TTree::Fill(); 
	}else{
		_py[0][i] = py_n;
		_a_tree[thread_id]->TTree::Fill();
	}
}

void forest::fill_pz(float pz_n, int i, int thread_id, bool thrown_){
	if(thrown_){
		_pz[1][i] = pz_n;
		_athrown_tree[thread_id]->TTree::Fill();
	}else{
		_pz[0][i] = pz_n;
		_a_tree[thread_id]->TTree::Fill();
	}
}
void forest::fill_p0(float p0_n, int i, int thread_id, bool thrown_){
	if(thrown_){
		_p0[1][i] = p0_n;
		_athrown_tree[thread_id]->TTree::Fill();
	}else{
		_p0[0][i] = p0_n;
		_a_tree[thread_id]->TTree::Fill();
	}
}
void forest::fill_pid(int pid_n, int i, int thread_id, bool thrown_){
	if(thrown_){
		_pid[1][i] = pid_n;
		_athrown_tree[thread_id]->TTree::Fill();
	}else{
		_pid[0][i] = pid_n;
		_a_tree[thread_id]->TTree::Fill();
	}
}
void forest::fill_hel(int hel_n, int thread_id, bool thrown_){
	if(thrown_){
		_hel[1] = hel_n;
		_athrown_tree[thread_id]->TTree::Fill();
	}else{
		_hel[0] = hel_n;
		_a_tree[thread_id]->TTree::Fill();
	}
}
void forest::fill_top(int top_n, int thread_id, bool thrown_){
	if(thrown_){
		_top[1] = top_n; 
		_athrown_tree[thread_id]->TTree::Fill();
	}else{
		_top[0] = top_n; 
		_a_tree[thread_id]->TTree::Fill();
	}
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

void forest::Grow_Forest(bool thrown_){
	std::cout<<"Growing a Forest\n";
	TList the_forest; 
	TList thr_forest;
	int num_ev[NUM_THREADS];
	int num_th[NUM_THREADS];
	for(int i = 0; i<NUM_THREADS ; i++){
		the_forest.TList::Add(_a_tree[i]);
	}
	_the_tree = TTree::MergeTrees(&the_forest,"");
	std::cout<<"Recon Forest Grown\n";
	if(thrown_){
		for(int j = 0; j<NUM_THREADS ; j++){
			thr_forest.TList::Add(_athrown_tree[j]);
		}
		_thrown_tree = TTree::MergeTrees(&thr_forest,"");
		std::cout<<"Thrown Forest Grown\n";
	}
}

void forest::Grow_Write_Forest(std::string tree_file_name, bool thrown_){
	forest::mkfile(tree_file_name, false);
	std::cout<<"Growing and Writing Forest:\n";
	TList the_forest; 
	TList thr_forest;
	int num_ev[NUM_THREADS];
	int num_th[NUM_THREADS];
	_tree_file->cd();
	TDirectory* _almanac1 = _tree_file->mkdir("r10");
	_almanac1->cd();
	for(int i = 0; i<NUM_THREADS ; i++){
		//the_forest->TList::Add(_a_tree[i]);
		the_forest.TList::Add(_a_tree[i]);
	}
	_the_tree = TTree::MergeTrees(&the_forest,"");
	//_new_tree1 = TTree::MergeTrees(_the_forest);
	std::cout<<"	Reconstructed Forest Grown\n";
	_the_tree->Write();
	_tree_file->Close();
	//_new_tree1->Write();
	std::cout<<"	Reconstructed Tree Written\n";
	if(thrown_){
		forest::mkfile(tree_file_name, thrown_);
		TDirectory* _almanac2 = _thr_tree_file->mkdir("t10");
		_almanac2->cd();
		for(int j = 0; j<NUM_THREADS ; j++){
			//thr_forest->TList::Add(_athrown_tree[j]);
			thr_forest.TList::Add(_athrown_tree[j]);
		}
		_thrown_tree = TTree::MergeTrees(&thr_forest,"");
		//_new_tree2 = TTree::MergeTrees(_thr_forest);
		std::cout<<"	Thrown Forest Grown\n";
		_thrown_tree->Write();
		//_new_tree2->Write();
		std::cout<<"	Thrown Tree Written\n";
		_thr_tree_file->Close();
	}
}


	/*
void forest::mkforest(){
  _the_tree = TTree::MergeTrees(_the_forest);
  _the_tree->TTree::Write();
}*/
