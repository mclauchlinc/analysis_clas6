#include "functions.hpp"

bool replace(std::string& str, const std::string& from, const std::string& to) {
    size_t start_pos = str.find(from);
    if(start_pos == std::string::npos)
        return false;
    str.replace(start_pos, from.length(), to);
    return true;
}

std::shared_ptr<TFile> Name_File(std::string a_file_name)
{
	std::string file_name = "$name.root";
	replace(file_name, "$name", a_file_name);
	return std::make_shared<TFile>(file_name.c_str(),"RECREATE");
}


std::vector<std::string> read_file_list(std::string path, int thread_num){
  std::ifstream infile(path.c_str()); // in file stream
  std::vector<std::string> result;
  std::string line;
  int t = 0;
  while(getline(infile,line)) { //getline sees if there is a line available
    if(thread_num == (t%NUM_THREADS)){
      result.push_back(line);//Gets the current line
    }
    t++;
  }
  return result;
}

void loadChain(std::shared_ptr<TChain> c, std::string file, int thread_id, int max)
{
  std::vector<std::string> filelist = read_file_list(file,thread_id);//read_file_list(file); //creates a vector of file names
  //If not specified will take in all the files in the text file
  int test = filelist.size();
  if(max > test)
  {
    std::cout<< "You tried to add too many files. This has been corrected" <<std::endl <<"Remember that you may only add " <<test <<" files" <<std::endl;
  }
  if(max == -1 || max > test) {//In case one tries to add too many files
    max = filelist.size();
  }
  //If specified then it will take in that number of files 
  for(int i = 0; i < max; i++) {
    c->AddFile(filelist[i].c_str());
  }
}


//Making a tree to contain all selected events
void mkttree(std::string tree_file_name, std::string tree_desc){
	Int_t evnt; //The #event for the given file
	Int_t apart;//The particle in each event  
	Float_t px; 
	Float_t py; 
	Float_t pz;
	Float_t p0;//Energy of the particle 
	Int_t pid;//particle ID 
	Int_t hel; //helicity
	Int_t top; //Topology {pmiss,pipmiss,pimmiss,zero} -> {1,2,3,4}

	std::make_shared<TFile> event_file(tree_file_name,"RECREATE");//The output rootfile containng the tree
	std::make_shared<TTree> event_tree = new TTree("t1",tree_desc);//
	t1->Branch("evnt",&evnt,"evnt/I");
	t1->Branch("apart",&apart"apart/I");
	t1->Branch("px",&evnt,"px[apart]/F");
	t1->Branch("py",&evnt,"py[apart]/F");
	t1->Branch("pz",&evnt,"pz[apart]/F");
	t1->Branch("p0",&evnt,"p0[apart]/F");
	t1->Branch("pid",&evnt,"pid[apart]/F");
	t1->Branch("hel",&hel,"hel/I");
	t1->Branch("top",&top,"top/I");
}

void output_event(std::shared_ptr<TTree> the_tree, std::shared_ptr<Event> _event, int event_n)
{
	//create a Tree file tree1.root
	//create the file, the Tree and a few branches
	evnt = event_n; 
	for(int i; i<4; i++){
		apart = i; 
		px[i] = _event.Get_px(i);
		py[i] = _event.Get_py(i);
		pz[i] = _event.Get_pz(i);
		p0[i] = _event.Get_p0(i);
		pid[i] = _event.Get_pid(i);
	}
	hel = _event.Get_hel();
	top = _event.Get_top();
	the_tree->Fill();

}