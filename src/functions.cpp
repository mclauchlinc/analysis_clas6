#include "functions.hpp"

bool fun::replace(std::string& str, const std::string& from, const std::string& to) {
    size_t start_pos = str.find(from);
    if(start_pos == std::string::npos)
        return false;
    str.replace(start_pos, from.length(), to);
    return true;
}

std::shared_ptr<TFile> fun::Name_File(std::string a_file_name)
{
	std::string file_name = "$name.root";
	replace(file_name, "$name", a_file_name);
	return std::make_shared<TFile>(file_name.c_str(),"RECREATE");
}

std::shared_ptr<TFile> fun::Name_Tree_File(std::string a_file_name, bool thrown_)
{
  std::string file_name;
  if(thrown_){
    file_name = "$name_thr_tree.root";
  }else{
    file_name = "$name_evnt_tree.root";
  }
	replace(file_name, "$name", a_file_name);
	return std::make_shared<TFile>(file_name.c_str(),"RECREATE");
}


std::vector<std::string> fun::read_file_list(std::string path, int thread_num){
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
void fun::removeTree(std::string file_name){

  TFile *file=new TFile((file_name).c_str(),"update");
  std::string object_to_remove="h10;1";
  //the object can be a tree, a histogram, etc, in this case "test1" is a TTree
  //notice the ";1" which means cycle 1; to remove all cycles do ";*"
  //if your object is not at the top directory, but in a directory in the .root file, called foo
  // you do first
  //file->cd("foo");
  //then continue with the Delete command which is only applied to the current gDirectory
  gDirectory->Delete(object_to_remove.c_str());
  file->Close();
}

void fun::loadChain(std::shared_ptr<TChain> c, std::string file, int thread_id, int max, int run_type)
{
  std::vector<std::string> filelist = fun::read_file_list(file,thread_id);//read_file_list(file); //creates a vector of file names
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
    //if(run_type ==3 || run_type == 4){//With some of the larger sim files this seems to be an issue where there are multiple trees in the sim files..?
    //  fun::removeTree(filelist[i]);
    //}
    c->AddFile(filelist[i].c_str());
  }
}

char* fun::appendCharToCharArray(char* array, char a)
{
    size_t len = strlen(array);

    char* ret = new char[len+2];

    strcpy(ret, array);    
    ret[len] = a;
    ret[len+1] = '\0';

    return ret;
}

bool fun::no_pro_pip_match(int idx1, int idx2[20]){//Designed to check to see if there was a double identification on a pip/proton
  bool pass = true;
  for(int i=0; i< 20; i++){
    if(idx1 == idx2[i]){
      pass = false;
    }
  }
  return pass; 
}

bool fun::hist_fitting(int species_, int cut_, int Wbin_, int pbin_, int fit_){
  bool pass = false; 
  //std::cout<<fit_;
  if(fit_ >= 0){
    if(fit_ == 1){
      if(species_==0){
        if(cut_ == 5 || cut_ == 6 || cut_ == 7){
          pass = true;
        }
      }else{
        if(cut_== 2 || cut_ == 3){
          pass = true;
        }
      }
    }else{
      pass = true;
    }
  }else if(fit_==-1 && Wbin_== 0 && pbin_ == 0){
    if(species_==0){
        if(cut_ != 5 && cut_ != 6 && cut_ != 7){
          pass = true;
        }
      }else{
        if(cut_!= 2 && cut_ != 3){
          pass = true;
        }
      }
  }
  return pass; 
}


