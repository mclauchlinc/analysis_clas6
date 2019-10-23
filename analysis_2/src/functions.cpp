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

std::shared_ptr<TFile> fun::Name_Tree_File(std::string a_file_name)
{
	std::string file_name = "$name_evnt_tree.root";
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

void fun::loadChain(std::shared_ptr<TChain> c, std::string file)
{
  c->AddFile(file.c_str());
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

