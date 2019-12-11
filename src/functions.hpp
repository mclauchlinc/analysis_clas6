#ifndef FUNCTIONS_H_GUARD
#define FUNCTIONS_H_GUARD

#include "constants.hpp"
#include "TFile.h"
#include "TChain.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include "branches.hpp"
#include "environment.hpp"

namespace fun {
bool replace(std::string& str, const std::string& from, const std::string& to);

std::shared_ptr<TFile> Name_File(std::string a_file_name);
std::shared_ptr<TFile> Name_Tree_File(std::string a_file_name);

std::vector<std::string> read_file_list(std::string path, int thread_num);

void loadChain(std::shared_ptr<TChain> c, std::string file, int thread_id, int max);

char* appendCharToCharArray(char* array, char a);



}

#endif