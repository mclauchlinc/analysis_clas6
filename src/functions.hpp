#ifndef FUNCTIONS_H_GUARD
#define FUNCTIONS_H_GUARD

#include "constants.hpp"
#include "TFile.h"
#include "TChain.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include "TROOT.h"
#include "TTree.h"
#include "TBrowser.h"
#include "event_class.hpp"


bool replace(std::string& str, const std::string& from, const std::string& to);

std::shared_ptr<TFile> Name_File(std::string a_file_name);


std::vector<std::string> read_file_list(std::string path, int thread_num);

void loadChain(std::shared_ptr<TChain> c, std::string file, int thread_id, int max);


//Making a tree to contain all selected events
void mkttree(std::string tree_file_name, std::string tree_desc);

//Outputs selected events to the TTree
void output_event(std::shared_ptr<TTree> the_tree, std::unique_ptr<Event> _event, int event_n);


#endif