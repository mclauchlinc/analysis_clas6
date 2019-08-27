#ifndef READ_IN_DATA_H_GUARD
#define READ_IN_DATA_H_GUARD

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <TChain.h>
#include "variables.h"
#include "TTree.h"
#constants.h


//Creates a vector of file names from a source file labeled "path"

std::vector<std::string> read_file_list(std::string path, int thread_id) //Choose a vector of strings. Vector because adding things onto the end works fine. Also don't know how to use a list
{
	std::ifstream infile(path.c_str()); // in file stream
	std::vector<std::string> result;
	std::string line;
	int num= 0; 
	while(getline(infile,line)) { //getline sees if there is a line available
		if(num % )
		result.push_back(line);//Gets the current line
	}
	return result;
}


//Create a data chain for a given data set
//Optional number of files to load
void loadChain(TChain* c, std::string file, int max=-1)
{
	std::vector<std::string> filelist = read_file_list(file);//read_file_list(file); //creates a vector of file names
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

void SetBranches(TChain* c, int witch)
{
	if(witch == 3){
		c->SetBranchAddress("evntclas2",&evntclas2);
		c->SetBranchAddress("q",&q_c);
        c->SetBranchAddress("gpart",&gpart);
        c->SetBranchAddress("sc_t",&sc_t_c);
        c->SetBranchAddress("sc_r",&sc_r_c);
        c->SetBranchAddress("sc",&sc_c);
        c->SetBranchAddress("p",&p_c);
        c->SetBranchAddress("cx",&cx_c);
        c->SetBranchAddress("cy",&cy_c);
        c->SetBranchAddress("cz",&cz_c);
        c->SetBranchAddress("stat",&stat_c);
        c->SetBranchAddress("dc_stat",&dc_stat_c);
        c->SetBranchAddress("vx",&vx_c);
        c->SetBranchAddress("vy",&vy_c);
        c->SetBranchAddress("vz",&vz_c);
        c->SetBranchAddress("dc",&dc_c);
        c->SetBranchAddress("cc",&cc_c);
        c->SetBranchAddress("ec",&ec_c);
        c->SetBranchAddress("etot",&etot_c);
        c->SetBranchAddress("id",&id_c);
        c->SetBranchAddress("cc_hit",&cc_hit_c);
        c->SetBranchAddress("cc_sect",&cc_sect_c);
        c->SetBranchAddress("cc_segm",&cc_segm_c);
        c->SetBranchAddress("cc_part",&cc_part_c);
        c->SetBranchAddress("nphe",&nphe_c);
        c->SetBranchAddress("q_l",&q_l_c);
	}else{
		c->SetBranchAddress("evntclas2",&evntclas2);
		c->SetBranchAddress("q",&q_b);
	    c->SetBranchAddress("gpart",&gpart);
	    c->SetBranchAddress("sc_t",&sc_t_b);
	    c->SetBranchAddress("sc_r",&sc_r_b);
	    c->SetBranchAddress("sc",&sc_b);
	    c->SetBranchAddress("p",&p_b);
	    c->SetBranchAddress("cx",&cx_b);
	    c->SetBranchAddress("cy",&cy_b);
	    c->SetBranchAddress("cz",&cz_b);
	    c->SetBranchAddress("stat",&stat_b);
	    c->SetBranchAddress("dc_stat",&dc_stat_b);
	    c->SetBranchAddress("vx",&vx_b);
	    c->SetBranchAddress("vy",&vy_b);
	    c->SetBranchAddress("vz",&vz_b);
	    c->SetBranchAddress("dc",&dc_b);
	    c->SetBranchAddress("cc",&cc_b);
	    c->SetBranchAddress("ec",&ec_b);
	    c->SetBranchAddress("etot",&etot_b);
	    c->SetBranchAddress("id",&id_b);
	    c->SetBranchAddress("cc_hit",&cc_hit_b);
	    c->SetBranchAddress("cc_sect",&cc_sect_b);
	    c->SetBranchAddress("cc_segm",&cc_segm_b);
	    c->SetBranchAddress("cc_part",&cc_part_b);
	    c->SetBranchAddress("nphe",&nphe_b);
	    c->SetBranchAddress("q_l",&q_l_b);
	}
}


/*
TChain LoadData(char q, TChain* chain){
	if(q == 'd')
	{
		loadChain(&chain,"/home/mclauchlinc/Desktop/e16/nick.txt"); //This is specifically for the desktop
	}
	return data;
}
*/
#endif

/*
  TChain chain("h10");
  loadChain(&chain,"/some/path/files.txt");
  loadChain(&chain,"/some/path/files.txt",100);
*/