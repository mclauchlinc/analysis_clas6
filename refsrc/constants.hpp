#ifndef CONSTANTS_H_GUARD
#define CONSTANTS_H_GUARD

#include <unordered_map>
#include "TMath.h"


static const int MAX_PARTS = 100; 
static const int NUM_THREADS = 4;


//Fun names for file lists
static const std::string list1 = "one";
static const std::string list2 = "two";
static const std::string list3 = "three";
static const std::string list4 = "four";
static const std::string list3p = "three+";
static const std::string list3n = "three-";

//Paths for file names
static const std::string path1 = "/home/mclauchlinc/Desktop/analysis/nick_convert_e16.txt";
static const std::string path2 = "/Users/cmc/Desktop/analysis/e16_10_18_17_ntple.txt";
static const std::string path3 = "/Users/cmc/Desktop/analysis/Nick_skim_e16.txt";
static const std::string path4 = "/Users/cmc/Desktop/analysis/arjun_sim.txt";
static const std::string path3p = "/Users/cmc/Desktop/analysis/NickSkim_e16_PlateIN.txt";
static const std::string path3n = "/Users/cmc/Desktop/analysis/NickSkim_e16_PlateOUT.txt";

//Run type for file lists {e16,e1f,e16sim,e1fsim} -> {1,2,3,4}
static const int type1 = 1;
static const int type2 = 1;
static const int type3 = 1;
static const int type4 = 3;
static const int type3p = 1;
static const int type3n = 1; 

//Map of file lists to run type
static std::unordered_map<std::string, std::string> filepath_map = 	{{list1,path1},
																					{list2,path2},
																					{list3,path3},
																					{list4,path4},
																					{list3p,path3p},
																					{list3n,path3n}};

static std::unordered_map<std::string, int> filetype_map = 	{{list1,type1},
																			{list2,type2},
																			{list3,type3},
																			{list4,type4},
																			{list3p,type3p},
																			{list3n,type3n}};

#endif
