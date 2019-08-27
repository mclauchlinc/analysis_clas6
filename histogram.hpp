#ifndef HISTOGRAM_H_GUARD
#define HISTOGRAM_H_GUARD

#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "constants.hpp"
//#include "variables.h"
//#include "CartesianGenerator.hh"

//To allow naming for things without putting ".root" in the initialization
bool replace(std::string& str, const std::string& from, const std::string& to) {
    size_t start_pos = str.find(from);
    if(start_pos == std::string::npos)
        return false;
    str.replace(start_pos, from.length(), to);
    return true;
}

TFile *Name_File(std::string a_file_name)
{
	std::string file_name = "$name.root";
	replace(file_name, "$name", a_file_name);
	return new TFile(file_name.c_str(),"RECREATE");
}

using TH2D_ptr = std::shared_ptr<TH2D>;
using TH1D_ptr = std::shared_ptr<TH1D>;


class Histogram {
protected:
	std::shared_ptr<TFile> RootOutputFile;
	TCanvas* def; 
public:
	Histogram(const std::string& output_file);
	~Histogram();
	void Write();
};





#endif