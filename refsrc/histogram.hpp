#ifndef HISTOGRAM_H_GUARD
#define HISTOGRAM_H_GUARD

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "constants.hpp"
#include "functions.hpp"
//#include "variables.h"
//#include "CartesianGenerator.hh"


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
	void WQ2_Write();
};





#endif