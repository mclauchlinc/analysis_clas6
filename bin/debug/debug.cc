#include <TTree.h>
#include <TFile.h>
#include <iostream>

int main() {
	const char pathname[]="/Users/cmc/Desktop/analysis/analysis_clas6/bin/09_06_2019_04_evnt_tree.root";
	TFile infile(pathname,"READ");
	TTree* tree;
	infile.GetObject("physics_phorest/TREEE",tree);
	std::cout << tree->GetEntries() << std::endl;
	return 0;
}