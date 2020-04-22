
#ifndef BRANCHES_H
#define BRANCHES_H
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include "TChain.h"
#include "constants.hpp"


//Create a data chain for a given data set
//Optional number of files to load

class Branches {
 private:
  std::shared_ptr<TChain> _tree;
  bool _MC = false;
  Int_t _evnt = -1; //The #event for the given file
  Int_t _apart = -1;//The particle in each event  
  Float_t _px[4] = {NAN,NAN,NAN,NAN}; 
  Float_t _py[4] = {NAN,NAN,NAN,NAN}; 
  Float_t _pz[4] = {NAN,NAN,NAN,NAN};
  Float_t _p0[4] = {NAN,NAN,NAN,NAN};//Energy of the particle 
  Int_t _pid[4] = {NAN,NAN,NAN,NAN};//particle ID 
  Int_t _hel = 2; //helicity
  Int_t _top = -1; //Topology {pmiss,pipmiss,pimmiss,zero} -> {1,2,3,4}
  Int_t _fc_tot = 0;
  Int_t _run_type = 0; ////{1,2,3,4}->{e16,e1f,e16sim,e1f sim} 


 public:
  Branches(std::shared_ptr<TChain> tree);
  ~Branches(){};
  int evnt();
  int apart();
  float px(int i);
  float py(int i);
  float pz(int i);
  float p0(int i);
  TLorentzVector Par_4Vec(int i);
  int hel();
  int top();
  int fc_tot();
  int run_type(); 
};

#endif
