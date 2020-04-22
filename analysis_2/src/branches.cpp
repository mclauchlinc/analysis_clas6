
#include "branches.hpp"




Branches::Branches(std::shared_ptr<TChain> tree) {
  _tree = tree;
  Branches::init();
}

void Branches::init() {
  _tree->SetBranchAddress("evnt",&_evnt);
  _tree->SetBranchAddress("apart",&_apart);
  _tree->SetBranchAddress("px",&_px);
  _tree->SetBranchAddress("py",&_py);
  _tree->SetBranchAddress("pz",&_pz);
  _tree->SetBranchAddress("p0",&_p0);
  _tree->SetBranchAddress("pid",&_pid);
  _tree->SetBranchAddress("hel",&_hel);
  _tree->SetBranchAddress("top",&_top);
  _tree->SetBranchAddress("fc_tot",&_fc_tot);
  _tree->SetBranchAddress("run_type",&_run_type);
}

int Branches::evnt() { return _evnt; }
int Branches::apart() { return _apart; }
float Branches::px(int i) { return _px[i]; }
float Branches::py(int i) { return _py[i]; }
float Branches::pz(int i) { return _py[i]; }
float Branches::p0(int i) { return _py[i]; }
TLorentzVector Branches::Par_4Vec(int i) {
      return physics::Make_4Vec(_px[i],_py[i],_pz[i],_p0[i]); 
    }
int Branches::hel() { return _hel; }
int Branches::top() { return _top; }
int Branches::fc_tot() { return _fc_tot; }
int run_type(){ return _run_type; }


