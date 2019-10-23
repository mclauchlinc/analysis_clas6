#ifndef MAIN_H_GUARD
#define MAIN_H_GUARD

#include <iostream>
#include "TFile.h"
#include "TH1.h"
#include "event_class.hpp"
#include "histogram.hpp" 
#include "constants.hpp"
#include "functions.hpp"
#include "branches.hpp"
#include "ntuple.hpp"

std::string comp; //Variable for choosing which data set will be used
char* output_name;//Variable for the output file name. This is reassigned through input parameters
int file_num = -1;//The initial assignment for the number of files in the program. -1 equates to all of them




size_t run(std::shared_ptr<TChain> _chain, std::shared_ptr<Histogram> _hists, int run_type){//, int &num_ppip){
	//Number of events in this thread
	size_t num_of_events = (int) _chain->GetEntries();

	//Make a data object which all the branches can be accessed from
	auto data = std::make_shared<Branches>(_chain);

	//Analysis loop
	for(size_t curr_event = 0; curr_event < num_of_events; curr_event++){
		//Get singular event
		_chain->GetEntry(curr_event);
		//For Progress just look at the 0th thread
		std::cout<<"\t" <<(100*curr_event/num_of_events) <<" %\r" <<std::flush;

		//Analyze the event and look for particle ID and topology 
		auto event = std::make_shared<Event>(data,_hists,run_type,plate_);//All event selection happens in here 
		//num_ppip += event->Event_Class::Get_ppip();
	}
}


size_t run_files(std::vector<std::string> inputs, std::string list_file, std::shared_ptr<Histogram> hists, std::shared_ptr<forest> bforest, int thread_id, int run_type, int max, int _case, int plate_){//, int &num_ppip){
	//Called once per thread
	//Make a new chain to process for this thread
	auto chain = std::make_shared<TChain>("h10");
	//Add every file to the chain
	if(_case==1){
		fun::loadChain(chain,list_file,thread_id,max);//Splits the list of files up into the individual thread chains
	}else if(_case==2){
		for(auto in:inputs) chain->Add(in.c_str());//This will have already split up the list of input files into individual lists which become chains
	}else{
		std::cout<<"Not a proper case" <<std::endl; 
	}

	
	//Run the function over each thread
	return run(chain,hists,bforest,thread_id,run_type,plate_);//,num_ppip);
}

#endif