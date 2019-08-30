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




size_t run(std::shared_ptr<TChain> _chain, std::shared_ptr<Histogram> _hists, int thread_id, int run_type){
	//Number of events in this thread
	size_t num_of_events = (int) _chain->GetEntries();
	//Print out information about the thread
	std::cout<<"Thread " <<thread_id <<": " <<num_of_events <<" Events" <<std::endl;

	//Make a data object which all the branches can be accessed from
	auto data = std::make_shared<Branches>(_chain);

	//Total number of events analyzed
	size_t total_com = 0; 
	size_t good_e = 0; 
	size_t good_event = 0; 

	//Analysis loop
	for(size_t curr_event = 0; curr_event<num_of_events; curr_event++){
		//Get singular event
		_chain->GetEntry(curr_event);
		//For Progress just look at the 0th thread
		if(thread_id == 0 && curr_event %1000 == 0){
			std::cout<<"\t" <<(100*curr_event/num_of_events) <<" %\r" <<std::flush;
		}

		total_com++; 
		//Sanity electron cuts
		//if(sanity_elec()) continue; *come back to this*

		//Make a reaction class from the data given
		//std::unique_ptr<Event> 
		auto event = std::make_unique<Event>(data,_hists,run_type);
		/*
		for(int part = 1; part < data->gpart(); part++){
			//Check Particle ID's and fill the reaction class

		}*/
		
	}
}


size_t run_files(std::vector<std::string> inputs, std::string list_file, std::shared_ptr<Histogram> hists, int thread_id, int run_type, int max, int _case){
	//Called once per thread
	//Make a new chain to process for this thread
	auto chain = std::make_shared<TChain>("h10");
	//Add every file to the chain
	if(_case==1){
		fun::loadChain(chain,list_file,thread_id,max);
	}else if(_case==2){
		for(auto in:inputs) chain->Add(in.c_str());
	}else{
		std::cout<<"Not a proper case" <<std::endl; 
	}
	
	
	//Run the function over each thread
	return run(chain,hists,thread_id,run_type);
}

#endif