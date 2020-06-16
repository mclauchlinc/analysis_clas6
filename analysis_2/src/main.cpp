#include "main.hpp"
#include <future>
#include <thread>
#include "TROOT.h"



//Main body
int main(int argc, char **argv){
	//Apparently one needs this to ensure R0t doesn't break with multiple threads
	ROOT::EnableThreadSafety();

	//start timer
	auto start = std::chrono::high_resolution_clock::now();


	int data_set = 0; //{e1-6,e1f, e1-6 sim, e1f sim} ->{1,2,3,4}

	/*
	This program is designed to take in TTrees from the first analysis program and spit out cross sections in the form of a THnSparse
	*/
	//Case 1
    comp = argv[1]; //variables.h
    output_name = argv[2]; //variables.h
    std::cout<<"Argument assignment complete" <<std::endl;

	// Make a set of threads (Futures are special threads which return a value)
  	std::future<size_t> threads[NUM_THREADS];
  	std::future<size_t> threads2[NUM_THREADS];

	//Define variable for total number of events
	size_t events = 0; 

	//Make histograms objects as a shared pointer that all threads will have
	auto hists = std::make_shared<Histogram>(output_name);//Check on this
	
	auto chain = std::make_shared<TChain>("h10");//Don't really need a chain, but it's what I know! 
	fun::loadChain(chain,fun::read_file_list(filepath_map[comp]));//Splits the list of files up into the individual thread chains

	run(chain,hists,filetype_map[comp]);

	
	std::cout<<std::endl;
	//Timer and Efficiency Counters
	std::cout.imbue(std::locale(""));//Puts commas in appropriately? Apparently?
	std::chrono::duration<double> elapsed_full = (std::chrono::high_resolution_clock::now()-start);
	std::cout<< elapsed_full.count() << " Sec" <<std::endl;//Total elapsed time
	std::cout<< events/elapsed_full.count() <<" Hz" <<std::endl; 
	return 0; 

}