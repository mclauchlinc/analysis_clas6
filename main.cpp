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

	//for threading
	std::vector<std::vector<std::string>> infilenames(NUM_THREADS);

	int data_set = 0; //{e1-6,e1f, e1-6 sim, e1f sim} ->{1,2,3,4}

	/*
	When running this analysis you have several options
	Case 1: ./analysis [file set] [num of files] [outputfile name]
	Case 2: ./analysis [output filename] [run type] [input file names]
	//Note, case 2 will only occur for inputing two or more files
	*/
	//Case 1
	if(argc ==4){
        comp = argv[1]; //variables.h
        file_num = std::atoi(argv[2]); //variables.h
        output_name = argv[3]; //variables.h
        std::cout<<"Argument assignment complete" <<std::endl;
       	//Should make a mapping of input paths to data types
        for(int i = 0; i <NUM_THREADS; i++){
       		infilenames[i] = read_file_list(filepath_map[argv[1]],i);//constants.hpp
       }
	}
	//Case 2
	if(argc >4){
		output_name = argv[1];
		data_set = std::atoi(argv[2]);
		for(int w = 3; w <argc; w++){//Putting all the files into 
			infilenames[w%NUM_THREADS].push_back(argv[w]);
		}
	}

	// Make a set of threads (Futures are special threads which return a value)
  	std::future<size_t> threads[NUM_THREADS];

	//Define variable for total number of events
	size_t events = 0; 

	//Make histograms objects as a shared pointer that all threads will have
	auto hists = std::make_shared<Histogram>(output_name);//Check on this

	//For each thread
	for(int i = 0; i<NUM_THREADS; i++){
		//Set the thread to run a tas asynchronously
		//The function running is the first argument
		//The functions arguments are all remaining arguments
		threads[i] = std::async(run_files, infilenames.at(i), hists, i, data_set); //This is running the analysis 
	}

	//For each thread to see how many events each thread successfully analyized
	for(int i = 0; i<NUM_THREADS; i++){
		events += threads[i].get();
	}

	//Timer and Efficiency Counters
	std::cout.imbue(std::locale(""));//Puts commas in appropriately? Apparently?
	std::chrono::duration<double> elapsed_full = (std::chrono::high_resolution_clock::now()-start);
	std::cout<< elapsed_full.count() << " Sec" <<std::endl;
	std::cout<< events/elapsed_full.count() <<" Hz" <<std::endl; 
	return 0; 

}