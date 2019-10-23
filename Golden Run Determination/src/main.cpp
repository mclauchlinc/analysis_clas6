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
	std::vector<std::vector<std::string>> infilenames2(NUM_THREADS);//For two lists of plate stuff

	int data_set = 0; //{e1-6,e1f, e1-6 sim, e1f sim} ->{1,2,3,4}
	int _case = 0; //For file entry
	int plate_stat = 0;//Status of the half wave plate
	bool hel_runs = false;

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
       	_case = 1; 
       	if(argv[1]!=list3h){
       		for(int i = 0; i <NUM_THREADS; i++){
	       		infilenames[i] = fun::read_file_list(filepath_map[argv[1]],i);//constants.hpp for map, functions.hpp for function
	       	}
	        if(argv[1] == list3p){
	        	std::cout<<"Half Wave Plate In" <<std::endl;
	       		plate_stat = 1; 
	        }else if(argv[1] == list3n){
	        	std::cout<<"Half Wave Plate Out" <<std::endl;
	       		plate_stat = -1; 
	        }
       	}else{
       		hel_runs = true;
       		std::cout<<"Looking at Both Plates" <<std::endl;
       		for(int i = 0; i <NUM_THREADS; i++){
	       		infilenames[i] = fun::read_file_list(filepath_map[list3p],i);//constants.hpp for map, functions.hpp for function
	       		infilenames2[i] = fun::read_file_list(filepath_map[list3n],i);//constants.hpp for map, functions.hpp for function
	       		plate_stat = 1;
	       	}
       	}
	}
	//Case 2
	if(argc >4){
		output_name = argv[1];
		data_set = std::atoi(argv[2]);
		for(int w = 3; w <argc; w++){//Putting all the files into 
			infilenames[(w-3)%NUM_THREADS].push_back(argv[w]);
		}
		_case = 2; 
	}
	std::cout<<"Brought Data in though Case: " <<_case <<std::endl; 

	// Make a set of threads (Futures are special threads which return a value)
  	std::future<size_t> threads[NUM_THREADS];
  	std::future<size_t> threads2[NUM_THREADS];

	//Define variable for total number of events
	size_t events = 0; 

	//Make histograms objects as a shared pointer that all threads will have
	auto hists = std::make_shared<Histogram>(output_name);//Check on this

	//Make relevant TTrees and Event Rootfile
	auto a_good_forest = std::make_shared<forest>(1); 
	a_good_forest->forest::mkfile(output_name);//Making the Tree File
	

	for(int i = 0; i<NUM_THREADS; i++){
		num_mixed_p_pip[i]=0;
	}

	//For each thread
	for(int i = 0; i<NUM_THREADS; i++){
		//Set the thread to run asynchronously
		//The function running is the first argument
		//The functions arguments are all remaining arguments
		threads[i] = std::async(run_files, infilenames.at(i), filepath_map[list3p], hists, a_good_forest, i, data_set, file_num, _case, plate_stat);//, num_mixed_p_pip[i]);
		//run_files(infilenames.at(i), filepath_map[argv[1]], hists, a_good_forest, i, data_set, file_num, _case);
		

		 //This is running the analysis 
	}
	for(int j = 0; j<NUM_THREADS; j++){
		threads[j].wait(); 
	}
	if(hel_runs){
		//For each thread for other plate config
		for(int i = 0; i<NUM_THREADS; i++){
			//Set the thread to run asynchronously
			//The function running is the first argument
			//The functions arguments are all remaining arguments
			threads2[i] = std::async(run_files, infilenames2.at(i), filepath_map[list3n], hists, a_good_forest, i, data_set, file_num, _case, -plate_stat);//, num_mixed_p_pip[i]);
			//run_files(infilenames.at(i), filepath_map[argv[1]], hists, a_good_forest, i, data_set, file_num, _case);
			

			 //This is running the analysis 
		}
		for(int j = 0; j<NUM_THREADS; j++){
			threads2[j].wait(); 
		}
	}
	
	a_good_forest->forest::Grow_Forest();//Combine all those Thread specific trees and output a root file with all event selected events with four vectors 

	//For each thread to see how many events each thread successfully analyized
	for(int i = 0; i<NUM_THREADS; i++){
		events += threads[i].get();
		events += threads2[i].get();
		//std::cout<<std::endl<<"Number of misIDed proton/pip particles: " <<num_mixed_p_pip[i] <<std::endl;
	}
	std::cout<<std::endl;
	//Timer and Efficiency Counters
	std::cout.imbue(std::locale(""));//Puts commas in appropriately? Apparently?
	std::chrono::duration<double> elapsed_full = (std::chrono::high_resolution_clock::now()-start);
	std::cout<< elapsed_full.count() << " Sec" <<std::endl;//Total elapsed time
	std::cout<< events/elapsed_full.count() <<" Hz" <<std::endl; 
	return 0; 

}