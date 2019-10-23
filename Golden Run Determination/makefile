#Get the proper includes for root
ROOTLIBS	= $(shell root-config --libs)
INCLUDE		= $(shell root-config --incdir)
ROOTSYS		= $(shell root-config --exec-prefix)
#Tell which compiler to use
CXX = g++
CXXFLAGS =      -O2 -fPIC -w -g -pthread -stdlib=libc++ -std=c++14 -m64 -I/usr/local/root/include
#Target is the name of the output executable
TARGET =	golden
#The name of the file you want to comilie usually something like main.cpp 
FILENAME =	main

SRC = $(wildcard src/*.cpp)
OBJ = $(patsubst src/%.cpp, obj/%.o, $(SRC))



#When you don't specifiy anything i.e. $ make the function all is called
#If you only want to clean the you can call the function clean with $ make clean
#if you only want to make the main program you can call the function $ make main
#all:	clean main

#main:	$(FILENAME).o 
#	$(CXX) $(FILENAME).o -L. $(CXXFLAGS) $(ROOTLIBS) -o $(TARGET) 

#clean:
#	-rm -f $(TARGET) $(FILENAME).o 

#This original bit was made by Nick, but proved to try to be too fancy. Gary helped me write this simplified version seen below
#Everything is explicit
#How to make analysis, what compiler to use etc. 

#To run a makefile with a name other than make 
#make -f <name of make file>
all: clean golden

golden: $(OBJ)
	$(CXX) $(CXXFLAGS)  $(ROOTLIBS) $(OBJ) -o bin/golden

obj/%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(TARGET) $(OBJ)

