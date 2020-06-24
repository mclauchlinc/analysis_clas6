#ifndef CARTESIANGENERATOR_HH
#define CARTESIANGENERATOR_HH

#include<vector>

// Cartesian product generator

class CartesianGenerator {	
   public:
   std::vector<long> space; //size of each dimension
   long index; //current combination index starting from 0
   std::vector<long> current; //current combination

   //Also for initialization
   void Reset() {
     //init index
   	  index = -1;
   	  //init current
   	  current.resize(space.size());
   	  for(int i = 0; i < space.size(); ++i) {
   	  	current[i] = 0;
   	  }
   	  current[0] = -1;
   }

   CartesianGenerator(const std::vector<long>& s) {
   	//init space
   	  space.resize(s.size());
   	  for(int i = 0; i < space.size(); ++i) {
   	  	space[i] = s[i];
   	  }
   	  Reset();
   }

   long operator[](int i) {
   	  return current[i];
   }

   bool GetNextCombination() {
   	  index++;
   	  bool success = false;
   	  int dim = 0;
   	  while(!success) {
   	  	//Check sanity
   	  	if(dim < current.size()) {
   	  	//Try to increment
   	  		current[dim]++;
   	  		if(current[dim] >= space[dim]) {
   	  			current[dim]=0;
   	  			dim++;
   	  		}
   	  		else success=true;
   	  	}
   	  	//Exhausted combinations
   	  	else return false;
   	  }
   	  //Successfully generated new combination
   	  return true;
   }
};

/*

//Setup space
  std::vector<long> space_dimensions(3);
  space_dimensions[0] = 6; //sector
  space_dimensions[1] = 4; //species
  space_dimensions[2] = 4; //cut status
  
//Create generator
  CartesianGenerator cart(space_dimensions);

  while(cart.GetNextCombination()) {
	long sector = cart[0] + 1;
	long species = cart[1];
	long cut_status = cart[2];
  }
*/

#endif