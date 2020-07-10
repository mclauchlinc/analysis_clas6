#ifndef SETUP_H_GUARD
#define SETUP_H_GUARD

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <ctime>
#include "functions.hpp"
#include "environment.hpp"


namespace Setup{
void set_envi(std::shared_ptr<Environment> setup, int run_type, int fit = 0);

void make_envi_file(const std::string& output_name, std::shared_ptr<Environment> envi);
}

#endif