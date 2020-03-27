#pragma once

#include <vector>
#include <complex>
#include <fstream>
#include <string>
#include <sstream>
#include <iostream>


template <typename A>
void Multiple_fileW_csv(std::string filename, std::vector<std::vector<A> > &data){
    int Nd = data.size();
    int Ld = data.front().size();
    std::ofstream os(filename);
    for(int i = 0; i < Ld; i++){
        for(int j = 0; j < Nd; j++){
            os << data[j][i] << ",";
        }
        os << std::endl;
    }
    os.close();
}

