#pragma once

#include <complex>
#include <vector>
#include <iterator>
#include <random>


std::vector<std::complex<double> > Rician_fading(std::vector<double> &Amp_h){
	std::vector<std::complex<double> > h(Amp_h.size());

	std::random_device seed1;
	std::random_device seed2;
	std::default_random_engine generator1(seed1());
	std::default_random_engine generator2(seed2());
	h[0] = std::complex<double>(Amp_h[0], 0.0);

	for(int i = 1; i < h.size(); i++){
		std::normal_distribution<double> dist1(0, Amp_h[i]/sqrt(2.0));
		std::normal_distribution<double> dist2(0, Amp_h[i]/sqrt(2.0));
		h[i] = std::complex<double>(dist1(generator1), dist2(generator2));
	}

	return h;
}


