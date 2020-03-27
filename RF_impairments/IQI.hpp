#pragma once

#include <vector>
#include <complex>
#include <random>
#include <iostream>

using COM = std::complex<double>;
using VCOM = std::vector<COM>;

VCOM Make_IQIB_random_const(int IRR_dB){
	std::random_device seed;
	std::default_random_engine generator(seed());

	double IRR = pow(10,-(double)IRR_dB / 20.0);
	double abs_coe1 = 1.0;
	double abs_coe2 = abs_coe1 * IRR;

	std::uniform_real_distribution<double> dist(0.0, 255.0);
	double rot2 = dist(generator);
	double phase2 = 2.0 * std::acos(-1.0) * rot2 / 256.0;
	COM eta2 = abs_coe2 * COM(std::cos(phase2), std::sin(phase2));
	
	VCOM eta(2);
	eta[0] = COM(abs_coe1, 0.0);
	eta[1] = eta2;

	return eta;
}