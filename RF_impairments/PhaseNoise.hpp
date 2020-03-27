#pragma once

#include <vector>
#include <complex>
#include <random>

using COM = std::complex<double>;
using V = std::vector<double>;
using VCOM = std::vector<COM>;

//Free Running Oscillator
void FRO_txrxaux_3dB(V &Phase_tx, V &Phase_aux, V &Phase_rx, int N, double _3dBBW, double Sampling_Interval, double tau_tx, double tau_aux){
	std::random_device seed1;
	std::random_device seed2;
	std::random_device seed3;
	std::default_random_engine generator1(seed1());
	std::default_random_engine generator2(seed2());
	std::default_random_engine generator3(seed3());

	double variance1 = 4.0 * std::acos(-1.0) * _3dBBW * Sampling_Interval;
	double variance2 = 4.0 * std::acos(-1.0) * _3dBBW * tau_tx;
	double variance3 = 4.0 * std::acos(-1.0) * _3dBBW * std::abs(tau_aux - tau_tx);
	std::normal_distribution<double> dist1(0.0, sqrt(variance1));
	std::normal_distribution<double> dist2(0.0, sqrt(variance2));
	std::normal_distribution<double> dist3(0.0, sqrt(variance3));

	Phase_rx.resize(N);
	Phase_tx.resize(N);
	Phase_aux.resize(N);
	Phase_rx[0] = 0.0;
	Phase_tx[0] = dist2(generator2);
	Phase_aux[0] = Phase_tx[0] + dist3(generator3);
	for(int i = 1; i < N; i++){
		Phase_rx[i] = Phase_rx[i - 1] + dist1(generator1);
		Phase_tx[i] = Phase_rx[i] + dist2(generator2);
		Phase_aux[i] = Phase_tx[i] + dist3(generator3);
	}

	for(int i = 0; i < N; i++){
		Phase_rx[i] = - Phase_rx[i];
	}
}

void Make_PN(V &Phase, int N, VCOM &PhaseNoise){
	PhaseNoise.resize(N);
	for(int i = 0; i < N; i++){
		PhaseNoise[i] = COM(std::cos(Phase[i]), std::sin(Phase[i]));
	}
}