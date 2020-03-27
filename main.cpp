#include <iostream>
#include <vector>
#include <string>
#include "tool/file_save.hpp"
#include "submain.hpp"

//$ g++ main.cpp -std=c++0x -lopenblas
//$ ./a.out


int main(){
	time_t start = time(NULL);

	std::vector<int> P{1,3,5,7};

	int K_dB = 30;
	// std::vector<int> K_dB;
	// for(int i = -20; i < 51; i+=5){
	// 	K_dB.push_back(i);
	// }

	int IRR_dB = 30;
	// std::vector<int> IRR_dB;
	// for(int i = 0; i < 51; i+=5){
	// 	IRR_dB.push_back(i);
	// }

	int INR_dB = 80;
	// std::vector<int> INR_dB;
	// for(int i = 0; i < 101; i+= 5){
	// 	INR_dB.push_back(i);
	// }

	double _3dBBW = 50;
	// std::vector<double> _3dBBW;
	// for(double i = -3.0; i < 3.1; i += 0.1){
	// 	_3dBBW.push_back(pow(10, i));
	// }
	
	// int BO = 10;
	std::vector<int> BO;
	for(int i = 0; i < 21; i++){
		BO.push_back(i);
	}

	// int Nsym = 100;
	int Nsym = 500;
	// std::vector<int> Nsym;
	// for(int i = 30; i < 101; i+=5){
	// 	Nsym.push_back(i);
	// }

	double tau_tx = 2.0e-9;
	// std::vector<double> tau_tx;
	// double val = 0;
	// for(double i = -3.0; i < 1.7; i += 0.1){
	// 	tau_tx.push_back(pow(10, i) * 1e-9);
	// }
	// tau_tx.push_back(49.999e-9);

	// int L = K_dB.size();
	// int L = IRR_dB.size();
	// int L = INR_dB.size();
	// int L = _3dBBW.size();
	int L = BO.size();
	// int L = Nsym.size();
	// int L = tau_tx.size();
	// std::cout << L << std::endl;

	std::vector<std::vector<double> > can(P.size(), std::vector<double>(L));
	for(int i = 0; i < L; i++){
		// std::cout << "K = " << K_dB[i] << " dB" << std::endl;
		// std::cout << "IRR = " << IRR_dB[i] << " dB" << std::endl;
		// std::cout << "INR = " << INR_dB[i] << " dB" << std::endl;
		// std::cout << "3dBBW = " << _3dBBW[i] << " Hz" << std::endl;
		std::cout << "BO = " << BO[i] << " dB" << std::endl;
		// std::cout << "Nsym = " << Nsym[i] << " こ" << std::endl;
		// std::cout << "tau tx = " << tau_tx[i]/1e-9 << " nsec" << std::endl;
		for(int j = 0; j < P.size(); j++){
			std::cout << "P = " << P[j] << std::endl;
			// can[j][i] = submain(P[j], K_dB[i], IRR_dB, INR_dB, _3dBBW, BO, Nsym, tau_tx);
			// can[j][i] = submain(P[j], K_dB, IRR_dB[i], INR_dB, _3dBBW, BO, Nsym, tau_tx);
			// can[j][i] = submain(P[j], K_dB, IRR_dB, INR_dB[i], _3dBBW, BO, Nsym, tau_tx);
			// can[j][i] = submain(P[j], K_dB, IRR_dB, INR_dB, _3dBBW[i], BO, Nsym, tau_tx);
			can[j][i] = submain(P[j], K_dB, IRR_dB, INR_dB, _3dBBW, BO[i], Nsym, tau_tx);
			// can[j][i] = submain(P[j], K_dB, IRR_dB, INR_dB, _3dBBW, BO, Nsym[i], tau_tx);
			// can[j][i] = submain(P[j], K_dB, IRR_dB, INR_dB, _3dBBW, BO, Nsym, tau_tx[i]);
		}
		std::cout << std::endl;
	}

	std::vector<std::vector<double> > file_power(P.size() + 1, std::vector<double>(L));
	for(int i = 0; i < L; i++){
		// file_power[0][i] = K_dB[i];
		// file_power[0][i] = IRR_dB[i];
		// file_power[0][i] = INR_dB[i];
		// file_power[0][i] = _3dBBW[i];
		file_power[0][i] = BO[i];
		// file_power[0][i] = Nsym[i];
		// file_power[0][i] = tau_tx[i]/1e-9;
		for(int j = 0; j < P.size(); j++){
			file_power[j + 1][i] = can[j][i];
		}
	}
	// std::string file = "K.csv";
	// std::string file = "IRR.csv";
	// std::string file = "INR.csv";
	// std::string file = "3dBBW.csv";
	std::string file = "BO.csv";
	// std::string file = "Nsym.csv";
	// std::string file = "tau.csv";
	Multiple_fileW_csv(file, file_power);

	time_t end = time(NULL);
	std::cout << "Duration = " << end - start << " sec." << std::endl;
}