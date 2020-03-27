#include <iostream>
#include <vector>
#include <complex>
#include <thread>
#include <string>
#include "ope2.hpp"

using COM = std::complex<double>;
using VCOM = std::vector<COM>;
using VVCOM = std::vector<VCOM>;
using V = std::vector<double>;
using VV = std::vector<V>;

void pro(VV &PSD_bfr, VV &PSD_aft, int l, int N, int Nact, int NN, int Ncp, int P, double BW, bool PN, int K_dB, int IRR_dB, int INR_dB, double _3dBBW, int BO, double tau_tx){
	VV result(2, V(N));
	result = main_operation(l, N, Nact, NN, Ncp, P, BW, PN, K_dB, IRR_dB, INR_dB, _3dBBW, BO, tau_tx);
	for(int j = 0; j < N; j++){
		PSD_bfr[l][j] = result[0][j];
		PSD_aft[l][j] = result[1][j];
	}
}

double submain(int P, int K_dB, int IRR_dB, int INR_dB, double _3dBBW, int BO, int Nsym, double tau_tx){
	bool PN = true; //phase noise
	int loop = 200; //# of trials
	int N = 512; //FFT size
	int Ncp = N/4; //Cyclic prefix
	int Nact = 52; //active subcarrier
	int NN = (N + Ncp) * Nsym; //Total samples
	double BW = 20e6; // bandwidth
	VV PSD_bfr(loop, V(N,0)); //power spectral density before cancellation
	VV PSD_aft(loop, V(N,0)); //power spectral density after cancellation

    //Loop
    int m = std::max(std::thread::hardware_concurrency(), 1u);
    std::vector<std::thread> worker;
    for (int i = 0; i < m; i++){
    	worker.emplace_back([&](int id){
    		int r0 = loop/m * id + std::min(loop%m, id);
    		int r1 = loop/m * (id + 1) + std::min(loop%m, id + 1);
    		for(int j = r0; j < r1; j++){
    			pro(PSD_bfr, PSD_aft, j, N, Nact, NN, Ncp, P, BW, PN, K_dB, IRR_dB, INR_dB, _3dBBW, BO, tau_tx);
    		}
    	}, i);
    }
    for(auto& t: worker) t.join();


    // calculation of cancellation amount
    double Power_bfr = 0.0;
    double Power_aft = 0.0;
    double can = 0.0;
    for(int i = 0; i < loop; i++){
        for(int j = 0; j < N; j++){
            Power_bfr += PSD_bfr[i][j];
            Power_aft += PSD_aft[i][j];
        }
        can += Power_bfr/Power_aft;
    }
    can /= (double)loop;

    std::cout << "Cancellation = " << 10.0 * log10(can) << " dB" << std::endl;

    return 10.0 * log10(can);
}