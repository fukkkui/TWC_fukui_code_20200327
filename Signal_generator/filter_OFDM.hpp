#include <complex>
#include <vector>
#include <iostream>

using COM = std::complex<double>;
using VCOM = std::vector<COM>;
using VVCOM = std::vector<VCOM>;
using V = std::vector<double>;
using VV = std::vector<V>;

class SpecialFilter{
public:
	V fil;
	SpecialFilter();
	void CP_Lap(VCOM &x, int N, int NN, int Ncp);
	void Make_filter(int Ncp);
	void CP_Lap_filtering(VCOM &x, int N, int NN, int Ncp);
};

SpecialFilter::SpecialFilter(){

}

void SpecialFilter::CP_Lap(VCOM &x, int N, int NN, int Ncp){
	VCOM output(NN);
	double coe;
	for(int i = 0; i < NN/(N + Ncp) - 1; i++){
		for(int j = N; j < N + Ncp; j++){
			coe = (double)(j - N) / (double)Ncp;
			x[j + i * (N + Ncp)] = x[j + i * (N + Ncp)] * (1.0 - coe) + x[j - Ncp + (i + 1) * (N + Ncp)] * coe;
		}
	}
}

void SpecialFilter::Make_filter(int Ncp){
	int L = Ncp;
	double Ts = Ncp;
	double a = 0.1;
	double pi = std::acos(-1.0);
	double val;
	fil.resize(L);
	fil[0] = 1.0;
	for(int i = 1; i < Ncp; i++){
		val = pi * (double)i / Ts;
		if(i%(int)Ts == 0){
			fil[i] = 0.0;
		}else{
			fil[i] = (std::sin(val)/(val)) * (std::cos(val * a)/(1 - pow(val * 2.0 * a / pi, 2)));
		}
	}
	// plot_V(fil);
}

void SpecialFilter::CP_Lap_filtering(VCOM &x, int N, int NN, int Ncp){
	Make_filter(Ncp);
	VCOM output(NN);
	for(int i = 0; i < NN/(N + Ncp) - 1; i++){
		for(int j = N; j < N + Ncp; j++){
			x[j + i * (N + Ncp)] = x[j + i * (N + Ncp)] * fil[j - N] + x[j - Ncp + (i + 1) * (N + Ncp)] * fil[fil.size() - 1 - (j - N)];
		}
	}
}