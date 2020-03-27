#pragma once

#include <vector>
#include <complex>

using COM = std::complex<double>;
using V = std::vector<double>;
using VV = std::vector<V>;
using VCOM = std::vector<COM>;
using VVCOM = std::vector<VCOM>;

void make_qam(int size, VVCOM &qam){
	int n = sqrt(size);
	//round(x) -> xを四捨五入
	int qam_bits = round(log2(n));
	uint16_t qam_mask = ((1 << qam_bits) - 1);
	double scale = sqrt((2 * (pow(2, qam_bits * 2) - 1)) / 3);
	qam = VVCOM(n, VCOM(n));

 	for(uint16_t i = 0; i < n; i++){
    	for(uint16_t j = 0; j < n; j++){
			double real = j;
      		double imag = n - i - 1;

      		real = real * 2 - qam_mask;
      		imag = imag * 2 - qam_mask;
      		real = real / scale;
      		imag = imag / scale;

      		qam[i][j] = COM(real,imag);
    	}
  	}
}

void OFDM_qam(VCOM &input, int N, int Nact, int M){
	int n = sqrt(M);
	VVCOM qam;
	make_qam(M, qam);

	int Column;
	int Row;

	std::random_device seed;
	std::default_random_engine generator(seed());
	std::uniform_int_distribution<> dist(0, n-1);

 	input.resize(N);
	for(int i = 0; i < N; i++){
		input[i] = 0;
	}
	double val = 1.0/sqrt(2.0);
	for(int i = 1; i < Nact/2; i++){
		Column = dist(generator);
		Row = dist(generator);
		input[i] = qam[Column][Row];

		if(i < Nact/2){
			Column = dist(generator);
			Row = dist(generator);
			input[N - (Nact/2) + i] = qam[Column][Row];
		}
	}
}