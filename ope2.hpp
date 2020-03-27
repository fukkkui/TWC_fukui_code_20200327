#include <iostream>
#include <vector>
#include <complex>
#include <stdlib.h> //for rand()
#include "RF_impairments/IQI.hpp"
#include "RF_impairments/AWGN.hpp"
#include "RF_impairments/channel.hpp"
#include "RF_impairments/amplifier.hpp"
#include "RF_impairments/PhaseNoise.hpp"
#include "Signal_generator/OFDM.hpp"
#include "Signal_generator/filter_OFDM.hpp"
#include "tool/fft.hpp"
#include "estimation.hpp"

using COM = std::complex<double>;
using VCOM = std::vector<COM>;
using VVCOM = std::vector<VCOM>;
using V = std::vector<double>;
using VV = std::vector<V>;

VV main_operation(int loop, int Ns, int Nact, int NN, int Ncp, int P, double Bandwidth, bool PN, int K_dB, int IRR_dB, int INR_dB, double _3dBBW, int BO, double tau_tx){
	//	------------------					
	//
	//	Preparation
	//
	//  ------------------


	srand((unsigned) time(NULL));

	//IQ Imbalance
	VCOM eta_tx, eta_aux, eta_rx;
	eta_tx = Make_IQIB_random_const(IRR_dB);
	eta_aux = Make_IQIB_random_const(IRR_dB);
	eta_rx = Make_IQIB_random_const(IRR_dB);


	//AWGN generation
	double PA_max = 0.0;
	double Floor = - INR_dB + (PA_max - BO);
	VCOM Noise(NN, 0);
	AWGN_no(Floor, Noise);

	//Phase noise generation
	V PN_tx(NN);
	V PN_rx(NN);
	V PN_aux(NN);
	VCOM PNR_tx(NN);
	VCOM PNR_rx(NN);
	VCOM PNR_aux(NN);
	if(PN == true){
		double Ts = 1.0 / Bandwidth;
		double tau_aux = 1.0e-9;		
		FRO_txrxaux_3dB(PN_tx, PN_aux, PN_rx, NN, _3dBBW, Ts, tau_tx, tau_aux);
		Make_PN(PN_rx, NN, PNR_rx);
		Make_PN(PN_tx, NN, PNR_tx);
		Make_PN(PN_aux, NN, PNR_aux);
	}
	if(PN == false){
		for(int i = 0; i < NN; i++){
			PNR_tx[i] = COM(1.0, 0.0);
			PNR_rx[i] = COM(1.0, 0.0);
			PNR_aux[i] = COM(1.0, 0.0);
		}
	}

	//Amp RX
	double smoo = 3.0;
	Amplifier amp_rx;
	V amp_coeff_rx(3, 0.0);
	double LNA_gain_rx = 0.0;
	amp_coeff_rx[0] = pow(10, LNA_gain_rx/20.0);//small signal gain
	double LNA_max_rx = 0.0;
	amp_coeff_rx[1] = pow(10, LNA_max_rx/20.0);//saturation power
	amp_coeff_rx[2] = smoo;//smoothness factor
	amp_rx.set_coeff(amp_coeff_rx);

	//Amp TX
	Amplifier amp_tx;
	V amp_coeff_tx(3, 0.0);
	double PA_gain = 0.0;
	amp_coeff_tx[0] = pow(10, PA_gain/20.0);//small signal gain
	amp_coeff_tx[1] = pow(10, PA_max/20.0);//saturation power
	amp_coeff_tx[2] = smoo;//smoothness factor
	amp_tx.set_coeff(amp_coeff_tx);

	//Signal generation
	VCOM wave_tx(NN,0);
	VCOM wave_aux(NN,0);
	VCOM wave_rx(NN,0);
	VCOM wave_part;
	for(int i = 0; i < NN/(Ns + Ncp); i++){
		OFDM_qam(wave_part, Ns, Nact, 16);
		fft_COMB(Ns,(double)2.0*std::acos(-1.0)/Ns, wave_part);
		for(int j = 0; j < Ns; j++){
			wave_part[j] /= (double)Ns;
		}
		for(int j = 0; j < Ns + Ncp; j++){
			wave_tx[j + i * (Ns + Ncp)] = wave_part[j%Ns];
		}
	}
	for(int i = 0; i < NN/(Ns + Ncp); i++){
		OFDM_qam(wave_part, Ns, Nact, 16);
		fft_COMB(Ns,(double)2.0*std::acos(-1.0)/Ns, wave_part);
		for(int j = 0; j < Ns; j++){
			wave_part[j] /= (double)Ns;
		}
		for(int j = 0; j < Ns + Ncp; j++){
			wave_aux[j + i * (Ns + Ncp)] = wave_part[j%Ns];
		}
	}
	SpecialFilter a;
	a.CP_Lap_filtering(wave_tx, Ns, NN, Ncp); //transition goes smoothly
	a.CP_Lap_filtering(wave_aux, Ns, NN, Ncp);

	//Power adjustment
	double ave_power_SI = 0.0;
	for(int i = 0; i < NN; i++){
		ave_power_SI += pow(abs(wave_tx[i]), 2);
	}
	ave_power_SI /= (double)NN;
	double fil_send_power = 0.0; //power when filtering
	double send_power2 = PA_max - (double)BO; //PA input power
	double fil_Amp_coe = pow(10, fil_send_power/20.0) / sqrt(ave_power_SI); //coefficient
	double Amp_coe2 = pow(10, (send_power2 - fil_send_power)/20.0); //coefficient
	VCOM wave_tx_cp(NN); //copy for estimation
	VCOM wave_aux_cp(NN); //copy for estimation
	for(int i = 0; i < NN; i++){
		wave_tx[i] = wave_tx[i] * fil_Amp_coe;
		wave_aux[i] = wave_aux[i] * fil_Amp_coe;
		wave_tx_cp[i] = wave_tx[i];
		wave_aux_cp[i] = wave_aux[i];
	}

	//channel generation
	COM haux = COM(std::cos(2.0 * std::acos(-1.0) * 3.0 / 1000.0), std::sin(2.0 * std::acos(-1.0) * 3.0 / 1000.0)); //random
	double K = pow(10, K_dB / 10.0); //K factor
	double Ohm_dB = 0.0; //channel power
	double Ohm = pow(10, Ohm_dB/10.0);
	double h0 = sqrt(K * Ohm / (1.0 + K)); //direct
	double h1 = sqrt((Ohm / (1.0 + K)) / (1.0 + 1.0e-1 + 1.0e-2)); //indirect
	V htx_amp{h0, h1, h1 * 1.0e-1, h1 * 1.0e-2};
	int L = htx_amp.size(); //channel length
	VCOM htx =  Rician_fading(htx_amp);


	//	------------------					
	//
	//	Estimation phase
	//
	//  ------------------
	for(int i = 0; i < NN; i++){
		wave_tx[i] *= Amp_coe2;
		wave_aux[i] *= Amp_coe2;

		//IQ imbalance
		wave_tx[i] = eta_tx[0] * wave_tx[i] + eta_tx[1] * std::conj(wave_tx[i]);
		wave_aux[i] = eta_aux[0] * wave_aux[i] + eta_aux[1] * std::conj(wave_aux[i]);

		//Phase Noise
		wave_tx[i] = wave_tx[i] * PNR_tx[i];
		wave_aux[i] = wave_aux[i] * PNR_aux[i];

		//Amplifier
        wave_tx[i] = amp_tx.get_output_one(wave_tx[i]);

		//Wireless channel
		for(int j = 0; j < L; j++){
			if(i >= j) wave_rx[i] += wave_tx[i - j] * htx[j];
		}
		wave_rx[i] += wave_aux[i] * haux;

		//AWGN
		wave_rx[i] += Noise[i];
        
        // LNA is bypassed
        // wave_rx[i] = amp_rx.get_output_one(wave_rx[i]);

		//Phase Noise
		wave_rx[i] = wave_rx[i] * PNR_rx[i];

		//IQ imbalance
		wave_rx[i] = eta_rx[0] * wave_rx[i] + eta_rx[1] * std::conj(wave_rx[i]);
	}

	//	------------------					
	//
	//	Estimation
	//
	//  ------------------

	//WLフィルタ
	int Lw = 4;
	int Nw = ((P+1)*(P+3))/4; // # of basis functions
	VVCOM Fltr(Nw, VCOM(Lw)); // FIR filter
	LS_estimation(NN, Ns, Lw, P, wave_rx, wave_tx_cp, wave_aux_cp, Fltr, loop);

	//signal generation
	wave_tx.resize(NN);
	wave_aux.resize(NN);
	wave_rx.resize(NN);
	for(int i = 0; i < NN; i++){
		wave_tx[i] = 0;
		wave_aux[i] = 0;
		wave_rx[i] = 0;
	}
	VCOM wave(NN, 0);
	for(int i = 0; i < NN/(Ns + Ncp); i++){
		OFDM_qam(wave_part, Ns, Nact, 16);
		//Time domain
		fft_COMB(Ns,(double)2.0*std::acos(-1.0)/Ns, wave_part);
		for(int j = 0; j < Ns; j++){
			wave_part[j] /= (double)Ns;
		}
		for(int j = 0; j < Ns + Ncp; j++){
			wave[j + i * (Ns + Ncp)] = wave_part[j%Ns];
		}
	}
	a.CP_Lap_filtering(wave, Ns, NN, Ncp);

	for(int i = 0; i < NN; i++){
		wave[i] = wave[i] * fil_Amp_coe;
	}


	//Filtering
	for(int i = 0; i < NN; i++){
		wave_tx[i] = wave[i];
		wave_aux[i] = 0.0;
		for(int l = 1; l <= P; l += 2){
			for(int m = 0; m <= l; m++){
				for(int j = 0; j < Lw; j++){
					int k = ((l * l) + 3) / 4 + m - 1;
					if(i >= j) wave_aux[i] += Fltr[k][j] * mypow(wave[i - j], m) * mypow(std::conj(wave[i - j]), l-m);
				}
			}
		}
	}

	//AWGN generation
	AWGN_no(Floor, Noise);

	//	------------------					
	//
	//	Cancellation phase
	//
	//  ------------------
	VVCOM all_result(2, VCOM(NN));
	for(int i = 0; i < NN; i++){
		wave_tx[i] *= Amp_coe2;
		wave_aux[i] *= Amp_coe2;

		//IQ imbalance
		wave_tx[i] = eta_tx[0] * wave_tx[i] + eta_tx[1] * std::conj(wave_tx[i]);
		wave_aux[i] = eta_aux[0] * wave_aux[i] + eta_aux[1] * std::conj(wave_aux[i]);

		//Phase Noise
		wave_tx[i] = wave_tx[i] * PNR_tx[i];
		wave_aux[i] = wave_aux[i] * PNR_aux[i];

		//Amplifier
        wave_tx[i] = amp_tx.get_output_one(wave_tx[i]);

		//Wireless channel
		wave_rx[i] = 0;
		for(int j = 0; j < L; j++){
			if(i >= j) wave_rx[i] += wave_tx[i - j] * htx[j];
		}
		//Save result before cancellation
		all_result[0][i] = wave_rx[i];
		
		//SI cancellation point
		wave_rx[i] += wave_aux[i] * haux;

		//AWGN
		wave_rx[i] += Noise[i];

		//LNA
        wave_rx[i] = amp_rx.get_output_one(wave_rx[i]);

        //Phase Noise
		wave_rx[i] *= PNR_rx[i];

		//IQ imbalance
		wave_rx[i] = eta_rx[0] * wave_rx[i] + eta_rx[1] * std::conj(wave_rx[i]);

		//Save result after cancellation
		all_result[1][i] = wave_rx[i];
	}


	//Averaging
	VVCOM all_result_part(2, VCOM(Ns));
	VV result(2, V(Ns, 0));
	for(int i = 0; i < NN; i++){
		if(i%(Ns + Ncp) < Ns){
			all_result_part[0][i%(Ns + Ncp)] = all_result[0][i + Ncp];
			all_result_part[1][i%(Ns + Ncp)] = all_result[1][i + Ncp];
			if(i%(Ns + Ncp) == Ns - 1){
				fft_COMB(Ns,(double)-2.0*std::acos(-1.0)/Ns, all_result_part[0]);
				fft_COMB(Ns,(double)-2.0*std::acos(-1.0)/Ns, all_result_part[1]);
				for(int j = 0; j < Ns; j++){
					result[0][j] += pow(abs(all_result_part[0][j]),2);
					result[1][j] += pow(abs(all_result_part[1][j]),2);
				}
			}
		}
	}
	for(int i = 0; i < Ns; i++){
		result[0][i] = result[0][i]/(double)(NN/(Ns + Ncp));
		result[1][i] = result[1][i]/(double)(NN/(Ns + Ncp));
	}

	//NaN check
	double nan_ch = 0.0;
	for(int i = 0; i < Ns; i++){
		nan_ch += result[1][i];
	}
	if(isnan(nan_ch))
		std::cout << "NAN" << std::endl;

	return result;
}
