#include <iostream>
#include <vector>
#include <complex>
#include <iterator>
#include <random>


template <typename A>
void AWGN_no(A N_dBm, std::vector<std::complex<A> > &data){
	std::random_device seed;
	std::default_random_engine generator(seed());

	A N = std::pow(10, ((N_dBm) / (A)20));
	std::normal_distribution<A> dist1(0.0, N/(A)sqrt(2.0));
	std::normal_distribution<A> dist2(0.0, N/(A)sqrt(2.0));

	for(int i = 0; i < data.size(); i++){
		data[i] = std::complex<A>(dist1(generator), dist2(generator));
	}

	// // power
	// double power = 0.0;
	// for(int i = 0; i < data.size(); i++){
	// 	power += pow(abs(data[i]),2);
	// }
	// std::cout << "Power" << 10.0 * log10(power/(double)data.size()) << " dB" << std::endl;
}

