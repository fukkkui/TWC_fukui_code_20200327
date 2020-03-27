#include <iostream>
#include <vector>
#include <complex>

using V = std::vector<double>;
using COM = std::complex<double>;
using VCOM = std::vector<COM>;
using VVCOM = std::vector<VCOM>;

extern "C"
int LAPACKE_zgelss(
        int matrix_order, int m, int n,
        int nrhs, COM* a,
        int lda, COM* b,
        int ldb, double* s, double rcond,
        int* rank);

COM mypow(COM x, int power){
	COM result;
	if(power == 0){
		result = COM(1.0, 0.0);
	}else{
		result = pow(x, power);
	}
	return result;
}

void LS_estimation(int NN, int Ns, int L, int P,
							VCOM &y,
							VCOM &data_tx_cp,
							VCOM &data_aux_cp,
							VVCOM &W,
							int loop){

	int Nw = ((P + 1)*(P + 3))/4;
	VVCOM x(NN, VCOM(L * Nw + 2));
	for(int i = 0; i < NN; i++){
		for(int j = 0; j < L; j++){
			for(int l = 1; l <= P; l += 2){
				for(int m = 0; m <= l; m++){
					if(i >= j){
						x[i][j + (((l * l) + 3)/4 + (m - 1)) * L] = mypow(data_tx_cp[i - j], m) * mypow(std::conj(data_tx_cp[i - j]), l-m);
					}else{
						x[i][j + (((l * l) + 3)/4 + (m - 1)) * L] = 0;
					}
				}
			}
		}
		x[i][L * Nw] = data_aux_cp[i];
		x[i][L * Nw + 1] = std::conj(data_aux_cp[i]);
	}

	VCOM x_line(NN * (L * Nw + 2));
	for(int j = 0; j < L * Nw + 2; j++){
		for(int i = 0; i < NN; i++){
			x_line[i + j * NN] = x[i][j];
		}
	}


	//LSM
	V workspace(std::min(NN, L * Nw + 2));
    int rankN;
    LAPACKE_zgelss(102, NN, L * Nw + 2, 1,
        &x_line[0], NN,
        &y[0],
        (int)std::max(NN, L * Nw + 2),
        &workspace[0],
        1e-8,
        &rankN);

	VVCOM c(Nw, VCOM(L));
	for(int i = 0; i < Nw; i++){
		for(int j = 0; j < L; j++){
			c[i][j] = y[j + L * i];
		}
	}
	COM Caux1 = y[L * Nw];
	COM Caux2 = y[L * Nw + 1];

	//filter
	W.resize(Nw, VCOM(L));
	for(int l = 1; l <= P; l += 2){
		for(int m = 0; m <= l; m++){
			for(int i = 0; i < L; i++){
				int j = ((l * l) + 3)/4 + m - 1;
				int k = ((l * l) + 3)/4 + l - m - 1;
				COM denomi = pow(abs(Caux2),2) - pow(abs(Caux1),2);
				W[j][i] = (std::conj(Caux1) * c[j][i] - Caux2 * std::conj(c[k][i]))/denomi;
			}
		}
	}

}

