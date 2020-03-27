#pragma once

#include <vector>
#include <complex>
#include <iostream>
#include <cmath>

using COM = std::complex<double>;
using VCOM = std::vector<COM>;
using V = std::vector<double>;

class Amplifier
{
	private:
		V a;
		/*
		--------------------	"poly"	-------------------
			a[0] = 1st order coeff
			a[1] = 3rd order coeff
			a[2] = 5th order coeff
			a[3] = 7th order coeff

		--------------------	"rapp"	-------------------
			a[0] = small signal gain
			a[1] = saturated output
			a[2] = smoothness factor
		*/
	public:
		Amplifier(void);

		void set_coeff(V &coefficients);
		V get_coeff(void);

		void get_output(VCOM &x);
		COM get_output_one(COM x);

		// void plot_am_am(float start, float end);
};

Amplifier::Amplifier(void)
{
	a = V(3, 0.0);
}

void Amplifier::set_coeff(V &x)
{
	a = x;
}

V Amplifier::get_coeff(void)
{
	return a;
}

void Amplifier::get_output(VCOM &x)
{
	VCOM y(x.size(), 0.0);

	for (int i = 0; i < x.size(); ++i)
	{
		double magnitude = abs(x[i]);
		if (magnitude == 0.0)
			x[i] = 0.0;
		else
		{
			double theta = std::arg(x[i]);
			magnitude = a[0] * magnitude / pow( 1.0 + pow(a[0] * magnitude / a[1], 2.0 * a[2]), 0.5/a[2]);
			x[i] = COM(magnitude*cos(theta), magnitude*sin(theta));
		}
	}
}

COM Amplifier::get_output_one(COM x)
{
	COM y = 0.0;
	double magnitude = abs(x);
	if (magnitude == 0.0)
		y = 0.0;
	else
	{
		double theta = std::arg(x);
		magnitude = a[0] * magnitude / pow( 1.0 + pow(a[0] * magnitude / a[1], 2.0 * a[2]), 0.5/a[2]);
		y = COM(magnitude*cos(theta), magnitude*sin(theta));
	}
	return y;
}