#include <fstream>
#include <iostream>
#include <cstring>
#include <cmath>
#include <omp.h>

const int n = 2;
const double t0 = 0;
const double t_max = 10;
const double tau = 0.0000001;
double y[n] = { 0.0, 0.05 };
double yy[n];

double f(double* y, double time, int i)
{
	switch (i)
	{
	case 0: return y[1];
	case 1: return (-0.19 * time * y[1] - pow(time,2.0) * y[0]);
	}
}

int main()
{
	setlocale(LC_ALL, "Russian");
	double time_begin, time_end, time_elapsed;

	time_begin = omp_get_wtime();
	double t = t0-tau;
	while ((t = t + tau) <= t_max)
	{
		for (int i = 0; i < n; i++)
			yy[i] = y[i] + tau * f(y, t, i);
		for (int i = 0; i < n; i++)
			y[i] = yy[i];
	}
	time_end = omp_get_wtime();
	std::cout << "Время расчетов: " << std::fixed << (time_end - time_begin) << " секунд"  << std::endl;
	for (int i = 0; i < n; ++i)
	{
		std::cout << std::fixed << "y[" << i+1 << "]: " << y[i] << std::endl;
	}

	
	return 0;
}

