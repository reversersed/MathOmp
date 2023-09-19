#include <fstream>
#include <iostream>
#include <cstring>
#include <cmath>
#include <omp.h>

const int n = 2;
const double t0 = 0;
const double t_max = 10;
const double tau = 0.000001;
double y[n] = { 0.0, 0.05 };
double yy[n];

double f(double* y, double time, int i)
{
	switch (i)
	{
	case 0: return y[1];
	case 1: return (-0.19 * time * y[1] - pow(time,2.0) * y[0]);
	default: return 0.0;
	}
}

int main()
{
	setlocale(LC_ALL, "Russian");
	double time_begin, time_end;
	
	//Euler
	std::cout << "Явный метод Эйлера:" << std::endl;
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
	std::cout << "Время расчетов: " << (time_end - time_begin) << " секунд"  << std::endl;
	for (int i = 0; i < n; ++i)
	{
		std::cout << "y[" << i+1 << "] = " << y[i] << std::endl;
	}

	//Runge-Kutt 2
	std::cout << std::endl << "Метод Рунге-Кутта 2:" << std::endl;

	y[0] = 0.0, y[1] = 0.05;
	yy[0] = 0, yy[1] = 0;
	double ff[n];

	time_begin = omp_get_wtime();
	for (double x = t0; x <= t_max; x += tau)
	{
		for (int i = 0; i < n; i++)
			yy[i] = y[i] + tau * 0.5 * f(y, x, i);
		for (int i = 0; i < n; i++)
			ff[i] = y[i] + tau * f(yy, x + tau * 0.5, i);
		for (int i = 0; i < n; i++)
			y[i] = ff[i];
	}
	time_end = omp_get_wtime();
	std::cout << "Время расчетов: " << (time_end - time_begin) << " секунд" << std::endl;

	for (int i = 0; i < n; i++)
		std::cout << "y[" << i+1 << "] = " << y[i] << "\n";

	//Predictor-Corrector
	std::cout << std::endl << "Метод Предиктор-Корректор:" << std::endl;

	y[0] = 0.0, y[1] = 0.05;
	yy[0] = 0, yy[1] = 0;
	ff[0] = 0, ff[1] = 0;

	time_begin = omp_get_wtime();
	for (double x = t0; x <= t_max; x += tau)
	{
		for (int i = 0; i < n; i++)
			yy[i] = y[i] + tau * f(y, x, i);
		for (int i = 0; i < n; i++)
			ff[i] = y[i] + tau * (f(y, x, i) + f(yy, x + tau, i)) / 2;
		for (int i = 0; i < n; i++)
			y[i] = ff[i];
	}
	time_end = omp_get_wtime();
	std::cout << "Время расчетов: " << (time_end - time_begin) << " секунд" << std::endl;

	for (int i = 0; i < n; i++)
		std::cout << "y[" << i+1 << "] = " << y[i] << "\n";

	//Runge-Kutt 4
	std::cout << std::endl << "Метод Рунге-Кутта 4:" << std::endl;

	y[0] = 0.0, y[1] = 0.05;
	yy[0] = 0, yy[1] = 0;
	double r1[n] = { 0.0 }, 
		r2[n] = { 0.0 }, 
		r3[n] = { 0.0 }, 
		r4[n] = { 0.0 }, 
		buf[n] = { 0.0 };


	time_begin = omp_get_wtime();
	for (double x = t0; x <= t_max; x += tau)
	{
		for (int i = 0; i < n; i++)
			r1[i] = tau * f(y, x, i);
		for (int i = 0; i < n; i++)
			buf[i] = y[i] + r1[i] * 0.5;
		for (int i = 0; i < n; i++)
			r2[i] = tau * f(buf, x + tau * 0.5, i);
		for (int i = 0; i < n; i++)
			buf[i] = y[i] + r2[i] * 0.5;
		for (int i = 0; i < n; i++)
			r3[i] = tau * f(buf, x + tau * 0.5, i);
		for (int i = 0; i < n; i++)
			buf[i] = y[i] + r3[i];
		for (int i = 0; i < n; i++)
			r4[i] = tau * f(buf, x + tau, i);
		for (int i = 0; i < n; i++)
			yy[i] = y[i] + (r1[i] + r2[i] * 2 + r3[i] * 2 + r4[i]) / 6;
		for (int i = 0; i < n; i++)
			y[i] = yy[i];
	}
	time_end = omp_get_wtime();
	std::cout << "Время расчетов: " << (time_end - time_begin) << " секунд" << std::endl;

	for (int i = 0; i < n; i++)
		std::cout << "y[" << i + 1 << "] = " << y[i] << "\n";

	//Euler
	std::cout << std::endl << "Неявный метод Эйлера:" << std::endl;

	y[0] = 0.0, y[1] = 0.05;
	yy[0] = 0, yy[1] = 0;
	double a[2][2] = { {-1 / 0.001,1},{0,-1 / 0.001} },
			fs[n] = { 0.0 };

	time_begin = omp_get_wtime();
	for (double x = t0; x <= t_max; x += tau)
	{
		double d[3] = {};
		for (int i = 0; i < n; i++)
			fs[i] = -f(y, x, i);

		d[0] = a[0][0] * a[1][1] - a[0][1] * a[1][0];
		d[1] = a[1][1] * fs[0] - fs[1] * a[0][1];
		d[2] = a[0][0] * fs[1] - fs[0] * a[1][0];
		for (int i = 0; i < n; i++)
			yy[i] = d[i + 1] / d[0];
		for (int i = 0; i < n; i++)
			y[i] = y[i] + yy[i];
	}
	time_end = omp_get_wtime();
	std::cout << "Время расчетов: " << (time_end - time_begin) << " секунд" << std::endl;

	for (int i = 0; i < n; i++)
		std::cout << "y[" << i + 1 << "] = " << y[i] << "\n";
	return 0;
}