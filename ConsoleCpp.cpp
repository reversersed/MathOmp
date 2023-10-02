#include <fstream>
#include <iostream>
#include <cstring>
#include <cmath>
#include <omp.h>

void ExplicitEuler(double* y, int size, double (*function[])(double*, double), double t0, double t_max, double tau)
{
	double t = t0 - tau;
	while ((t = t + tau) <= t_max)
	{
		for (int i = 0; i < size; i++)
			y[i] += tau * function[i](y, t);
	}
}
void RungeKutt2(double* y, int size, double (*function[])(double*, double), double t0, double t_max, double tau)
{
	double* temporary_y = new double[size];

	for (double x = t0; x <= t_max; x += tau)
	{
		for (int i = 0; i < size; i++)
			temporary_y[i] = y[i] + tau * 0.5 * function[i](y, x);
		for (int i = 0; i < size; i++)
			y[i] += tau * function[i](temporary_y, x + tau * 0.5);
	}

	delete[] temporary_y;
}
void PredictorCorrector(double* y, int size, double (*function[])(double*, double), double t0, double t_max, double tau)
{
	double* temporary_y = new double[size];
	double* function_value = new double[size];

	for (double x = t0; x <= t_max; x += tau)
	{
		for (int i = 0; i < size; i++)
		{
			function_value[i] = function[i](y, x);
			temporary_y[i] = y[i] + tau * function_value[i];
		}
		for (int i = 0; i < size; i++)
			y[i] += tau * (function_value[i] + function[i](temporary_y, x + tau)) / 2.0;
	}

	delete[] temporary_y;
	delete[] function_value;
}
void RungeKutt4(double* y, int size, double (*function[])(double*, double), double t0, double t_max, double tau)
{
	double 
		*temporary_r1 = new double[size],
		*temporary_r2 = new double[size],
		*temporary_r3 = new double[size],
		*temporary_r4 = new double[size],
		*buffer = new double[size];

	for (double x = t0; x <= t_max; x += tau)
	{
		for (int i = 0; i < size; i++)
			temporary_r1[i] = tau * function[i](y, x);
		for (int i = 0; i < size; i++)
			buffer[i] = y[i] + temporary_r1[i] * 0.5;
		for (int i = 0; i < size; i++)
			temporary_r2[i] = tau * function[i](buffer, x + tau * 0.5);
		for (int i = 0; i < size; i++)
			buffer[i] = y[i] + temporary_r2[i] * 0.5;
		for (int i = 0; i < size; i++)
			temporary_r3[i] = tau * function[i](buffer, x + tau * 0.5);
		for (int i = 0; i < size; i++)
			buffer[i] = y[i] + temporary_r3[i];
		for (int i = 0; i < size; i++)
			temporary_r4[i] = tau * function[i](buffer, x + tau);
		for (int i = 0; i < size; i++)
			y[i] += (temporary_r1[i] + temporary_r2[i] * 2.0 + temporary_r3[i] * 2.0 + temporary_r4[i]) / 6.0;
	}

	delete[] temporary_r1;
	delete[] temporary_r2;
	delete[] temporary_r3;
	delete[] temporary_r4;
	delete[] buffer;
}
void ImplicitEuler(double* y, int size, double (*function[])(double*, double), double t0, double t_max, double tau)
{
	double *temporary_f = new double[size];
	double a[2][2] = { {-1.0 / tau,1},{0,-1.0 / tau} };
	double d[3];

	for (double x = t0; x <= t_max; x += tau)
	{
		a[1][0] = -x*x;
		a[1][1] = -0.19*x-1.0/tau;

		for (int i = 0; i < size; i++)
			temporary_f[i] = -function[i](y, x);

		d[0] = a[0][0] * a[1][1] - a[0][1] * a[1][0];
		d[1] = a[1][1] * temporary_f[0] - temporary_f[1] * a[0][1];
		d[2] = a[0][0] * temporary_f[1] - temporary_f[0] * a[1][0];

		if (d[0] == 0)
		{
			std::cout << "Неудается расчитать значения. Во время расчета определителей происходит деление на 0" << std::endl;
			return;
		}
		for (int i = 0; i < size; i++)
			y[i] += d[i + 1] / d[0];
	}
	delete[] temporary_f;
}

const int n = 2;
const double tau = 0.000001;
const double t0 = 0;
const double t_max = 10;
double y[n] = { 0.0, 0.05 };


double f1(double* y, double time) {
	return y[1];
}
double f2(double* y, double time) {
	return (-0.19 * time * y[1] - pow(time, 2.0) * y[0]);
}

int main()
{
	setlocale(LC_ALL, "Russian");
	double time_begin, time_end;
	double (*function[])(double*, double) = { f1, f2 };
	
	//Euler
	std::cout << "Явный метод Эйлера:" << std::endl;
	time_begin = omp_get_wtime();
	ExplicitEuler(y, n, function, t0, t_max, tau);
	time_end = omp_get_wtime();

	std::cout << "Время расчетов: " << (time_end - time_begin) << " секунд"  << std::endl;
	for (int i = 0; i < n; ++i)
		std::cout << "y[" << i+1 << "] = " << y[i] << std::endl;
	
	//Runge-Kutt 2
	std::cout << std::endl << "Метод Рунге-Кутта 2:" << std::endl;

	y[0] = 0.0, y[1] = 0.05;
	time_begin = omp_get_wtime();
	RungeKutt2(y, n, function, t0, t_max, tau);
	time_end = omp_get_wtime();

	std::cout << "Время расчетов: " << (time_end - time_begin) << " секунд"  << std::endl;
	for (int i = 0; i < n; ++i)
		std::cout << "y[" << i+1 << "] = " << y[i] << std::endl;

	//Predictor-Corrector
	std::cout << std::endl << "Метод Предиктор-Корректор:" << std::endl;

	y[0] = 0.0, y[1] = 0.05;
	time_begin = omp_get_wtime();
	PredictorCorrector(y, n, function, t0, t_max, tau);
	time_end = omp_get_wtime();

	std::cout << "Время расчетов: " << (time_end - time_begin) << " секунд"  << std::endl;
	for (int i = 0; i < n; ++i)
		std::cout << "y[" << i+1 << "] = " << y[i] << std::endl;

	//Runge-Kutt 4
	std::cout << std::endl << "Метод Рунге-Кутта 4:" << std::endl;

	y[0] = 0.0, y[1] = 0.05;
	time_begin = omp_get_wtime();
	RungeKutt4(y, n, function, t0, t_max, tau);
	time_end = omp_get_wtime();

	std::cout << "Время расчетов: " << (time_end - time_begin) << " секунд"  << std::endl;
	for (int i = 0; i < n; ++i)
		std::cout << "y[" << i+1 << "] = " << y[i] << std::endl;

	//Euler
	std::cout << std::endl << "Неявный метод Эйлера:" << std::endl;

	y[0] = 0.0, y[1] = 0.05;
	time_begin = omp_get_wtime();
	ImplicitEuler(y, n, function, t0, t_max, tau);
	time_end = omp_get_wtime();

	std::cout << "Время расчетов: " << (time_end - time_begin) << " секунд"  << std::endl;
	for (int i = 0; i < n; ++i)
		std::cout << "y[" << i+1 << "] = " << y[i] << std::endl;
	return 0;
}