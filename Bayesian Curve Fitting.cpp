//Predictive_function
#include "pch.h"
#include <iostream>
#include <armadillo>
#include <string>
#include <fstream>

using namespace std;
using namespace arma;
double A = 0.1;		// value of alpha
double B =12;		//known value of Beta
const int N=20;		// known value of N, leave one value for testing 21 value
const int M = 4;				//5 or 10 times of dataset size, degree
void initialize_t(mat&, double[],string);
mat initialize_phi(mat X, int x, int M);
mat initialize_phi_t(mat,int,mat,int,int);			//where x =1,2,3,4,5......N
void initialize_X_data(mat&);
mat evaluate_s_inverse(mat, mat,mat,mat);			//beta, alpha,X_data and phi of x
mat sum_phi_T(mat, mat);
double read_real_value(string);

int main()
{
	//*********************BETA AND ALPHA****************************
	mat Beta = randu<mat>(1,1);
	Beta = 0.0;
	mat alpha = randu<mat>(1, 1);
	alpha = 0.0;
	alpha(0, 0) = 0.1;
	Beta(0, 0) = B;
	cout <<"BETA:"<< Beta;
	cout << "ALPHA:"<<alpha << endl;
	//********************initialize data_t*****************
	double t_data_copy[N];
	string symbol_name;
	cout << "Enter the symbol of the company:";
	cin >> symbol_name;
	string filename = symbol_name + ".txt";
	mat t_data = randu<mat>(N, 1);
	initialize_t(t_data,t_data_copy,filename);
	//********************initialize X**********************
	mat X_data = randu<mat>(N, 1);
	initialize_X_data(X_data);
	//cout << X_data << endl;
	//********************find sum of phi and t*************
	mat summation_phi = randu<mat>(M, 1);
	summation_phi = 0.0;
	summation_phi = sum_phi_T(X_data,t_data);
	//cout <<summation_phi << endl;
	//****************FIND X_PREDICT***************
	double x_predict = 21;
	//***************phi and phi_t********************
	mat phi_X_predict = randu<mat>(M, 1);
	phi_X_predict = initialize_phi(X_data, x_predict, M);
	mat phi_X_transpose = phi_X_predict.t();
	//****************find_inverse****************
	mat S = randu<mat>(M, M);
	S = evaluate_s_inverse(Beta,alpha,X_data, phi_X_transpose);
	//cout << S << endl;
	//*************CALCULATE MEAN****************************
	mat mean = randu<mat>(1, 1);
	mean = B*phi_X_transpose*(S*summation_phi);				//which is also the mode, therefore maximum probable value
	double predicted_value = as_scalar(mean);				//convert to scalar
	cout << "Mean and Mode of the probability distribution function:" << predicted_value << endl;
	//**************CALCULATE VARIANCE************************
	mat variance = randu<mat>(1, 1);
	variance = (1 / B) + phi_X_transpose * S*phi_X_predict;
	double variance_value = as_scalar(variance);
	cout << "Variance of the probability distribution function:" << variance_value << endl;
	//**************CALCULATE ERROR***************************
	double actual_value = read_real_value(symbol_name);
	cout << "The actual value for the company from historical data is:" << actual_value << endl;
	double absolute_error = predicted_value - actual_value;
	cout << "The absolute error is:" << absolute_error << endl;;
	double relative_error = 100 * (predicted_value - actual_value) / actual_value;
	cout << "The relative error is:" << relative_error << "%" << endl;
	return 0;
}

void initialize_t(mat& t, double t_data[],string filename)
{
	int i = 0;
	ifstream read;
	read.open(filename);
	string item;
	while (!read.eof())
	{
		read >> item;
		t_data[i] = stod(item);	
		i = i + 1;
	}
	read.close();
	for (int j = 0; j < N;j++)
	{
		t(j, 0) = t_data[j];
	}
}
void initialize_X_data(mat& x_data)
{
	for (int i = 0; i < N; i++)
	{
		x_data(i, 0) = i;
	}
}
mat initialize_phi(mat X, int x,int M)     //initalize phir and get tthe value for phi*t for specific x and t
{
	int k;
	mat phi_x = randu<mat>(M, 1);
	if (x < N)
	{
		k = X(x, 0);
	}
	else
		k = x;
	for (int i = 0; i < M; i++)
	{
		phi_x(i, 0) = pow(k, i);
	}
	return phi_x;
}
mat initialize_phi_t(mat X,int x,mat T,int t,int M)     //initalize phir and get tthe value for phi*t for specific x and t
{
	mat phi_x = randu<mat>(M, 1);
	int k = X(x,0);
	int k2 = T(t, 0);
	for (int i = 0; i < M; i++)
	{
		phi_x(i,0) = pow(k, i)*k2;
	}
	return phi_x;
}
mat sum_phi_T(mat X, mat t)
{
	mat phi_t_result = randu<mat>(M, 1);
	for (int i = 0; i < N; i++)
	{
		mat phi_multiplied = initialize_phi_t(X,i,t,i,M);
		phi_t_result = phi_t_result + phi_multiplied;
	}
	return phi_t_result;
}
mat evaluate_s_inverse(mat b, mat a,mat X_data, mat phi_x)				//evaluate s inverse but send inverse of that because we need S
{
	mat identity = eye<mat>(M,M);
	mat sum_phi_phiT = randu<mat>(M, M);
	for (int i = 0; i < N; i++)
	{
		sum_phi_phiT += initialize_phi(X_data, i, M)*(phi_x);
	}
	mat s_inverse = randu<mat>(M, M);
	s_inverse = (A * identity) + (B * sum_phi_phiT);
	return s_inverse.i();
}
double read_real_value(string symbol)
{
	ifstream read;
	string item;
	double value;
	read.open("ActualValues_21st_day.txt");
	while (!read.eof())
	{
		read >> item >> value;
		if (item == symbol)
			break;
	}
	read.close();
	return value;
}
































