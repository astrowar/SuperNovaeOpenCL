// Rotina de integracao multi dimensional
#ifndef INTEGRATEFUNCTIONS_HPP
#define INTEGRATEFUNCTIONS_HPP
#include <vector>
#include <cassert>
using namespace  std;

//====================================================================================
 
double Integra_1n(double *y, vector<double> &x); 
double Integra_1n(double *y, int i);
  

//====================================================================================
template <  typename... Rest>
double Integra_2n(double **y, vector<double> &x, Rest... rest)
{
	double soma = 0.0;
	const int ifinal = x.size();
	//if (ifinal == 1) return  Integra_1n(y[0], rest...);
	assert(ifinal > 1);
	if (ifinal > 3)
	{
		double ya = Integra_1n(y[0], rest...);
		double yb = 0;
		double yc = 0;
		for (int i = 0; i < ifinal - 3; i += 2)
		{
			yb = Integra_1n(y[i + 1], rest...);
			yc = Integra_1n(y[i + 2], rest...);
			soma += (ya +4.0* yb + yc)*(x[i + 2] - x[i]);
			ya = yc;
		}
		//assert(soma > DBL_EPSILON);
		return soma / 6.0;
	}
	for (int i = 0; i < ifinal - 1; ++i)
	{
		soma += Integra_1n(y[i], rest...)*(x[i + 1] - x[i]);
	}
	//assert(soma > DBL_EPSILON);
	return soma;
}

template <  typename... Rest>
double Integra_2n(double **y, int i, Rest... rest)
{
	return Integra_1n(y[i], rest...);
}

//====================================================================================
template <  typename... Rest>
double Integra_3n(double ***y, vector<double> &x, Rest... rest)
{
	double soma = 0.0;
	const int ifinal = x.size();
	//if (ifinal == 1) return Integra_2n(y[0], rest...);
	assert(ifinal > 1);
	if (ifinal > 3)
	{
		double ya = Integra_2n(y[0], rest...);
		double yb = 0;
		double yc = 0;
		for (int i = 0; i < ifinal - 3; i += 2)
		{
			yb = Integra_2n(y[i + 1], rest...);
			yc = Integra_2n(y[i + 2], rest...);
			soma += (ya +4.0* yb + yc)*(x[i + 2] - x[i]);
			ya = yc;
		}
		//assert(soma > DBL_EPSILON);
		return soma / 6.0;
	}
	for (int i = 0; i < ifinal - 1; ++i)
	{
		soma += Integra_2n(y[i], rest...)*(x[i + 1] - x[i]);
	}
	//assert(soma > DBL_EPSILON);
	return soma;
}

template <  typename... Rest>
double Integra_3n(double ***y, int i, Rest... rest)
{
	return Integra_2n(y[i], rest...);
}


//====================================================================================
template <  typename... Rest>
double Integra_4n(double **** y, vector<double> &x, Rest... rest)
{
	double soma = 0.0;
	const int ifinal = x.size();
	//if (ifinal == 1) return Integra_3n(y[0], rest...);
	assert(ifinal > 1);
	if (ifinal > 3)
	{
		double ya = Integra_3n(y[0], rest...);
		double yb = 0;
		double yc = 0;
		for (int i = 0; i < ifinal - 3; i += 2)
		{
			yb = Integra_3n(y[i + 1], rest...);
			yc = Integra_3n(y[i + 2], rest...);
			soma += (ya +4.0* yb + yc)*(x[i + 2] - x[i]);
			ya = yc;
		}
		//assert(soma > DBL_EPSILON);
		return soma / 6.0;
	}
	for (int i = 0; i < ifinal - 1; ++i)
	{
		soma += Integra_3n(y[i], rest...)*(x[i + 1] - x[i]);
	}
	//assert(soma > DBL_EPSILON);
	return soma;
}

template <  typename... Rest>
double Integra_4n(double **** y, int i, Rest... rest)
{
	return Integra_3n(y[i], rest...);
}



//====================================================================================
template <  typename... Rest>
double Integra_5n(double ***** y, vector<double> &x, Rest... rest)
{
	double soma = 0.0;
	const int ifinal = x.size();
	//if (ifinal == 1) return Integra_4n(y[0], rest...);
	assert(ifinal > 1);
	if (ifinal > 3)
	{
		double ya = Integra_4n(y[0], rest...);
		double yb = 0;
		double yc = 0;
		for (int i = 0; i < ifinal - 3; i += 2)
		{
			yb = Integra_4n(y[i + 1], rest...);
			yc = Integra_4n(y[i + 2], rest...);
			soma += (ya +4.0* yb + yc)*(x[i + 2] - x[i]);
			ya = yc;
		}
		//assert(soma > DBL_EPSILON);
		return soma / 6.0;
	}
	for (int i = 0; i < ifinal - 1; ++i)
	{
		soma += Integra_4n(y[i], rest...)*(x[i + 1] - x[i]);
	}
	//assert(soma > DBL_EPSILON);
	return soma;
}

template <  typename... Rest>
double Integra_5n(double ***** y, int i, Rest... rest)
{
	return Integra_4n(y[i], rest...);
}


//====================================================================================

template <  typename... Rest>
double Integra_6n(double ****** y, vector<double> &x, Rest... rest)
{
	double soma = 0.0;
	const int ifinal = x.size();
	//if (ifinal == 1) return Integra_5n(y[0], rest...);
	assert(ifinal > 1);
	 
	if (ifinal > 3)
	{
		double ya = Integra_5n(y[0], rest...);
		double yb = 0;
		double yc = 0;
		for (int i = 0; i < ifinal - 3; i += 2)
		{
		 
			yb = Integra_5n(y[i + 1], rest...);
			yc = Integra_5n(y[i + 2], rest...);
			soma += (ya + 4.0 * yb + yc)*(x[i + 2] - x[i]);
			ya = yc;			
		}
		//assert(soma > DBL_EPSILON);
		return soma / 6.0;
	}
  
		for (int i = 0; i < ifinal - 1; ++i )
		{ 
			soma += Integra_5n(y[i ], rest...)*(x[i + 1] - x[i]); 
		}
		//assert(soma > DBL_EPSILON);
		return soma;
	 
	 
}
 
template <  typename... Rest>
double Integra_6n(double ****** y, int i, Rest... rest)
{
	return Integra_5n(y[i], rest...); 
}


#endif // INTEGRATEFUNCTIONS_HPP
