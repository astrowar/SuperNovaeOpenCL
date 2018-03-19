#include "IntegrateFunctions.hpp"
#include <vector>
#include <cassert>
using namespace  std;

double Integra_1n(double *y, vector<double> &x)
{
	double soma = 0.0;
	const int ifinal = x.size();
	//if (ifinal == 1) return y[0];
	assert(ifinal > 1); 
	if (ifinal > 3)
	{
 
		for (int i = 0; i < ifinal - 3; i += 2)
		{
			double ya = y[i ];
			double yb = y[i + 1];
			double yc = y[i + 2];
			soma += (ya + 4.0 * yb + yc)*(x[i + 2] - x[i]); 
		}
		//assert(soma > DBL_EPSILON);
		return soma / 6.0;
	}
	for (int i = 0; i < ifinal -1  ; ++i)
	{
		soma += y[i] * (x[ i+1] - x[i]);
	}
	//assert(soma > DBL_EPSILON);
	return soma;
}


double Integra_1n(double *y, int i)
{
	return  y[i];
}