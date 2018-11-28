// snova.cpp : Defines the entry point for the console application.
//

 

#define _CRT_SECURE_NO_WARNINGS


#include <cstdio>
#include <cstdlib>
#include <cmath>


#include <vector>
#include <ctime>
#include <algorithm>
#include <time.h>


#include "BidimensionalMap.hpp"

#include  "IntegrateFunctions.hpp"
#include "generate_matrix.hpp"

#define NP 400

#define Q 1.29
#define me 0.511
#define eta 0.0
#define Mk  2.14     //  Kamiokande detector mass (kton).
#define Mimb 6.8     // IMB detector mass (kton).
#define Mbaksan 0.28 // Baksan detector mass (kton).
#define fk 1.0
#define fimb 0.9055
#define fb 1.0
#define Cn    (1.0/(4 * 3.1416))
#define lnVk   14.56
#define lnVimb  15.73
#define delt 1.0//2.3
#define timeK 20.0//10.43
#define timeIMB 16.0 //5.9



 


#define real double


using namespace std;

#define Integra( r, ff , x, x1 , x2 , dx )  for((x) = (x1) ; (x) <= (x2) ; (x)+=0  ){ (r) = (r) + (dx) * ( ff )/6.0 ; (x)+=(dx)/2.0 ; (r) = (r) + 4*(dx) * ( ff )/6.0; (x)+=(dx)/2.0 ;(r) = (r) + (dx) * ( ff )/6.0; };


#define DUMP(  ff , x, x1 , x2 , dx )  for((x) = (x1) ; (x) <= (x2) ;(x)+=(dx)   ){ printf(">%f %f \n",x, ( ff ) );};

#define H  6.0

//regiao de calculo da massa dos distibuicoes T
#define Tmin 2.0
#define Tmax 8.0
#define ddT  0.02

#define sigma_tp 1.0

//  tempo
#define tpmin 0.01 
#define tpmax 0.1
#define ddtp  1.03 


// alpha
#define Amin 0.2
#define Amax 20.0
#define ddA  0.1

//  eps
#define epsmin 1.0
#define epsmax 70.0  
#define ddeps  5.0  


// tau1

#define  tau1min 0.1
#define  tau1max 15.0   
#define  ddtau1  0.3  


// tau2

#define  tau2min 0.0
#define  tau2max 4.0
#define  ddtau2  5.0


//Scale parameter
#define apmin  0.01
#define apmax  0.4
#define ddap   0.4


// Noise's parameters

#define pk     1.0
#define pimb   1.0
#define pb     1.0
#define effK   1e-5
#define effIMB 1e-5


//double PARMS[5];
//(FILE*)FILES[5];
//int NUMPARAM=5;





double Meff;


//vector <double> vAlpha, vAlpha2 ,vTemp , vDeltaTp , vAp ,vG  ;  //vetores que armazenam os parametros

vector <double> vAlpha, vAlpha2, vTemp, vtp, vAp, vtau1, vtau2;  //vetores que armazenam os parametros



double *****L;

double ******Lk; //vetor que armazena a likehood  K base para todos

double ****Limb; // vetor que armazena a likelihood IMB base para todos

double ******L_s12; //vetor que armazena a likehood S1 e S2

double *****L_g; //vetor que armazena a likehood g

double ***L_a; //vetor que armazena a likehood a

double Lmax;
//********************************************************** Data ***********************************************************************




 


typedef struct LikelihoodParameter
{
	real alpha, T, ap, tp, tau1, tau2;
	real result;

	LikelihoodParameter(real _alpha, real _T, real _ap, real _tp, real _tau1, real _tau2);
} LikelihoodParameter;
LikelihoodParameter::LikelihoodParameter(real _alpha, real _T, real _ap, real _tp, real _tau1, real _tau2)
{

	alpha = _alpha;
	T = _T;
	ap = _ap;
	tp = _tp;
	tau1 = _tau1;
	tau2 = _tau2;
	result = 0.0;
}




//**************************************************** Priori functions  ****************************************************************


double priori1(double alpha, double T, double tp, double ap, double tau1, double tau2) {


	
	// return 1.0 / (  pow(T, 2.0));

	//	return  1.0/ pow(alpha,2)  ;
	const double pAlpha = 1.0 / pow(alpha, 1.0);
	const double pTemp =  1.0 / pow(T, 2.0) ;
	//double pTau = 1.0 / (pow(tau1, 0.2) * pow(tau2, 0.2));


	return pTemp  * pAlpha;


}



//**************************************************** Likelihood's Calculus *************************************************************


//**************************************************  Predictive's calculus  **********************************************************


// Prior 1 --> two gaussians:


double predictive1(void) {


	double soma = 1.0;

	return soma;

}

// Prior 2 --> 1 / M:

double predictive2(void) {

	return 0.0;
}

//********************************************** Probabilities calculus **************************************************************

//Probability 1:

double Probability1(int i1, int i2, int i3, int i4, int i5, int i6, int i7) {
	 
	return 1;

}

void computeParams(std::vector<LikelihoodParameter> &params);




float Linter(float x, int ki ,std::vector<float>& xp , int k1, int k2)
{
	float s = 1.0;
	for (int km = k1; km <= k2; ++km)
	{
		if (ki != km) { s = s*(x - xp[km]) / (xp[ki] - xp[km]); }
	}
	return s;
}

float interpolation(float x, std::vector<float>& xp, std::vector<float>& yp)
{
	int i;
	int n = xp.size();
	if (x <= xp[0]) return yp[0];
	if (x >= xp[n-1]) return yp[n-1];

	for (i = 0; i < n-1; ++i)
	{
		if (x >= xp[i] && x <= xp[i + 1])
		{
			break;
		}
	}

	int k1 = std::max(0, i - 2);
	int k2 = std::min(n-1, i + 2);

	float y = 0;
	for (int k = k1; k <= k2; ++k)
	{
		y += Linter( x, k , xp, k1,k2) *yp[k];
	}
	return y;
}











std::pair<float,float> max_value(std::vector<float>& xp, std::vector<float>& yp)
{
	float ymax = yp[0];
	int jmax = 0;
	//encontra o ponto mais alto
	for (size_t j = 0; j < yp.size(); ++j)
	{
		if (yp[j] > ymax) { ymax = yp[j]; jmax = j; }
	}
	int j1 = std::max(0, jmax - 2);
	int j2 = std::min( jmax + 2, int(yp.size()-1));
	//procura entorno desse valor o maximo com precisao
	float x1 = xp[j1];
	float x2 = xp[j2];
	float dx = (x2 - x1) / 20.0f;

	ymax = yp[0];
	float xmax = xp[0];

	for (int loop = 0; loop < 3; ++loop)
	{
		for (float x = x1; x <= x2; x += dx)
		{
			float y = interpolation(x, xp, yp);
			if (y > ymax) { ymax = y; xmax = x; }
		}
		x1 = xmax - 2 * dx;
		x2 = xmax + 2 * dx;
		dx = (x2 - x1) / 20.0f;
	}
	return std::pair<float, float>(xmax,ymax);
}

int build_Program();












int main() {
	double i;
	double x;
	double y;
	double z;

	double x1, x2, x3, x4, x5;
	double lmax1, lmax2, lmax3, lzmax;
	double predi;
	double ei1, ei2, ei3, ei4, ei5, ei6, ei7;

	double tmp_max = 0;

	printf("Start ! \n");

	Lmax = 0.0;


	ei1 = 0;
	ei2 = 0;
	ei3 = 0;
	ei4 = 0;
	ei5 = 0;
	ei6 = 0;
	ei7 = 0;

	srand(time(NULL));

	for (x = Amin; x <= Amax; x = x + ddA) vAlpha.push_back(x);

	for (x = Tmin; x <= Tmax; x = x + ddT) vTemp.push_back(x);
 

	for (x = tpmin; x <= tpmax; x = x + ddtp) vtp.push_back(x);

	for (x = apmin; x <= apmax; x = x + ddap) vAp.push_back(x);

	for (x = tau1min; x <= tau1max; x = x + ddtau1) vtau1.push_back(x);

	for (x = tau2min; x <= tau2max; x = x + ddtau2) vtau2.push_back(x);



	printf("Pts =%i  %i %i %i %i %i  \n", vAlpha.size(), vTemp.size(), vAp.size(), vtp.size(), vtau1.size(), vtau2.size());
	printf("Size =%i \n", vAlpha.size()* vTemp.size()* vtp.size() *vAp.size() *vtau1.size() * vtau2.size());

	Lk = gen_matrix6(vAlpha.size(), vTemp.size(), vAp.size(), vtp.size(), vtau1.size(), vtau2.size());
 


	time_t t_inicial;
	time_t t_atual;
	t_inicial = time(NULL);


	std::vector<LikelihoodParameter> params;
	std::vector<LikelihoodParameter> paramsBuffer ;
	LikelihoodParameter max_likelihood_parameter(0,0,0,0,0,0);
	for (size_t i1 = 0; i1 < vAlpha.size() ; i1++)
	{

		 
		 
			{  

				params.clear(); //Mantem um cache de memoria do ciclo anterior
				paramsBuffer.clear();
				for (size_t i2 = 0; i2 < vTemp.size(); i2++)
				for (size_t i3 = 0; i3 < vAp.size(); i3++)
				for (size_t i4 = 0; i4 < vtp.size(); i4++)
					for (size_t i5 = 0; i5 < vtau1.size(); i5++)
						for (size_t i6 = 0; i6 < vtau2.size(); i6++)
						{
							//printf(" %i %i %i %i %i \n", i1,i2,i3,i4,i5);
							//Lk[i1][i2][i3][i4][i5][i6] = Likelihood_combined(vAlpha[i1], vTemp[i2], vAp[i3], vtp[i4], vtau1[i5], vtau2[i6]);
							//printf("%6.2f %6.2f %6.2f %6.2f  %6.2f %6.2f : %15.12e \n", vAlpha[i1], vTemp[i2], vAp[i3], vtp[i4], vtau1[i5], vtau2[i6], Lk[i1][i2][i3][i4][i5][i6]);
							 //                 _alpha,      _T,       _ap,       _tp,     _tau1,      _tau2
							paramsBuffer.emplace_back(vAlpha[i1], vTemp[i2], vAp[i3], vtp[i4], vtau1[i5], vtau2[i6]);
							if (paramsBuffer.size() >= 1024 * 16 )
							{
								computeParams(paramsBuffer);params.insert(params.end(), paramsBuffer.begin(), paramsBuffer.end()); paramsBuffer.clear();
							}

						}
				computeParams(paramsBuffer); params.insert(params.end(), paramsBuffer.begin(), paramsBuffer.end()); paramsBuffer.clear();
				 
				int iParam = 0;
				for (size_t i2 = 0; i2 < vTemp.size(); i2++)
				for (size_t i3 = 0; i3 < vAp.size(); i3++)
				for (size_t i4 = 0; i4 < vtp.size(); i4++)
					for (size_t i5 = 0; i5 < vtau1.size(); i5++)
						for (size_t i6 = 0; i6 < vtau2.size(); i6++)
						{
							const double i_priori1 = priori1(vAlpha[i1], vTemp[i2], vAp[i3], vtp[i4], vtau1[i5], vtau2[i6]);
							//params[iParam].result =  params[iParam].result;
							Lk[i1][i2][i3][i4][i5][i6] = i_priori1 *params[iParam].result;
																				
							if (params[iParam].result > Lmax)
							{
								Lmax =  params[iParam].result;
								max_likelihood_parameter = params[iParam];
							}
							iParam++;

						 	// if((iParam%300) ==0 ) printf("%6.2f %6.2f %6.2f %6.2f  %6.2f %6.2f : %15.12e \n", vAlpha[i1], vTemp[i2], vAp[i3], vtp[i4], vtau1[i5], vtau2[i6], Lk[i1][i2][i3][i4][i5][i6]);
						}
				 
			}

		t_atual = time(NULL);
		printf(" %i  %i   restam %4.1f  minutos \n", i1, vAlpha.size(), (vAlpha.size() - (i1 + 1))*(((t_atual - t_inicial) / 60.0) / (i1 + 1)));
	 
	}
 
 

 
 

	// Normaliza a preditivTemp, aplica a priori, e translada as matrizes


	lzmax = 0.0;
	for (size_t i1 = 0; i1 < vAlpha.size(); i1++) {
		for (size_t i2 = 0; i2 < vTemp.size(); i2++)
			for (size_t i3 = 0; i3 < vAp.size(); i3++)

			{
				for (size_t i4 = 0; i4 < vtp.size(); i4++)
					for (size_t i5 = 0; i5 < vtau1.size(); i5++)
						for (size_t i6 = 0; i6 < vtau2.size(); i6++)
							break;

							//Lk[i1][i2][i3][i4][i5][i6] = Lk[i1][i2][i3][i4][i5][i6] /  predi;

							//Lk[i1][i2][i3][i4][i5][i6] = Lk[i1][i2][i3][i4][i5][i6] * priori1(vAlpha[i1], vTemp[i2], vAp[i3], vtp[i4], vtau1[i5], vtau2[i6]);
				//	Limb[i1][i3][i4][i5] = Limb[i1][i3][i4][i5] * priori1(vAlpha[i1],vTemp[i3],vDeltaTp[i4], vAp[i5] );

			}
	}

	lzmax = 0.0;
	lmax1 = vAlpha[0];

	lmax3 = vtp[0];

	FILE *fmm;


 

	printf("Finalizado\n ");
	vAlpha.push_back(vAlpha.back() + ddA);
	vTemp.push_back(vTemp.back() + ddT);
	vAp.push_back(vAp.back() + ddap);
	vtp.push_back(vtp.back() + ddtp);
	vtau1.push_back(vtau1.back() + ddtau1);
	vtau2.push_back(vtau2.back() + ddtau2);



	predi = Integra_6n(Lk, vAlpha, vTemp, vAp, vtp, vtau1, vtau2);
	printf("predi Value %g \n", predi);


	fmm = fopen("Lat.dat", "w+"); 
	for (size_t i1 = 0; i1 < vAlpha.size()-1; i1++)
	{
		for (size_t i2 = 0; i2 < vTemp.size()-1; i2++)
		{
			const double soma = Integra_6n(Lk ,i1, i2, 0, vtp, vtau1, vtau2)/ predi;
			fprintf(fmm, "%f %f %g \n", vAlpha[i1], vTemp[i2], soma);
		} 
		fprintf(fmm, "\n");
	}
	fclose(fmm);


	//===========================================================================


	BidimensionalMap bm = BidimensionalMap([&](double x, double y)
	{
		const auto x0 = vTemp[0];
		const auto dx = vTemp[1] - vTemp[0];
		int i_temp = (x - x0) / dx;
        i_temp = std::min(i_temp, int(vTemp.size() - 2));

		const auto y0 = vAlpha[0];
		const auto dy = vAlpha[1] - vAlpha[0]; 
		int i_alpha = (y - y0) / dy;
		i_alpha = std::min(i_alpha, int(vAlpha.size() - 2));		
	 

		return Integra_6n(Lk, i_alpha, i_temp, 0, vtp, vtau1, vtau2) / predi;
	}, vTemp, vAlpha);


	bm.integrate_limite(0.5);

	fmm = fopen("Lta.dat", "w+");
	for (size_t i2 = 0; i2 < vTemp.size() - 1; i2++)		
	{
		for (size_t i1 = 0; i1 < vAlpha.size() - 1; i1++)
		{
			const double soma = Integra_6n(Lk, i1, i2, 0, vtp, vtau1, vtau2) / predi;
			fprintf(fmm, "%f %f %g \n",vTemp[i2] , vAlpha[i1],   soma);
		}
		fprintf(fmm, "\n");
	}
	fclose(fmm);


	//============================================================================




	double alpha_max = 0;
	double alpha_val = -1;
	fmm = fopen("Lalpha.dat", "w+");
	for (size_t i1 = 0; i1 < vAlpha.size()-1; i1++)
	{

		const double nsoma = Integra_6n(Lk, i1, vTemp, vAp, vtp, vtau1, vtau2)/ predi;
		fprintf(fmm, "%f  %g \n", vAlpha[i1],   nsoma  );
		if (nsoma > alpha_max)
		{
			alpha_max = nsoma;
			alpha_val = vAlpha[i1];
		}
	}
	fclose(fmm);

	//constroe a interpolacao de alpha


	std::vector<float> yAlpha;
	std::vector<float> xAlpha;
	for (size_t i1 = 0; i1 < vAlpha.size() - 1; i1++)
	{ 
		yAlpha.push_back(static_cast<float>(Integra_6n(Lk, i1, vTemp, vAp, vtp, vtau1, vtau2) / predi));
	    xAlpha.push_back((float)vAlpha[i1]);
	} 
	auto maxAlpha = max_value(xAlpha, yAlpha);


	std::vector<float> yTemp;
	std::vector<float> xTemp;
	for (size_t i2 = 0; i2 < vTemp.size() - 1; i2++)
	{
		yTemp.push_back(static_cast<float>(Integra_6n(Lk, vAlpha, i2, vAp, vtp, vtau1, vtau2) / predi));
		xTemp.push_back((float)vTemp[i2]);
	}
	auto maxTemp = max_value(xTemp, yTemp);


	std::pair<float, float> maxTau1 = std::pair<float, float>(0.0f, 0.0f);
	if (vtau1.size() > 2)
	{
		std::vector<float> yTau1;
		std::vector<float> xTau1;
		for (int i5 = 0; i5 < vtau1.size() - 1; i5++)
		{
			yTau1.push_back((float)(Integra_6n(Lk, vAlpha, vTemp, vAp, vtp, i5, vtau2) / predi));
			xTau1.push_back((float)vtau1[i5]);
		}
		maxTau1 = max_value(xTau1, yTau1);
	}


	std::pair<float, float> maxTp = std::pair<float, float>(0.0f, 0.0f);
	if (vtp.size() > 2)
	{
		std::vector<float> yTp;
		std::vector<float> xTp;
		for (int i4 = 0; i4 < vtp.size() - 1; i4++)
		{
			yTp.push_back((float)(Integra_6n(Lk, vAlpha, vTemp, vAp, i4, vtau1, vtau2) / predi));
			xTp.push_back((float)vtp[i4]);
		}
		maxTp = max_value(xTp, yTp);
	}






	double temp_max = 0;
	double temp_val = -1;
	fmm = fopen("LTemp.dat", "w+");	
	for (size_t i2 = 0; i2 < vTemp.size()-1; i2++)
	{	
		const double nsoma = Integra_6n(Lk, vAlpha, i2, vAp, vtp, vtau1, vtau2)/ predi;
		fprintf(fmm, "%f  %g \n", vTemp[i2], nsoma );
		if (nsoma > temp_max)
		{
			temp_max = nsoma;
			temp_val = vTemp[i2];
		}
	}
	fclose(fmm);




	double tau1_max = 0;
	double tau1_val = -1;
	fmm = fopen("LTau1.dat", "w+"); 
	for (size_t i5 = 0; i5 < vtau1.size() - 1; i5++)
	{
		const double nsoma = Integra_6n(Lk, vAlpha, vTemp, vAp, vtp, i5, vtau2)/predi;
		fprintf(fmm, "%f   %g \n", vtau1[i5], nsoma);
		if (nsoma > tau1_max)
		{
			tau1_max = nsoma;
			tau1_val = vtau1[i5];
		}
	}
	fclose(fmm);







	double tp_max = 0;
	double tp_val = -1;
	fmm = fopen("Ltp.dat", "w+");
	for (size_t i4 = 0; i4 < vtp.size() - 1; i4++)
	{
		const double nsoma = Integra_6n(Lk, vAlpha, vTemp, vAp, i4, vtau1, vtau2) / predi;
		fprintf(fmm, "%f   %g \n", vtp[i4], nsoma);
		if (nsoma > tp_max)
		{
			tp_max = nsoma;
			tp_val = vtp[i4];
		}
	}
	fclose(fmm);




	fmm = fopen("parameters.dat", "w+");
	fprintf(fmm, "Max Likehood Parameter:%g\n alpha %6.2f\n T %6.2f\n Tp %6.2f\nAmplitute %6.2f\n Tau 1 %6.2f\n Tau 2%6.2f \n", max_likelihood_parameter.result,
		max_likelihood_parameter.alpha, max_likelihood_parameter.T,
		max_likelihood_parameter.tp,
		max_likelihood_parameter.ap,
		max_likelihood_parameter.tau1, max_likelihood_parameter.tau2);

	fprintf(fmm, "best-fit  alpha  %f %f \n", alpha_val, maxAlpha.first);
	fprintf(fmm, "best-fit  T      %f %f\n", temp_val,  maxTemp.first);
	fprintf(fmm, "best-fit  Tau1   %f %f\n", tau1_val,  maxTau1.first);
	fprintf(fmm, "best-fit  Tbust  %f %f\n", tp_val, maxTp.first);
	fprintf(fmm, "Integration    %g \n", predi);
	fclose(fmm);



	fmm = fopen("tau12.dat", "w+"); 

	for (size_t i5 = 0; i5 < vtau1.size()-1; i5++)
	{
		for (size_t i6 = 0; i6 < vtau2.size()-1; i6++)

		{
			double soma = 0;
			for (size_t i1 = 0; i1 < vAlpha.size()-1; i1++)
				for (size_t i2 = 0; i2 < vTemp.size()-1; i2++)
					for (size_t i3 = 0; i3 < vAp.size()-1; i3++)
						for (size_t i4 = 0; i4 < vtp.size()-1; i4++)
						{
							soma += ddap * ddtp * ddA* ddT*   Lk[i1][i2][i3][i4][i5][i6];
						}
			fprintf(fmm, "%f %f %g \n", vtau1[i5], vtau2[i6], soma);
		}
		fprintf(fmm, "\n");
	}
	fclose(fmm);


 





	fmm = fopen("LTT.dat", "w+"); 
	for (size_t i2 = 0; i2 < vTemp.size()-1; i2++)
	{ 
		for (size_t i5 = 0; i5 < vtau1.size()-1; i5++)
		{
			const double nsoma = Integra_6n(Lk, vAlpha, i2, vAp, vtp, i5, vtau2)/predi;
			fprintf(fmm, "%f %f %g  \n", vTemp[i2], vtau1[i5], nsoma );

		}
		fprintf(fmm, "\n");
	}
	fclose(fmm);


	printf("Likelihood = %g \n", Lmax);





	//fmm = fopen("rcol.dat", "w+");
	//for (double e = 1; e < 50; e += 0.2)
	//{
	//	for (double Te = 0.1; Te < 20; Te += 0.2)
	//	{
	//		fprintf(fmm, " %f %f %e \n", e, Te, etabarK(e)*Rcol(e, max_likelihood_parameter.alpha, 0, Te, max_likelihood_parameter.ap, max_likelihood_parameter.tp, max_likelihood_parameter.tau1, max_likelihood_parameter.tau2));
	//	}
	//	fprintf(fmm, "\n");
	//}
	//fclose(fmm); 


	//fmm = fopen("temp.dat", "w+");	 
	//{
	//	for (double tt = 0; tt < 20; tt += 0.02)
	//	{
	//		double Te = Temp(tt, max_likelihood_parameter.T, max_likelihood_parameter.ap, max_likelihood_parameter.tp, max_likelihood_parameter.tau1, max_likelihood_parameter.tau2);
	//		fprintf(fmm, " %f %f \n", tt, Te);
	//	}
	// 
	//}
	//fclose(fmm);
	//

	//fmm = fopen("alpha.dat", "w+");
	//{
	//	for (double tt = 0; tt < 20; tt += 0.02)
	//	{
	//		double Te = get_alpha(max_likelihood_parameter.alpha, tt, max_likelihood_parameter.T, max_likelihood_parameter.ap, max_likelihood_parameter.tp, max_likelihood_parameter.tau1, max_likelihood_parameter.tau2);
	//		fprintf(fmm, " %f %f \n", tt, Te);
	//	}

	//}
	//fclose(fmm);

	return 0;

}






