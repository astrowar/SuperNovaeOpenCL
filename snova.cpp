// snova.cpp : Defines the entry point for the console application.
//

 

#define _CRT_SECURE_NO_WARNINGS


#include <cstdio>
#include <cstdlib>
#include <cmath>


#include <vector>
#include <ctime>
#include <algorithm>


#include  "IntegrateFunctions.hpp"


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
#define Cn    1/(4 * 3.1416)
#define lnVk   14.56
#define lnVimb  15.73
#define delt 1.0//2.3
#define timeK 20.0//10.43
#define timeIMB 16.0 //5.9



 


#define real double


using namespace std;

#define Integra( r, ff , x, x1 , x2 , dx )  for(x = (x1) ; x <= (x2) ; x+=0  ){ r = r + (dx) * ( ff )/6.0 ; x+=(dx)/2.0 ; r = r + 4*(dx) * ( ff )/6.0; x+=(dx)/2.0 ;r = r + (dx) * ( ff )/6.0; };


#define Dump(  ff , x, x1 , x2 , dx )  for(x = (x1) ; x <= (x2) ;x+=(dx)   ){ printf(">%f %f \n",x, ( ff ) );};

#define H  6.0

//regiao de calculo da massa dos distibuicoes T
#define Tmin 2.0
#define Tmax 8.0
#define ddT  0.1

#define sigma_tp 1.0

//  tempo
#define tpmin 0.2 
#define tpmax 0.3
#define ddtp 1.0 


// alpha
#define Amin 0.2
#define Amax 20.0
#define ddA  0.2

//  eps
#define epsmin 1.0
#define epsmax 70.0  
#define ddeps  5.0  


// tau1

#define  tau1min 0.1
#define  tau1max 15.0   
#define  ddtau1  0.2  


// tau2

#define  tau2min 0.0
#define  tau2max 4.0
#define  ddtau2  4.0


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


//Kamikande events:


double tk[17] = { 0.0, 0.0, 0.107, 0.303, 0.324, 0.507, 0.686, 1.541, 1.728, 1.915, 9.219, 10.433, 12.439, 17.641, 20.257, 21.355, 23.814 }; // times of events;
double Ek[17] = { 0.0, 20.0, 13.5, 7.5, 9.2, 12.8, 6.3, 35.4, 21, 19.8, 8.6, 13, 8.9, 6.5, 5.4, 4.6, 6.5 };   // energy of events;
double Sigmak[17] = { 0.0, 2.9, 3.2, 2.0, 2.7, 2.9, 1.7, 8, 4.2, 3.2, 2.7, 2.6, 1.9, 1.6, 1.4, 1.3, 1.6 }; // standard deviation by events;
double Bk[17] = { 1, 1.6e-5, 1.9e-3, 2.9e-2, 1.2e-2, 2.1e-3, 3.7e-2, 4.5e-5, 8.2e-5, 1.5e-5,            // detector's noise;
1.5e-2, 1.9e-3, 1.6e-2, 3.8e-2, 2.9e-2, 2.8e-2, 3.8e-2 };


// IMB events:

double timb[9] = { 0.0, 0.0, 0.412, 0.650, 1.141, 1.562, 2.684, 5.010, 5.582 };
double Eimb[9] = { 0,38,37,28,39,36,36,19,22 }; // energy of events;
double Sigmaimb[9] = { 0,7,7 ,6,7,9, 6,5,5 };      // standard deviation by events;
double Bimb[9] = { 0,0,0,0,0,0,0,0,0 };         // detector's noise;


												// Defined functions and matrixes:



double**  gen_matrix2(int n, int m)
{
	double **d;

	d = (double**)malloc(sizeof(double*)*n);

	for (int i3 = 0; i3<n; i3++)
	{
		d[i3] = (double*)malloc(sizeof(double)*m);
	}

	return d;
}
double*** gen_matrix3(int n, int m, int m2)
{
	double ***d;

	d = (double***)malloc(sizeof(double**)*n);


	for (int i3 = 0; i3<n; i3++)
	{
		d[i3] = gen_matrix2(m, m2);
	}

	return d;
}

double**** gen_matrix4(int n, int m, int m2, int m3)
{
	double ****d;

	d = (double****)malloc(sizeof(double***)*n);

	for (int i3 = 0; i3<n; i3++)
	{
		d[i3] = gen_matrix3(m, m2, m3);
	}

	return d;
}

double***** gen_matrix5(int n, int m, int m2, int m3, int m4)
{
	double *****d;

	d = (double*****)malloc(sizeof(double****)*n);

	for (int i4 = 0; i4<n; i4++)
	{
		d[i4] = gen_matrix4(m, m2, m3, m4);
	}

	return d;
}

double****** gen_matrix6(int n, int m, int m2, int m3, int m4, int m5)
{
	double ******d;

	d = (double******)malloc(sizeof(double*****)*n);

	for (int i5 = 0; i5<n; i5++)
	{
		d[i5] = gen_matrix5(m, m2, m3, m4, m5);
	}

	return d;
}

double******* gen_matrix7(int n, int m, int m2, int m3, int m4, int m5, int m6)
{
	double *******d;

	d = (double*******)malloc(sizeof(double******)*n);

	for (int i6 = 0; i6 < n; i6++)
	{
		d[i6] = gen_matrix6(m, m2, m3, m4, m5, m6);
	}

	return d;
}



double  Integra_1(double *y, vector<double> x1)
{
	double soma = 0.0;
	int i;
	int ifinal = x1.size();
	for (i = 0; i< ifinal - 2; i += 2)
	{
		soma += (y[i] + 4 * y[i + 1] + y[i + 2])*(x1[i + 2] - x1[i]);
	}
	return soma / 6.0;
}


double  Integra_2(double **yf, vector<double> x1, vector<double> x2)
{
	double soma = 0.0;
	double   ya, yb, yc;
	int i;
	int ifinal = x1.size();
	ya = Integra_1(yf[0], x2);
	for (i = 0; i< ifinal - 2; i += 2)
	{
		yb = Integra_1(yf[i + 1], x2);
		yc = Integra_1(yf[i + 2], x2);
		soma += (ya + 4 * yb + yc)*(x1[i + 2] - x1[i]);
		ya = yc;
	}
	return soma / 6.0;
}

double  Integra_3(double ***y, vector<double> x1, vector<double> x2, vector<double> x3)
{
	double soma = 0.0;
	double  ya, yb, yc;
	int i;
	int ifinal = x1.size();
	ya = Integra_2(y[0], x2, x3);

	for (i = 0; i< ifinal - 2; i += 2)
	{
		yb = Integra_2(y[i + 1], x2, x3);
		yc = Integra_2(y[i + 2], x2, x3);
		soma += (ya + 4 * yb + yc)*(x1[i + 2] - x1[i]);
		ya = yc;
	}
	return soma / 6.0;


}

double  Integra_4(double ****y, vector<double> x1, vector<double> x2, vector<double> x3, vector<double> x4)
{
	double soma = 0.0;
	double   ya, yb, yc;
	int i;
	int ifinal = x1.size();
	ya = Integra_3(y[0], x2, x3, x4);

	for (i = 0; i< ifinal - 2; i += 2)
	{
		yb = Integra_3(y[i + 1], x2, x3, x4);
		yc = Integra_3(y[i + 2], x2, x3, x4);
		soma += (ya + 4 * yb + yc)*(x1[i + 2] - x1[i]);
		ya = yc;
	}
	return soma / 6.0;



}

double  Integra_5(double *****y, vector<double> x1, vector<double> x2, vector<double> x3, vector<double> x4, vector<double> x5)
{
	double soma = 0.0;
	double   ya, yb, yc;
	int i;
	int ifinal = x1.size();
	ya = Integra_4(y[0], x2, x3, x4, x5);

	for (i = 0; i< ifinal - 2; i += 2)
	{
		yb = Integra_4(y[i + 1], x2, x3, x4, x5);
		yc = Integra_4(y[i + 2], x2, x3, x4, x5);
		soma += (ya + 4 * yb + yc)*(x1[i + 2] - x1[i]);
		ya = yc;
	}
	return soma / 6.0;



}


double  Integra_6(double ******y, vector<double> x1, vector<double> x2, vector<double> x3, vector<double> x4, vector<double> x5, vector<double> x6)
{
	double soma = 0.0;
	double   ya, yb, yc;
	int i;
	int ifinal = x1.size();
	ya = Integra_5(y[0], x2, x3, x4, x5, x6);
	for (i = 0; i< ifinal - 2; i += 2)
	{
		yb = Integra_5(y[i + 1], x2, x3, x4, x5, x6);
		yc = Integra_5(y[i + 2], x2, x3, x4, x5, x6);
		soma += (ya + 4 * yb + yc)*(x1[i + 2] - x1[i]);
		ya = yc;
	}
	return soma / 6.0; 
}


double gaussian(double  x, double x0, double sigma) {

	if (fabs(x - x0) < sigma * H) return (1.0 / (sqrt(2 * 3.1415928* sigma*sigma))) * exp(-0.5*pow((x - x0) / sigma, 2.0));

	return 0.0;

}

double gaussian_n(double  x, double x0, double sigma) {

	if (fabs(x - x0) < sigma * H) return  exp(-0.5*pow((x - x0) / sigma, 2.0));

	return 0.0;

}




double obs_sample_statistic_err(double mu)
{
	double soma;
	soma = 0;


	return soma;
}

double obs_sample_statistic(double mu, double dm)
{
	double soma = 0;

	return soma;
}

//*************************************************** Defined functions ******************************************************************

// Function type Press - Schecter
double Temp(double tu, double T, double ap, double tp, double tau1, double tau2) {

	double durac = 1.0;

	// return T;
	if (tu > tp )  return  T*  exp(-min(tp / tau1, 10.0));
	return  T*  exp(-min(tu / tau1, 10.0));
 

	if (tu <= tp) return T * exp(-tu / tau1) + 0.01;
	return T * ap*  exp(-(tu - tp) / tau2) + 0.01;

	double TR = T;
	if (tu <= tp) return  std::max( T + tau1 * tu, 0.01);	 
	return  std::max(T + tau1 * tp, 0.01);


	return   T * exp(-min(tu / tau1, 10.0));
	double T1 = T * exp(-min(tu / tau1, 10.0));
	if (tu <= tp) return   T1 + 0.01;
	double T2 = ap * T * exp(-min((tu - tp) / tau2, 10.0));
	return  T2 + 0.01;


	//Linear
	//double T1 = T * exp(-std::min( tu / tau1  , 10.0 )) ;
	//if (tu <= tp) return   T1  + 0.01;
	//double T2 = ap * T * exp(-std::min((tu-tp) / tau2 , 10.0 ));
	//return  T2 + 0.01;
	


	//Mudanca

	if (tu > tp + 4*tau1 )   return 0.01;
	if (tu >= tp) return T *  exp(-(tu - tp) / tau1);
	return T;	
	



	  // return T;

		  if(tp >= 0.0 && tp <= 2.0) return T;

		  if(tp >= 2.0 && tp <= 2.0 + tu) return T * exp(-tu / (1.0 * tau1))  + 0.1;

		  if (tp >=2.0 + tu  && tp <= 2.0 + durac + tu ) return  ap * T *  exp(-(tu - 1.0 * tp) / (1.0 * tau2) );

		  if (tp >= 2.0 + durac + tu) return 0.01;


	double t_saida = 0 ;

	//	t_saida = T * exp(-tu / 4.0 * tau1)  + 0.1 ;

	//	if (tu >= tp) t_saida +=  T * 100 * pow(ap,2) * exp(-(tu - 1.0 * tp) / 1.0 * tau2) ;

	//	return t_saida;

	//	if (tu > tp) and if(tu < tp + durac)

	//		return T;

	if (tu >= tp) //and if(tu < tp + durac + durac) 

		return t_saida = T * exp(-tu / 1.0 * tau1) + 0.1;

	//	if(tu >= tp + durac) 

	//		return T;



	if (tu >= tp + durac) t_saida += T * ap * exp(-(tu - 1.0 * tp) / 1.0 * tau2);

	return t_saida;



}


real get_alpha(double alpha, real tu, real T, real ap, real tp, real tau1, real tau2)
{
 

	double tt = Temp(tu, T, ap, tp, tau1, tau2);
	return alpha * pow(tt / T, -tau2);



	if (tu > tp) return alpha * exp(min(tp / tau2, 10.0));
	return alpha * exp(min(tu / tau2, 10.0));

	return alpha;





	//5*exp(-0.8*x)+3.2 -y  w l
	if (tu > tp) return  ap* alpha;
	return   alpha ;


 

}


double r(double ap, double tp, double tau1, double tau2) {

	//return exp(-Delta_tp/(tp+0.1));

	return 1.0;

	//return exp(- (tu - tp) / 15 * tau2) + 0.1;

}

double conv_gauss_fast( double x,   double xa0, double sa,  double xb0, double sb)
{
	double rr = (0.1994711309046004*
		exp((-0.5*(xa0*xa0) + 1.*xa0*xb0 - 0.5*(xb0* xb0)) / ((sa*sa) + (sb*sb)))
		*sa*sb*	erf(((sb*sb)*(0.7071067811865475*x - 0.7071067811865475*xa0) +
		(sa*sa)*(0.7071067811865475*x - 0.7071067811865475*xb0)) /
			(sa*sb*sqrt((sa* sa) + (sb*sb))))) /
			(sqrt((sa*sa))*sqrt((sb*sb))*sqrt((sa*sa) + (sb* sb)));
	return rr;
}

double cov_gauss(double x1, double q1, double x2, double q2)
{
	// calcula a covolucao de duas gaussianas
	double a1, a2;
	double ss;
	double u;
	a1 = std::max(x1 - H*q1, x2 - H*q2);
	a2 = std::min(x1 + H*q1, x2 + H*q2);
	if (a2<a1) return 0;
	//ss = 0;
	//Integra(ss, (gaussian(u, x1, q1) * gaussian(u, x2, q2)), u, a1, a2, (a2 - a1) / 4.0);

	ss = conv_gauss_fast(a2, x1, q1, x2,q2 ) - conv_gauss_fast(a1, x1, q1,x2, q2);

	return ss;

}


double f_fermi(double eps, double tu, double T, double ap, double tp, double tau1, double tau2) {

	// if(eps/T > 5){return 0.0;}
	//  if(eps/T < -5){return 1;}
	
	double To = Temp(tu, T, ap, tp, tau1, tau2);
	 
	if (To < 0.1) return 0;
	double f_ret =      1.0 / (exp((eps + Q) / To) + 1.0);
	return f_ret;

}
//  Corrections function

double kappa(double eps) {
	double a, b, c, E;
	E = eps + Q;
	a = 1.0 - Q / E;
	b = 1.0 - (2 * Q / E);
	c = (pow(Q, 2) - pow(me, 2)) / pow(E, 2);

	if ((b + c) < 0.0) return 0.0;

	return a * sqrt(b + c);
}

// Rate's neutrino - cooling component=

double Rcol(double eps, double alpha, double tu, double T, double ap, double tp, double _tau1, double _tau2) {

	//double fm;
	double kp;
	double saida = 0;
	//  if ((alpha == last_alpha) &&(T==last_T) &&(eps==last_eps) &&(tp==last_tp))         { return last_out; }

	double fm = f_fermi(eps, tu, T, ap, tp, _tau1, _tau2);
 
	 
	if (fm == 0.0) { return 0.0; }
	kp = kappa(eps);
	if (kp == 0.0) return 0.0;
	saida = (1.22e-5) * pow(alpha, 2) * Meff * (pow((eps + Q), 4)) * fm * kp * pow(r(ap, tp, _tau1, _tau2), 2);
	if (saida <0)  printf("ERROR Rcol \n");
	return saida;

}


//----------------------------  Kamiokande  -------------------------------------------------



double  noiseK(double eps)
{
	double gs;
	//  gs = effK*(gaussianK(eps, 6 ,1)+ 0.01875*gaussianK(eps,20,10));
	gs = effK* (gaussian(eps, 6, 1) + 0.01);
	// <---------------------- noise parameter


	//  if (eps > 45 ) return (0.0 + gs)*pk ;
	//  if (eps < 5) return (0.0001 +gs)*pk ;
	//  if (eps < 20 ) return  (gs+ 0.0016 * (eps/20.0))*pk;
	//  return   (0.0016 *(eps - 20.0)/( 30.0 ) + gs)*pk  ;
	// if (eps < 10.0) return 0.01;


	return gs;
}

double etabarK(double eps) {

	double c;

	c = 0.95*(1.0 - exp(-pow((eps - Q) / 9.3, 4)));

	if (c < 0.0) return 0.0;
	return c;




	// return 80*atan(eps/80);
	//return tanh(eps / 10);
	// return 0.9;
	double eK;
	double y1, y2, x1, x2;
	y1 = 0;
	y2 = 0.9;
	//  x1 = 4.0;
	//x2 = 10.0;
	x1 = 10.0;
	x2 = 12.0;


	if (eps < 4.0) return 0.0;
	if (eps > 10.0) return 0.9;
	return (y2 - y1) / (x2 - x1) + y1;
	//  return 0.9 ;


	if (eps < 5.0) { return 0.0; }
	eK = (0.35) * atan(0.49 * eps - 4.0) + 0.43;
	if (eK < 0.0) { return 0.0; }
	return eK;
}

// Step function

double StepK(double eps) {


	if (eps >  5.0) {
		return 1;
	}
	else { return 0.0; }
}
//     IMB



double LikelihoodK(double alpha, double T, double tp, double ap, double tau1, double tau2) {

	int  i, j;
	double soma;


	double termo1, termo2;
	double prod;
	double eps;
	double ti;

	double e1, e2;

	soma = 0.0;
	Meff = Mk;

	termo1 = 0.0;
	//         Integra( termo1 ,   0.73 * StepK( eps - 5  )* (Rcol(alpha , T, eps, tp ) + noiseK( eps )) , eps  ,  epsmin, epsmax, ddeps )  ;



	for (ti = 0; ti <= 30.0; ti = ti + ddtp)
	{
		Integra(termo1, ddtp * etabarK(eps) *(Cn*Rcol(eps, alpha, ti, T, ap, tp, tau1, tau2) + 1.0 * noiseK(eps)), eps, epsmin, epsmax, ddeps);
	}

	if (termo1<0)printf("termo1= %f \n", termo1);

	prod = 0.0;
	//printf("-----------------------------------------------\n");
	for (i = 1; i <= 12; i++)
	{
		if (i == 6) continue;


		// if(Ek[i] < 10.0) continue ;
		termo2 = 0.0;

		e1 = std::max(epsmin, Ek[i] - H * Sigmak[i]);
		e2 = std::min(epsmax, Ek[i] + H * Sigmak[i]);

		Integra(termo2, StepK(eps)* Cn * lnVk * gaussian(eps, Ek[i], Sigmak[i])*Cn*(Rcol(eps, alpha, tk[i], T, ap, tp, tau1, tau2) + 1.0*noiseK(eps)), eps, e1, e2, Sigmak[i] / 5.0);

		//      for( ti= tk[i] -3*sigma_tp ; ti <= tk[i]+3*sigma_tp ; ti= ti+ddtp  )


		//Integra( termo2 , ddtp* 0.9 * etabarK(eps)*  gaussian( ti, tk[i] , sigma_tp  )* gaussian( eps, Ek[i] , Sigmak[i]  )*(Rcol( eps, alpha,  T,  ap, ti, Delta_tp ) )    , eps  ,  epsmin , Ek[i]+3*Sigmak[i] , ddeps )  ;


		// if (termo2<0 )
		//printf("termo2= %g \n", termo2 );


		prod = prod + log(termo2);
	}
	// soma = exp(-23.0*Integral_rateK(alpha,T,tp)) * ProdKa(alpha,T,tp);
	// if (prod < 1e-7  )printf("prod= %f \n", prod );
	double e_term = -1.0 * delt * timeK * termo1 + prod;
	soma = exp(e_term);


	//printf("integral rate=%e \n", Integral_rateK(alpha,T,tp) );
	//printf("Exp integral rate=%e \n",  exp(-Integral_rateK(alpha,T,tp))   );
	// printf("prod K= %e \n",ProdKa(alpha,T,tp)  );
	// printf("Like K , a=%e T=%e = %e\n", alpha,T, soma  );

	if (soma > Lmax) Lmax = soma;
	return soma;
}





//-------------------------------------------  IMB  ----------------------------------------------------------------

double  noiseIMB(double eps)

{
	double gs = 1.0;
	gs = effIMB* (gaussian(eps, 6, 1) + 0.001);


	return pimb*gs;
}

double etabarIMB(double eps) {
	double c;
	//           double eIMB;


	c = (1.0 - 3.0*exp(-(eps - Q) / 16.0));

	if (c < 0.0) return 0.0;
	return c;


	c = 1.0;
	if (eps > 30) c = 5.0;//2.885;
	if (eps > 60.0) { return 1.0; }
	if (eps < 20.0) { return 0.0; }


	return  c*(0.1*sqrt(6 * eps + 1) - 1.0);

}




double StepIMB(double eps) {


	if (eps > 19.0) {
		//if (eps > 17.11){
		return 1;
	}
	else { return 0.0; }

}

double LikelihoodIMB(double alpha, double T, double tp, double ap, double tau1, double tau2) {

	int  i;
	double soma;

	double termo1, termo2;
	double prod;
	double eps;
	double ti;


	soma = 0.0;
	Meff = Mimb;

 
	prod = 1.0;
	for (i = 1; i <= 8; i++)
	{
		//  if (i==7) continue;
		// if (i==8) continue;
		 

		// if (i==8) continue;
		// if(Ek[i] < 10.0) continue ;
		termo2 = 0.0;

		//for( ti= timb[i] - 3 * sigma_tp ; ti <= timb[i] + 3*sigma_tp ; ti= ti + ddtp  )
		// termo2=0.0 ;
		//Integra(termo2, StepIMB(eps) *Cn*lnVimb* gaussian(eps, Eimb[i], Sigmaimb[i])*(Cn*Rcol(eps, alpha, timb[i], T, ap, tp, tau1, tau2) + noiseIMB(eps)), eps, epsmin, epsmax, Sigmaimb[i] / 2.0);
		Integra(termo2, StepIMB(eps) *Cn*lnVimb* gaussian(eps, Eimb[i], Sigmaimb[i])*(Cn*Rcol(eps, alpha, timb[i], T, ap, tp, tau1, tau2) + noiseIMB(eps)), eps, epsmin, epsmax, ddeps);

		prod = prod * termo2;
	}

	termo1 = 0.0;

	for (ti = 0; ti <= 20; ti = ti + ddtp)
	{
		Integra(termo1, ddtp*etabarIMB(eps) * (Cn*Rcol(eps, alpha, ti, T, ap, tp, tau1, tau2) + noiseIMB(eps)), eps, epsmin, epsmax, ddeps);
	}
	if (termo1<0)printf("termo1= %f \n", termo1);


	double sExpTerm = -1.0 * delt * timeIMB * termo1;
	soma = exp(sExpTerm) * prod;
	if (soma > Lmax) Lmax = soma;
	return soma;

}

//--------------------------------------------------- Combined likelihood ---------------------------------------------------------------

double Likelihood_combined(double alpha, double T, double ap, double tp, double tau1, double tau2) {

	//double Lcomb =  LikelihoodK(alpha, T, ap, tp, tau1, tau2) * LikelihoodIMB(alpha, T, ap, tp, tau1, tau2);
	double Lcomb = LikelihoodIMB(alpha, T, ap, tp, tau1, tau2)  ;

	return Lcomb;
	//   return LikelihoodK(alpha,T,ap, tp,tau1,tau2);

	// return LikelihoodIMB(alpha,T,ap,Delta_tp);

}


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


int build_Program();

int main() {

	 

	FILE *arquivo;
	double i, x, y, z;

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
	//Limb = gen_matrix4(  vAlpha.size(), vTemp.size(), vDeltaTp.size() ,vAp.size() );


	//L_s12 = gen_matrix6(vAlpha.size(), vTemp.size(), vAp.size(), vtp.size(), vtau1.size(), vtau2.size());

	//L_s12= gen_matrix5( vAp.size(), vAlpha.size(), vTemp.size(), vtau1.size(),vtau2.size() );


	time_t t_inicial;
	time_t t_atual;
	t_inicial = time(NULL);


	std::vector<LikelihoodParameter> params;
	std::vector<LikelihoodParameter> paramsBuffer ;
	LikelihoodParameter max_likelihood_parameter(0,0,0,0,0,0);
	for (int i1 = 0; i1 < vAlpha.size() ; i1++)
	{

		 
		 
			{  

				params.clear(); //Mantem um cache de memoria do ciclo anterior
				paramsBuffer.clear();
				for (int i2 = 0; i2 < vTemp.size(); i2++)
				for (int i3 = 0; i3 < vAp.size(); i3++)
				for (int i4 = 0; i4 < vtp.size(); i4++)
					for (int i5 = 0; i5 < vtau1.size(); i5++)
						for (int i6 = 0; i6 < vtau2.size(); i6++)
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
				for (int i2 = 0; i2 < vTemp.size(); i2++)
				for (int i3 = 0; i3 < vAp.size(); i3++)
				for (int i4 = 0; i4 < vtp.size(); i4++)
					for (int i5 = 0; i5 < vtau1.size(); i5++)
						for (int i6 = 0; i6 < vtau2.size(); i6++)
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
	for (int i1 = 0; i1 < vAlpha.size(); i1++) {
		for (int i2 = 0; i2 < vTemp.size(); i2++)
			for (int i3 = 0; i3 < vAp.size(); i3++)

			{
				for (int i4 = 0; i4 < vtp.size(); i4++)
					for (int i5 = 0; i5 < vtau1.size(); i5++)
						for (int i6 = 0; i6 < vtau2.size(); i6++)
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
	for (int i1 = 0; i1 < vAlpha.size()-1; i1++)
	{
		for (int i2 = 0; i2 < vTemp.size()-1; i2++)

		{
 

			const double soma = Integra_6n(Lk ,i1, i2, 0, vtp, vtau1, vtau2);
			fprintf(fmm, "%f %f %g \n", vAlpha[i1], vTemp[i2], soma);

		} 
		fprintf(fmm, "\n");
	}
	fclose(fmm);







	double alpha_max = 0;
	double alpha_val = -1;
	fmm = fopen("Lalpha.dat", "w+");
	for (int i1 = 0; i1 < vAlpha.size()-1; i1++)
	{
		double soma = 0; 
			for (int i2 = 0; i2 < vTemp.size()-1; i2++)
			  for (int i3 = 0; i3 < vAp.size()-1; i3++)
				for (int i4 = 0; i4 < vtp.size()-1; i4++)
					for (int i5 = 0; i5 < vtau1.size()-1; i5++)
						for (int i6 = 0; i6 < vtau2.size()-1; i6++)							 
							soma += ddT* ddap * ddtp * ddtau1 * ddtau2 * Lk[i1][i2][i3][i4][i5][i6];

		const double nsoma = Integra_6n(Lk, i1, vTemp, vAp, vtp, vtau1, vtau2)/ predi;
		fprintf(fmm, "%f  %g %g\n", vAlpha[i1],   nsoma , soma);
		if (nsoma > alpha_max)
		{
			alpha_max = nsoma;
			alpha_val = vAlpha[i1];
		}
	}
	fclose(fmm);

	double temp_max = 0;
	double temp_val = -1;
	fmm = fopen("LTemp.dat", "w+");	
	for (int i2 = 0; i2 < vTemp.size()-1; i2++)
	{
		double soma = 0;
		for (int i1 = 0; i1 < vAlpha.size()-1; i1++)
			for (int i3 = 0; i3 < vAp.size()-1; i3++)
				for (int i4 = 0; i4 < vtp.size()-1; i4++)
					for (int i5 = 0; i5 < vtau1.size()-1; i5++)
						for (int i6 = 0; i6 < vtau2.size()-1; i6++)
							soma += ddA * ddap * ddtp * ddtau1 * ddtau2 * Lk[i1][i2][i3][i4][i5][i6];
		
		
		const double nsoma = Integra_6n(Lk, vAlpha, i2, vAp, vtp, vtau1, vtau2)/ predi;
		fprintf(fmm, "%f  %g  %g\n", vTemp[i2], nsoma, soma);
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
	for (int i5 = 0; i5 < vtau1.size() - 1; i5++)
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


	fmm = fopen("parameters.dat", "w+");
	fprintf(fmm, "Max Likehood Parameter:%g\n alpha %6.2f\n T %6.2f\n Tp %6.2f\nAmpliture %6.2f\n Tau 1 %6.2f\n Tau 2%6.2f \n", max_likelihood_parameter.result,
		max_likelihood_parameter.alpha, max_likelihood_parameter.T,
		max_likelihood_parameter.tp,
		max_likelihood_parameter.ap,
		max_likelihood_parameter.tau1, max_likelihood_parameter.tau2);

	fprintf(fmm, "best-fit  alpha %f \n", alpha_val);
	fprintf(fmm, "best-fit  T     %f \n", temp_val);
	fprintf(fmm, "best-fit  Tau1  %f \n", tau1_val);
	fclose(fmm);



	fmm = fopen("tau12.dat", "w+"); 

	for (int i5 = 0; i5 < vtau1.size()-1; i5++)
	{
		for (int i6 = 0; i6 < vtau2.size()-1; i6++)

		{
			double soma = 0;
			for (int i1 = 0; i1 < vAlpha.size()-1; i1++)
				for (int i2 = 0; i2 < vTemp.size()-1; i2++)
					for (int i3 = 0; i3 < vAp.size()-1; i3++)
						for (int i4 = 0; i4 < vtp.size()-1; i4++)
						{
							soma += ddap * ddtp * ddA* ddT*   Lk[i1][i2][i3][i4][i5][i6];
						}
			fprintf(fmm, "%f %f %g \n", vtau1[i5], vtau2[i6], soma);
		}
		fprintf(fmm, "\n");
	}
	fclose(fmm);


 





	fmm = fopen("LTT.dat", "w+"); 
	for (int i2 = 0; i2 < vTemp.size()-1; i2++)
	{ 
		for (int i5 = 0; i5 < vtau1.size()-1; i5++)
		{
			const double nsoma = Integra_6n(Lk, vAlpha, i2, vAp, vtp, i5, vtau2);
			fprintf(fmm, "%f %f %g  \n", vTemp[i2], vtau1[i5], nsoma );

		}
		fprintf(fmm, "\n");
	}
	fclose(fmm);


	printf("Likelihood = %g \n", Lmax);

	fmm = fopen("rcol.dat", "w+");
	for (double e = 1; e < 50; e += 0.2)
	{
		for (double Te = 0.1; Te < 20; Te += 0.2)
		{
			fprintf(fmm, " %f %f %e \n", e, Te, etabarK(e)*Rcol(e, max_likelihood_parameter.alpha, 0, Te, max_likelihood_parameter.ap, max_likelihood_parameter.tp, max_likelihood_parameter.tau1, max_likelihood_parameter.tau2));
		}
		fprintf(fmm, "\n");
	}
	fclose(fmm); 


	fmm = fopen("temp.dat", "w+");	 
	{
		for (double tt = 0; tt < 20; tt += 0.02)
		{
			double Te = Temp(tt, max_likelihood_parameter.T, max_likelihood_parameter.ap, max_likelihood_parameter.tp, max_likelihood_parameter.tau1, max_likelihood_parameter.tau2);
			fprintf(fmm, " %f %f \n", tt, Te);
		}
	 
	}
	fclose(fmm);
	

	fmm = fopen("alpha.dat", "w+");
	{
		for (double tt = 0; tt < 20; tt += 0.02)
		{
			double Te = get_alpha(max_likelihood_parameter.alpha, tt, max_likelihood_parameter.T, max_likelihood_parameter.ap, max_likelihood_parameter.tp, max_likelihood_parameter.tau1, max_likelihood_parameter.tau2);
			fprintf(fmm, " %f %f \n", tt, Te);
		}

	}
	fclose(fmm);

	return 0;

}






