
#pragma OPENCL EXTENSION cl_amd_fp64 : enable

#define real double
#define number float

#if CL_VERSION_1_2
#else
  #include <math.h>
real min(real a, real b);
real max(real a, real b);
  int printf(const char *   format, ...);

#endif



#define NP 400

#define Q 1.29
#define me 0.511
#define eta 0.0
#define Mk  2.14     //  Kamiokande detector mass (kton).
#define Mimb 6.8     // IMB detector mass (kton).
#define Mbaksan 0.28 // Baksan detector mass (kton).
#define fk (1.0)
#define fimb (0.9055)
#define fb (1.0)
#define Cn    (1.0/(4.0 * 3.1416))
#define lnVk   14.56
#define lnVimb  15.73
#define delt 1.0  //2.3
//#define delt 2.3
#define timeK 20.0  //10.43
#define timeIMB 16.0   //5.9
#define timeEnd 25.0   




 

#define Integra( r, ff , x, x1 , x2 , dx )  for(x = (x1) ; x <= (x2) ; x+=0  ){ r = r + (dx) * ( ff )/6.0 ; x+=(dx)/2.0 ; r = r + 4*(dx) * ( ff )/6.0; x+=(dx)/2.0 ;r = r + (dx) * ( ff )/6.0; };
 

#define Dump(  ff , x, x1 , x2 , dx )  for(x = (x1) ; x <= (x2) ;x+=(dx)   ){ printf(">%f %f \n",x, ( ff ) );};

#define H  6.0
#define NEPS 128

#define epsmin 0.1
#define epsmax 80.0 // 50.0
#define ddeps   0.1 // 5.0



// Noise's parameters

#define pk     1.0
#define pimb   1.0
#define pb     1.0
#define effK   1e-5
#define effIMB 1e-5


//number PARMS[5];
//(FILE*)FILES[5];
//int NUMPARAM=5;





 
//Kamikande events:


//const number    tk[17] = { 0.0, 0.0, 0.107, 0.303, 0.324, 0.507, 0.686, 1.541, 1.728, 1.915, 9.219, 10.433, 12.439, 17.641, 20.257, 21.355, 23.814 }; // times of events;
//const number    Ek[17] = { 0.0, 20.0, 13.5, 7.5, 9.2, 12.8, 6.3, 35.4, 21, 19.8, 8.6, 13, 8.9, 6.5, 5.4, 4.6, 6.5 };   // energy of events;
//const number    Sigmak[17] = { 0.0, 2.9, 3.2, 2.0, 2.7, 2.9, 1.7, 8, 4.2, 3.2, 2.7, 2.6, 1.9, 1.6, 1.4, 1.3, 1.6 }; // standard deviation by events;
//const number    Bk[17] = { 1, 1.6e-5, 1.9e-3, 2.9e-2, 1.2e-2, 2.1e-3, 3.7e-2, 4.5e-5, 8.2e-5, 1.5e-5,            // detector's noise;
//1.5e-2, 1.9e-3, 1.6e-2, 3.8e-2, 2.9e-2, 2.8e-2, 3.8e-2 };


// IMB events:

//const number    timb[9] = { 0.0, 0.0, 0.412, 0.650, 1.141, 1.562, 2.684, 5.010, 5.582 };
//const number    Eimb[9] = { 0,38,37,28,39,36,36,19,22 }; // energy of events;
//const number     Sigmaimb[9] = { 0,7,7 ,6,7,9, 6,5,5 };      // standard deviation by events;
//const number    Bimb[9] = { 0,0,0,0,0,0,0,0,0 };         // detector's noise;
//

												// Defined functions and matrixes:

 
 

 

#define Integra( r, ff , x, x1 , x2 , dx )  for(x = (x1) ; x <= (x2) ; x+=0  ){ r = r + (dx) * ( ff )/6.0 ; x+=(dx)/2.0 ; r = r + 4*(dx) * ( ff )/6.0; x+=(dx)/2.0 ;r = r + (dx) * ( ff )/6.0; };



number gaussian(number  x, number x0, number sigma) {

	if (fabs(x - x0) < sigma * H) return (1.0 / (sqrt((number)2 * 3.1415928* sigma*sigma))) * exp(-0.5*pow((x - x0) / sigma, (number)2.0));


	//return (1.0 / (sqrt((number)2 * 3.1415928))) * exp(-0.5*pow(H, (number)2.0));

	return 0.0;

}

number gaussian_n(number  x, number x0, number sigma) {

	if (fabs(x - x0) < sigma * H) return  exp(-0.5*pow((x - x0) / sigma, (number)2.0));

	return 0.0;

}




number obs_sample_statistic_err(number mu)
{
	number soma;
	soma = 0;


	return soma;
}

number obs_sample_statistic(number mu, number dm)
{
	number soma = 0;

	return soma;
}

//*************************************************** Defined functions ******************************************************************

 
typedef struct PressSchecter
{
	number ap,  tp,  tau1 , tau2;
}PressSchecter;


// Function type Press - Schecter
number Temp(number tu, number T, PressSchecter tmp) {

	  
	if (tu > tmp.tp)
	{ 
		number u = (tu - tmp.tp) / (4.0*tmp.tau1);
		if (u > 10.0) return 0.0;
		return  T* exp(-u);
	}
	return  T;

}



number get_alpha(double alpha,   number tu  )
{
	//return 1.0 / pow(1.0 + tu / tau1, 1.5);
	return alpha;
}

number r(number tu, number T, PressSchecter tmp   ) {

	//return exp(-1.0 * min(10.0, tu / (2.0*tmp.tau1)));
	return 1.0;

}

number conv_gauss_fast(number x, number xa0, number sa, number xb0, number sb)
{
	number rr = (0.1994711309046004*
		exp((-0.5*(xa0*xa0) + 1.*xa0*xb0 - 0.5*(xb0* xb0)) / ((sa*sa) + (sb*sb)))
		*sa*sb*	erf(((sb*sb)*(0.7071067811865475*x - 0.7071067811865475*xa0) +
		(sa*sa)*(0.7071067811865475*x - 0.7071067811865475*xb0)) /
			(sa*sb*sqrt((sa* sa) + (sb*sb))))) /
			(sqrt((sa*sa))*sqrt((sb*sb))*sqrt((sa*sa) + (sb* sb)));
	return rr;
}

number cov_gauss(number x1, number q1, number x2, number q2)
{
	// calcula a covolucao de duas gaussianas
	number a1, a2;
	number ss;
	// number u ;
	a1 = max(x1 - H*q1, x2 - H*q2);
	a2 = min(x1 + H*q1, x2 + H*q2);
	if (a2<a1) return 0;
	ss = conv_gauss_fast(a2, x1, q1, x2, q2) - conv_gauss_fast(a1, x1, q1, x2, q2);

	return ss;

}


number f_fermi(number eps,   number T )
{
	 if (T < 0.1) return 0;
	//return     1.0 / (exp((E) / T) + 1.0);
	return 1.0 / (exp((eps + Q) / T) + 1.0);

}
//  Corrections function

number kappa(number eps) {
	
	number a, b, c, E;
	E = eps + Q;
	a = 1.0 - Q / E;
	b = 1.0 - (2 * Q / E);
	c = (pow((number)Q, (number)2) - pow((number)me, (number)2)) / pow((number)E, (number)2);

	if ((b + c) < 0.0) return 0.0;

	return a * sqrt(b + c);
}

// Rate's neutrino - cooling component=

number Rcol(number eps, number alpha,   number T , number MMeff ,number radius ) {

	number fm;
	number kp;
	number saida = 0.0;
	 

	fm = f_fermi(eps,  T   );
	if (fm <= 0.0) { return 0.0; }
	 
	number alpha_t =   (get_alpha(alpha, T));
	 

	kp = kappa(eps);
	if (kp <= 0.0 ) return 0.0;
	saida = (1.22e-5) * alpha_t*alpha_t * MMeff * (pow((number)(eps + Q), (number)4.0)) * fm * kp * pow(radius, (number)2.0);
	if (saida <0)  printf((__constant char *)"ERROR Rcol \n");
	return max(saida, (number)0.0);

}


//----------------------------  Kamiokande  -------------------------------------------------



number  noiseK(number eps)
{
	 
	return  effK* (gaussian(eps, 6.0, 1.0) +0.001 )  ;

	
}

number etabarK(number eps) {

	number c;

	c = 0.95*(1.0 - exp(-pow((number)((eps - Q) / 9.3), (number)4)));

	if (c < 0.0) return 0.0;
	return c;


}

// Step function

number StepK(number eps) {


	if (eps >  5.0) {
		return 1;
	}
	else { return 0.0; }
}
//     IMB


 

 
 

 
 
real LikelihoodK(number alpha, number T, PressSchecter tmp, real *LMax, __local number eps_value[NEPS], __local number      etabarK_value[NEPS], __local number      noiseK_value[NEPS]) {

	
 
	const number    tk[17] = { 0.0, 0.0, 0.107, 0.303, 0.324, 0.507, 0.686, 1.541, 1.728, 1.915, 9.219, 10.433, 12.439, 17.641, 20.257, 21.355, 23.814 }; // times of events;
	const number    Ek[17] = { 0.0, 20.0, 13.5, 7.5, 9.2, 12.8, 6.3, 35.4, 21, 19.8, 8.6, 13, 8.9, 6.5, 5.4, 4.6, 6.5 };   // energy of events;
	const number    Sigmak[17] = { 0.0, 2.9, 3.2, 2.0, 2.7, 2.9, 1.7, 8, 4.2, 3.2, 2.7, 2.6, 1.9, 1.6, 1.4, 1.3, 1.6 }; // standard deviation by events;
	const number    Bk[17] = { 1, 1.6e-5, 1.9e-3, 2.9e-2, 1.2e-2, 2.1e-3, 3.7e-2, 4.5e-5, 8.2e-5, 1.5e-5,            // detector's noise;
		1.5e-2, 1.9e-3, 1.6e-2, 3.8e-2, 2.9e-2, 2.8e-2, 3.8e-2 };


	//if (get_local_id(0) == 0)
	//{
	//	number de = (epsmax - epsmin) / NEPS;
	//	for (int i = 0; i < NEPS; i++)
	//	{
	//		number eps = epsmin + de *i;
	//		etabarK_value[i] = etabarK(eps);
	//		noiseK_value[i] = noiseK(eps);
	//		eps_value[i] = eps;
	//	}
	//}
	//barrier(CLK_LOCAL_MEM_FENCE)


	int  i, j;
	number soma;


	real termo1, termo2;
	real prod;
	number eps;
	number ti;

	number e1, e2;

	soma = 0.0;

	number jddtp =  max(0.01, tmp.tau1 /5.0);
	number time_end = min(tmp.tp + 3*H *tmp.tau1, timeEnd);
	termo1 = 0.0;


	 number de = (epsmax - epsmin) / NEPS;
	//{
	//	number Tj = Temp(0, T, tmp); 
	//	//Integra(termo1, (etabarK(eps)*(Cn*Rcol(eps + Q, alpha, Tj, Mk) + noiseK(eps))), eps, epsmin, epsmax, ddeps);
	//	 for (int i = 0; i < NEPS; ++i)
	//	 {
	//	 	eps = eps_value[i];
	//	 	termo1 += de*  (etabarK_value[i] * (Cn*Rcol(eps , alpha, Tj, Mk) + noiseK_value[i]));
	//	 }

	//	termo1 = termo1 * tmp.tp;
	//}
	//for (ti = tmp.tp; ti <= time_end; ti = ti + jddtp)
	//{
	//	number Tj = Temp(ti , T,tmp);
	//	 for (int i = 0; i < NEPS; ++i)
	//	 {
	//	 	eps = eps_value[i];
	//	 	termo1 += de* jddtp*(etabarK_value[i] * (Cn*Rcol(eps , alpha, Tj, Mk) + noiseK_value[i]));
	//	 }
	//	//Integra(termo1, jddtp*(etabarK(eps)*(Cn*Rcol(eps + Q,alpha,Tj,Mk) + noiseK(eps))), eps, epsmin, epsmax, ddeps);		 
	//}
	// 
 //
	//{
	//	number termo_ending = 0.0;
	//	number Tj = Temp(time_end, T, tmp);
	//	//Integra(termo1, (etabarK(eps)*(Cn*Rcol(eps + Q, alpha, Tj, Mk) + noiseK(eps))), eps, epsmin, epsmax, ddeps);
	//	for (int i = 0; i < NEPS; ++i)
	//	{
	//		eps = eps_value[i];
	//		termo_ending += de*  (etabarK_value[i] * (Cn*Rcol(eps , alpha, Tj, Mk) + noiseK_value[i]));
	//	}

	//	termo_ending = termo_ending * fabs(time_end - timeEnd);
	//	termo1 += termo_ending;
	//}





	 number ddtp = tmp.tp / 2.0;
	 for (ti = 0; ti <= tmp.tp ; ti = ti + ddtp)
	 {
		 number Tj = Temp(0, T, tmp);
		 number radius = r(0, Tj, tmp);
		 eps = eps_value[0];
		 number y1 = (etabarK_value[0] * (Cn*Rcol(eps, alpha, Tj, Mk, radius) + noiseK_value[0]));
		 for (int i = 0; i < NEPS-1; ++i)
		 {
			 eps = eps_value[i];
			 //number y1 = (etabarK_value[i] * (Cn*Rcol(eps, alpha, Tj, Mk , radius) + noiseK_value[i]));
			 number y2 = (etabarK_value[i+1] * (Cn*Rcol(eps+de, alpha, Tj, Mk, radius) + noiseK_value[i+1]));
			 termo1 += de* (ddtp)* (y1 + y2) / 2.0;
			 y1 = y2;
		 }

		 // Integra(termo1, ddtp * etabarK(eps) *(Cn*Rcol(eps, alpha, Tj, Mk) + 1.0 * noiseK(eps)), eps, epsmin, epsmax, ddeps);
		  
	 }


	ddtp = 0.2;
	//for (ti = tmp.tp ; ti <= max(0* timeK,   (tmp.tp + 4.0*tmp.tau1)) ; ti = ti + ddtp)
	for (ti = tmp.tp; ti <= timeEnd; ti = ti + ddtp)
	{
		number Tj = Temp(ti, T, tmp);
		number radius = r(ti, Tj, tmp);
		eps = eps_value[0];
		number y1 = (etabarK_value[0] * (Cn*Rcol(eps, alpha, Tj, Mk, radius) + noiseK_value[0]));

		for (int i = 1; i < NEPS; ++i)
		{
			eps = eps_value[i];
			//number y1 = (etabarK_value[i] * (Cn*Rcol(eps, alpha, Tj, Mk, radius) + noiseK_value[i]));
			number y2 = (etabarK_value[i  ] * (Cn*Rcol(eps  , alpha, Tj, Mk,radius) + noiseK_value[i ]));
			termo1 += de* (ddtp)* (y1 + y2) / 2.0;
			y1 = y2;
		//	termo1 += de* ddtp*(etabarK_value[i] * (Cn*Rcol(eps, alpha, Tj, Mk) + noiseK_value[i]));
			
		}
		
		int i = 0;
		//Integra(termo1, ddtp * etabarK(eps) *(Cn*Rcol(eps, alpha, Tj, Mk) + 1.0 * noiseK(eps)), eps, epsmin, epsmax, ddeps);
		 
	}


	prod = 0.0;	 
	for (i = 1; i <= 12; i++)
	{
		if (i == 6) continue;  
		termo2 = 0.0;
		e1 = max(epsmin, Ek[i] - H * Sigmak[i]);
		e2 = min(epsmax, Ek[i] + H * Sigmak[i]);
		number dSigma = Sigmak[i] / 2.0;
		number Tj = Temp(tk[i], T, tmp);
		number radius = r(tk[i], Tj, tmp);
		//Integra(termo2, RDet(eps) , eps, e1, e2, dSigma);
		
		//Integra(termo2, (StepK(eps)* Cn * lnVk * gaussian(eps, Ek[i], Sigmak[i])*Cn*(Rcol(eps, alpha,   Tj,Mk) + 1.0*noiseK(eps))), eps, e1, e2, dSigma);
		Integra(termo2, (StepK(eps)* Cn * lnVk * gaussian(eps, Ek[i], Sigmak[i])*Cn*(Rcol(eps, alpha, Tj, Mk,radius) + 1.0*noiseK(eps))), eps, e1, e2, dSigma);
		prod = prod + log(max( 1e-20, Bk[i] + termo2));
		//prod = prod  + log (Bk[i] + termo2 + 0.001);
	}
	 
	//double e_term = -1.0 * delt * termo1 + prod;
	//soma = exp(e_term);
	double e_term = -1.0   *delt *  5*(termo1)      ;
	soma = exp(e_term + prod)    ;
 
	 

	if (soma > (*LMax)) *LMax = soma;
	return soma;
}





//-------------------------------------------  IMB  ----------------------------------------------------------------

number  noiseIMB(number eps)

{
	number gs;
	gs = effIMB* (gaussian(eps, 6.0, 2.0) + 0.001);
	return pimb*gs;
}

number etabarIMB(number eps) {
	number c;
	 


	c = (1.0 - 3.0*exp(-(eps - Q) / 16.0));

	if (c < 0.0) return 0.0;
	return c;


	//c = 1.0;
	//if (eps > 30) c = 5.0;//2.885;
	//if (eps > 60.0) { return 1.0; }
	//if (eps < 20.0) { return 0.0; }
	//return  c*(0.1*sqrt(6.0 * eps + 1) - 1.0);

}




number StepIMB(number eps) {


	if (eps > 19.0) {
	 
		return 1.0;
	}
    return 0.0; 

}

number perParticleIMB(int i ,number eps , number alpha, number T, PressSchecter tmp )
{
	const number    timb[9] = { 0.0, 0.0, 0.412, 0.650, 1.141, 1.562, 2.684, 5.010, 5.582 };
	const number    Eimb[9] = { 0,38,37,28,39,36,36,19,22 }; // energy of events;
	const number     Sigmaimb[9] = { 0,7,7 ,6,7,9, 6,5,5 };      // standard deviation by events;
	 const number    Bimb[9] = { 0,0,0,0,0,0,0,0,0 };         // detector's noise;


	number Tj = Temp(timb[i], T, tmp);
	number radius = r(timb[i], Tj, tmp);
	number y1 = StepIMB(eps);
	//number y1 = etabarIMB(eps);
	number y2 = Cn*lnVimb* gaussian(eps, Eimb[i], Sigmaimb[i]);
	number y3 = (Cn*Rcol(eps, alpha,  Tj,  Mimb, radius) + noiseIMB(eps));

	//printf((__constant char *)"%g %g %g \n",y1,y2,y3);

	return y1*y2*y3;
}


number LikelihoodIMB(number alpha, number T, PressSchecter tmp, real* LMax, __local number eps_value[NEPS], __local number      etabarIMB_value[NEPS], __local number      noiseIMB_value[NEPS]) {


	const number    timb[9] = { 0.0, 0.0, 0.412, 0.650, 1.141, 1.562, 2.684, 5.010, 5.582 };
	const number    Eimb[9] = { 0,38,37,28,39,36,36,19,22 }; // energy of events;
	const number    Sigmaimb[9] = { 0,7,7 ,6,7,9, 6,5,5 };      // standard deviation by events;
	const number    Bimb[9] = { 0,0,0,0,0,0,0,0,0 };         // detector's noise;




	int  i;
	number soma;

	number termo1, termo2;
	number prod;
	number eps;
	number ti;


	soma = 0.0;
	//Meff = Mimb;

	termo1 = 0.0;
	number jddtp = max(0.01, min(0.5 / tmp.tau1, 1.0));

	number time_end = timeEnd;


	//for (ti = 0; ti <= time_end; ti = ti + jddtp)
	//{
	//	number Tj = Temp(ti, T, tmp);
	//	//Integra(termo1, jddtp * StepIMB(eps) * (Cn*Rcol(eps, alpha, ti, T, ap, tp, tau1, tau2, Mimb) + noiseIMB(eps)), eps, epsmin, epsmax, ddeps / 2.0);
	//	Integra(termo1,  etabarIMB(eps) * (Cn*Rcol(eps + Q, alpha,  Tj, Mimb) + noiseIMB(eps)), eps, epsmin, epsmax, ddeps );
	//}
	//termo1 = termo1 *jddtp;


	number de = (epsmax - epsmin) / NEPS;
	for (ti = 0; ti <= tmp.tp; ti = ti + tmp.tp / 2.0)
	{
		number Tj = Temp(0, T, tmp);
		number radius = r(timb[i], Tj, tmp);
		eps = eps_value[0];
		number y1 = (etabarIMB_value[0] * (Cn*Rcol(eps, alpha, Tj, Mk, radius) + noiseIMB_value[0]));
		for (int i = 0; i < NEPS; ++i)
		{
			eps = eps_value[i];
			//number y1 = (etabarIMB_value[i] * (Cn*Rcol(eps, alpha, Tj, Mk, radius) + noiseIMB_value[i]));
			number y2 = (etabarIMB_value[i + 1] * (Cn*Rcol(eps + de, alpha, Tj, Mk, radius) + noiseIMB_value[i + 1]));
			termo1 += de* (tmp.tp / 2.0)* (y1 + y2) / 2.0;
			y1 = y2;
			//termo1 += de* tmp.tp / 2.0*(etabarIMB_value[i] * (Cn*Rcol(eps, alpha, Tj, Mimb) + noiseIMB_value[i]));
		}
	 
	}


	number ddtp = 0.2;
	//for (ti = tmp.tp; ti <=  max( 0*timeIMB,  (tmp.tp + 4.0*tmp.tau1)); ti = ti + ddtp)
	for (ti = tmp.tp; ti <= timeEnd; ti = ti + ddtp)
	{
		number Tj = Temp(ti, T, tmp);
		number radius = r(ti, Tj, tmp);

		eps = eps_value[0];
		number y1 = (etabarIMB_value[0] * (Cn*Rcol(eps, alpha, Tj, Mk, radius) + noiseIMB_value[0]));
		 

		for (int i = 0; i < NEPS; ++i)
		{
			eps = eps_value[i];
			//number y1 = (etabarIMB_value[i] * (Cn*Rcol(eps, alpha, Tj, Mk, radius) + noiseIMB_value[i]));
			number y2 = (etabarIMB_value[i + 1] * (Cn*Rcol(eps + de, alpha, Tj, Mk, radius) + noiseIMB_value[i + 1]));
			termo1 += de* (ddtp)* (y1 + y2) / 2.0;

			y1 = y2;
			//termo1 += de* ddtp*(etabarIMB_value[i] * (Cn*Rcol(eps, alpha, Tj, Mimb) + noiseIMB_value[i]));
		}
	 
	}



	prod = 0.0;
	

	for (i = 1; i <= 8; i++)
	{
 
		termo2 = 0.0; 
		double eps_range_begin = Eimb[i] - H*Sigmaimb[i];
		double eps_range_end = Eimb[i] + H*Sigmaimb[i];
		eps_range_begin = max(eps_range_begin, epsmin);
		eps_range_end = min(eps_range_end, epsmax);

		Integra(termo2, perParticleIMB(i, eps ,     alpha,   T, tmp) , eps, eps_range_begin, eps_range_end, Sigmaimb[i] / 3.0);
		//Integra(termo2, StepIMB(eps) *Cn*lnVimb* gaussian(eps, Eimb[i], Sigmaimb[i])*(Cn*Rcol(eps, alpha, timb[i], T, ap, tp, tau1, tau2,Mimb) + noiseIMB(eps)), eps, epsmin, epsmax, Sigmaimb[i] / 2.0);

		// prod = prod * termo2;
		number dprod = Bimb[i] + termo2;
		dprod = max(dprod, (number)1e-20);
		prod = prod + log(dprod);
		//prod = prod + log( 1e-4 +  Bimb[i] + termo2)  ;
	}
	 
	

 
	//soma = exp(-1.0 * delt * timeIMB * termo1) * prod;
	soma = exp(-1.0 * delt *  timeIMB* termo1 + prod);
	//soma = exp(-termo1 +  prod);

	if (soma > (*LMax) ) *LMax = soma;
	return soma;

}

//--------------------------------------------------- Combined likelihood ---------------------------------------------------------------


typedef struct LikelihoodParameter
{
	
	real alpha, T, ap, tp, tau1, tau2;
	real result;
} LikelihoodParameter;

real Likelihood_combined(real alpha, real T, real ap, real tp, real tau1, real tau2, __local number eps_value[NEPS], __local number      etabarK_value[NEPS], __local number      noiseK_value[NEPS], __local number      etabarIMB_value[NEPS], __local number      noiseIMB_value[NEPS])
{

	PressSchecter psch = { ap, tp, tau1,tau2 };
	real LMax = 0.0;

	real LK= LikelihoodK(alpha, T, psch, &LMax,eps_value,etabarK_value, noiseK_value);
	 
	 real LIMB = LikelihoodIMB(alpha, T, psch, &LMax, eps_value, etabarIMB_value, noiseIMB_value);
	//return LIMB;
	  return LK*LIMB;
	
 
}


 kernel void LikelihoodList(global const LikelihoodParameter *inputParams, global real *output,   const int nArgs)
{
	 __local number      eps_value[NEPS];
	 __local number      etabarK_value[NEPS];
	 __local number      noiseK_value[NEPS];
	 __local number      etabarIMB_value[NEPS];
	 __local number      noiseIMB_value[NEPS];

	int tid = get_global_id(0);
	//if ( tid ==0) printf((__constant char *)" id \n");



	if (get_local_id(0) == 0)
	{
		number de = (epsmax - epsmin) / NEPS;
		for (int i = 0; i < NEPS; i++)
		{
			number eps = epsmin + de *i;
			etabarK_value[i] = etabarK(eps);
			noiseK_value[i] =  noiseK(eps);
			etabarIMB_value[i] = etabarIMB(eps);
			noiseIMB_value[i] =  noiseIMB(eps);
			eps_value[i] = eps;
		}
	}
	barrier(CLK_LOCAL_MEM_FENCE);

 


	if (tid < nArgs)
	{
		real res = Likelihood_combined(inputParams[tid].alpha, inputParams[tid].T,   inputParams[tid].ap, inputParams[tid].tp, inputParams[tid].tau1, inputParams[tid].tau2,   eps_value , etabarK_value , noiseK_value, etabarIMB_value, noiseIMB_value);
		output[tid] = res;
	}
}



kernel void vectorAdd(global const real *inputA, global const real *inputB, global real *output, const real x)
{
	 
	output[get_global_id(0)] = inputA[get_global_id(0)] + inputB[get_global_id(0)] + x;
}
