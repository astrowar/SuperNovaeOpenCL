#pragma once

#include<functional>
#include<vector>
#include <climits>
#include <algorithm>
#include <float.h>

class BidimensionalMap
{
public:
	std::vector<double> x;
	std::vector<double> y;
	std::vector<std::vector<double>> val;
	BidimensionalMap( std::function<double(double, double)>  func, std::vector<double> _x, std::vector<double> _y):x(_x),y(_y)
	{		
		size_t nx = x.size();
		size_t ny = y.size();
		val = std::vector<std::vector<double>>(ny);
		for (size_t j = 0;j< ny;++j)
		{
			val[j]= std::vector<double>(nx);
			for (size_t i = 0; i < nx; ++i)
			{
				auto r0 = x[i];
				auto r1 = y[j];
				val[j][i] = func(r0,r1);
			}
		}
	}

	double max_value()
	{
		double vmax = -FLT_MAX;
		for(auto vy : val)
			for (auto vx : vy)
			{
				vmax = std::max(vmax, vx);
			}
		return vmax;
	}

	double min_value()
	{
		double vmin = FLT_MAX;
		for (auto vy : val)
			for (auto vx : vy)
			{
				vmin = std::min(vmin, vx);
			}
		return vmin;
	}

	double integrate(float min_value)
	{
		double s = 0.0;
		for (auto &vy : val)
			for (auto &vx : vy)
			{
				if (vx >= min_value)
				{
					s = s + vx;
				}
			}
		return s;
	}

	double integrate_limite(double percentil)
	{
		 if (min_value()<0)
		 {
			 throw "only positive define functions";
		 }

		 double umax = max_value();  //maior Z
		 double smax = integrate(0.0); //integral de todos os valores >= 0 

		 double z0 = 0.0;
		 double z1 = umax;
		 double zn = umax / 2.0;
		for(int loop =0; loop < 20;loop++)
		{
			zn = 0.5*(z1 + z0); //mid
			double sn = integrate(zn);

			double dz = 0.25*(z1 - z0);

			if (sn <= percentil*smax)
			{
				z1 = zn + dz;
			}
			else
			{
				z0 = zn- dz ;
			}
			printf("Bisec Inner %f  => %f (%f)\n", zn, sn , percentil*smax  );
		}

		for(double   z = 0.0f * umax ; z< 1.5*umax ; z = z +  0.1*umax)
		{
			double s = integrate(z);
			 printf("Inner %f  => %f \n",z, s);
		}
		return  0.1;
	}

};
