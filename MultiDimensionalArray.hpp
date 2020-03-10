
#ifndef MULTIDIMENSIONALARRAY_HPP
#define MULTIDIMENSIONALARRAY_HPP
#include <utility>
#include <vector>


template <  typename N>
class Interval
{
	const double start;
	const double end;
	const double delta;
	const size_t num_divs;
	Interval(double _start, double _end , double _delta): start(_start),end(_end), delta(_delta), num_divs(size_t(fabs(end - start) / delta))
	{
		 
	}
	int i(double x )
	{
		int ii =  (x - start) / delta;
		if (ii < 0) return 0;
		if (ii > num_divs-1) return num_divs-1;
		return ii;
	}

};


template <  int N>
class MultiDimensionalArray
{
public:
	std::vector< MultiDimensionalArray< N - 1> > data;
	std::vector< double> x;

	MultiDimensionalArray(std::vector<double> _x ) : x(std::move(_x))
	{
		
	}
};


template< > class MultiDimensionalArray<1>
{
public:
	std::vector<double> data;
	std::vector< double> x;
	MultiDimensionalArray(std::vector<double> _x) : x(std::move(_x))
	{
		data.resize(x.size(), 0.0);
	}
};



template < int N ,typename... Rest>
MultiDimensionalArray<N>  make_mdarray(std::vector< double> x, Rest... rest)
{
 
	if constexpr (sizeof...(rest) > 0)
	{		
		auto m = MultiDimensionalArray<N>(x);
		//m.data.reserve(x.size());
		for (auto _x : x)
		{
			m.data.push_back(make_mdarray<N - 1>(rest...));
		}
		return m;
	}
	else
	{
		return MultiDimensionalArray<1>(x);
	}
}


template <typename F, int N,  typename... Rest>
void  iterator_function( F func,  MultiDimensionalArray<N> &m , Rest... rest)
{
	if constexpr (N > 1)
	{
	for(int i =0;i<m.x.size();++i)
	  
	  {
		  double x = m.x[i];
		  iterator_function(func,m.data[i], rest..., x);
	  }
	}
	else
	{		
		for (auto x : m.x)
		{
			func(  rest... , x );
		}
	}
}

//MultiDimensionalArray<1> make_mdarray(std::vector< double> x )
//{	
//	return MultiDimensionalArray<1>(x);
//}


#endif // MULTIDIMENSIONALARRAY_HPP
