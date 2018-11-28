
#include <vector>


template <size_t N>
class MultiDimensionalArray
{
public:
	std::vector< MultiDimensionalArray< N - 1> > data;
	std::vector< float> x;


	MultiDimensionalArray(std::vector< float> _x ) : x(_x)
	{

	}
};


template<> class MultiDimensionalArray<1>
{
public:
	std::vector<double> data;
	std::vector< float> x;
	
};


template <size_t N,  typename... Rest>
MultiDimensionalArray<N> make(std::vector< float> x, Rest... rest)
{
	auto m = MultiDimensionalArray<N>(x);
	for (auto _x : x)
	{
		m.data.push_back(make<N - 1>(rest...));
	}
	return m;
}


template <typename... Rest>
MultiDimensionalArray<1> make(std::vector< float> x )
{	
	return MultiDimensionalArray<1>(x);
}

