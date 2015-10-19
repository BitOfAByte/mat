#ifndef MATH_IMPL_H
#define MATH_IMPL_H

#ifndef MAT_IMPL_NO_MATH
# include <math.h>
# include <float.h>
#endif

namespace matrix
{
	template<>
	inline float mathutils::pi<float>()
	{
		return 3.1415926535897932384626433832795;
	}

	template<>
	inline double mathutils::pi<double>()
	{
		return 3.1415926535897932384626433832795;
	}

	template<typename Tp>
	inline Tp mathutils::todeg(const Tp& x)
	{
		return x * 180. / pi<Tp>();
	}
	
	template<typename Tp>
	inline Tp mathutils::torad(const Tp& x)
	{
		return x / 180. * pi<Tp>();
	}
	
	template<typename Tp>
	inline void mathutils::swap(Tp& a, Tp& b)
	{
		Tp c(a);
		a = b;
		b = c;
	}
	
	template<typename Tp>
	inline Tp mathutils::abs(const Tp& x)
	{
		if(x < 0)
			return -x;
		else
			return x;
	}

#ifndef MAT_IMPL_NO_MATH

	template<>
	inline float mathutils::epsilon(const float& dumb)
	{
		return FLT_EPSILON;
	}
	
	template<>
	inline double mathutils::epsilon(const double& dumb)
	{
		return DBL_EPSILON;
	}
	
	template<>
	inline float mathutils::sqrt<float>(const float& x)
	{
		return ::sqrt(x);
	}
	
	template<>
	inline double mathutils::sqrt<double>(const double& x)
	{
		return ::sqrt(x);
	}
	
	template<>
	inline float mathutils::sin<float>(const float& x)
	{
		return ::sin(x);
	}
	
	template<>
	inline double mathutils::sin<double>(const double& x)
	{
		return ::sin(x);
	}
	
	template<>
	inline float mathutils::cos<float>(const float& x)
	{
		return ::cos(x);
	}
	
	template<>
	inline double mathutils::cos<double>(const double& x)
	{
		return ::cos(x);
	}
#endif
}
#endif

