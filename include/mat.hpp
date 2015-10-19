#ifndef MAT_H
#define MAT_H

#include <math.hpp>
#include <impl.hpp>

namespace matrix
{
	typedef MAT_SCALAR_TYPE scalar_t;
	
	enum actions {
		T = 1
	};
	
	template<typename Tp = scalar_t, typename U = mathutils>
	class mat
	{
		typedef impl<Tp, U> I;
	public:
		mat() : m_cx(0), m_cy(0), m_data(NULL) {}

		mat(Tp* data, const unsigned int& cx, const unsigned int& cy, const bool& transpose = false) : m_cx(0), m_cy(0), m_data(NULL)
		{
			set(data, cx, cy, transpose);
		}

		mat(const mat& m) : m_cx(0), m_cy(0), m_data(NULL)
		{
			*this = m;
		}
		
		static mat rot(const Tp& theta)
		{
			mat t;
			t.set(NULL, 2, 2);
			I::rot(theta);
			return t;
		}
		
		static mat id(const int& cx = -1)
		{
			mat t;
			t.set(NULL, cx, cx, false, false);
			I::id(data(), cx);
			return t;
		}

		~mat()
		{
			set(NULL, 0, 0,false);
		}

		inline const unsigned int& cx() const
		{
			return m_cx;
		}

		inline const unsigned int& cy() const
		{
			return m_cy;
		}

		inline const Tp& at(const unsigned int& x, const unsigned int& y) const
		{
			return m_data[y * m_cx + x];
		}

		inline Tp& at(const unsigned int& x, const unsigned int& y)
		{
			return m_data[y * m_cx + x];
		}

		inline Tp* data()
		{
			return m_data;
		}

		inline const Tp* data() const
		{
			return m_data;
		}

		mat& operator=(const mat& m)
		{
			return set(m.data(), m.cx(), m.cy());
		}

		mat& operator+=(const mat& x)
		{
			assert(cx() == x.cx() && cy() == x.cy());
			I::add(data(), m_cx, m_cy, data(), x.data());
			return *this;
		}

		mat& operator-=(const mat& x)
		{
			assert(cx() == x.cx() && cy() == x.cy());
			I::sub(data(), m_cx, m_cy, data(), x.data());
			return *this;
		}

		//m+m;
		mat operator+(const mat& x)
		{
			mat t(x);
			t += *this;
			return t;
		}

		//m-m;
		mat operator-(const mat& x) const
		{
			mat t(x);
			t -= *this;
			return t;
		}

		//m*x
		mat operator*(const Tp& x) const
		{
			mat t(*this);
			t *= x;
			return t;
		}

		//m*=x
		mat& operator*=(const Tp& x)
		{
			I::smul(data(), m_cx, m_cy, data(), x);
			return *this;
		}

		//m*m;
		mat operator*(const mat& x) const
		{
			//s must mul x
			assert(cx() == x.cy());
			mat t;
			t.set(NULL, x.cx(), cy());
			I::mul(t.data(), data(), cx(), cy(), x.data(), x.cx(), x.cy());
			return t;
		}

		//m/m;
		mat operator/(const mat& x) const
		{
			//x must be square
			assert(x.cx() == x.cy());
			assert(cx() == x.cx());
			mat t;
			t.set(NULL, x.cx(), cy());
			if (!I::div(t.data(), data(), cx(), cy(), x.data(), x.cx()))
				t.set(NULL, 0, 0);
			return t;
		}
		
		mat operator^(const int& a) const
		{
			switch(a) {
				case 1:
					return transpose();
				case -1:
					return inv();
				default:
					return mat();
			}
		}
		
		mat transpose() const
		{
			mat t;
			t.set(data(), cx(), cy(), true);
			return t;
		}

		mat cholesky() const
		{
			//assert(cx() == cy());
			mat t;
			if(cx() != cy())
				return t;
			t.set(NULL, cx(), cy(), false, false);
			if (!I::cholesky(t.data(), cx(), data()))
				t.set(NULL, 0, 0);
			return t;
		}
		
		mat inv() const
		{
			//x must be square
			//assert(cx() == cy());
			mat t;
			if (cx() != cy())
				return t;
			t.set(NULL, cx(), cx());
			if (!I::inv(t.data(), cx(), data()))
				t.set(NULL, 0, 0);
			return t;
		}
		
		//cholesky invert
		mat cinv() const
		{
			//x must be square
			//assert(cx() == cy());
			mat t;
			if (cx() != cy())
				return t;
			t.set(NULL, cx(), cx());
			if (!I::cinv(t.data(), cx(), data()))
				t.set(NULL, 0, 0);
			return t;
		}
		
		//invert
		mat operator!() const
		{
			return inv();
		}

		inline mat& zero(const int& cx = -1)
		{
			return set(NULL, cx, cx, false, true);
		}

		mat& set(const Tp* d, const unsigned int& cx, const unsigned int& cy, const bool& transpose = false, const bool& zero = false)
		{
			unsigned int sz = mat::cx() * mat::cy();
			m_cx = transpose ? cy : cx;
			m_cy = transpose ? cx : cy;
			if (cx * cy != sz) {
				if (m_data)
					delete[] m_data;
				m_data = NULL;
				if (cx == 0 || cy == 0)
					return *this;
				m_data = new Tp[cx * cy];
			}
			if (d) {
				if (transpose)
					I::trans(data(), cx, cy, d);
				else
					I::cpy(data(), cx, cy, d);
			} else {
				if(zero)
					I::zero(data(), cx, cy);
			}
			return *this;
		}

		Tp dot(const mat& x) const
		{
			assert(cx() == x.cx());
			assert(cy() == 1 && x.cy() == 1);
			Tp s = 0;
			for (unsigned int i = 0; i < cx(); i++)
				s += m_data[i] * x.m_data[i];
			return s;
		}

		Tp len() const
		{
			assert(cy() == 1);
			return I::sqrt(dot(*this));
		}

		mat normal() const
		{
			mat t(*this);
			t *= 1. / t.len();
			return t;
		}

		mat& normalize()
		{
			return *this *= 1. / len();
		}
	private:
		unsigned int m_cx, m_cy;
		Tp* m_data;
	};
}

//1-m
template<typename Tp, typename U>
matrix::mat<Tp, U>& operator-(const int& x, const matrix::mat<Tp, U>& m)
{
	matrix::mat<Tp, U> t;
	t.id() -= m;
	return t;
}

//x*m
template<typename Tp, typename U>
matrix::mat<Tp, U> operator*(const Tp& x, const matrix::mat<Tp, U>& m)
{
	return m * x;
}

#endif
