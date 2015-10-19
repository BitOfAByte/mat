#include "mat.hpp"

#include <iostream>

template<typename Tp, typename U>
std::ostream& operator<<(std::ostream& s, const matrix::mat<Tp, U>& m)
{
	size_t x = 0, y = 0;
	for (y = 0; y < m.cy(); y++) {
		if(y)
			s << std::endl;
		s << "[";
		for(x = 0; x < m.cx(); x++) {
			if(x)
				s << ", ";
			s << m.at(x, y);
		}
		s << "]";
	}
	if(x == 0 || y == 0)
		s << "[]";
	return s;
}

int main(int argc, char** argv)
{
	using namespace std;
	using namespace matrix;
	//random test matrix
	float aa[] = {
		0.60, 1.20, 1.40, 1.00, 0.60, 1.00, 1.20, 0.40, 1.80, 0.20,
		0.40, 1.40, 0.00, 1.80, 0.60, 1.20, 0.00, 1.20, 0.40, 1.20,
		0.20, 1.60, 1.40, 1.80, 0.40, 0.00, 0.40, 0.60, 1.40, 1.00,
		1.80, 0.40, 0.40, 1.60, 1.80, 1.40, 0.60, 1.20, 0.20, 0.40,
		1.80, 0.60, 0.20, 1.80, 0.80, 1.40, 1.60, 0.80, 1.00, 0.00,
		0.60, 1.20, 0.20, 0.00, 1.20, 0.60, 0.40, 0.00, 1.20, 0.20,
		1.00, 1.00, 0.80, 1.40, 1.20, 1.00, 1.20, 1.80, 0.60, 1.40,
		0.80, 1.00, 0.40, 1.00, 0.80, 1.40, 0.80, 0.80, 0.60, 0.00,
		1.40, 1.60, 1.20, 1.60, 1.60, 0.80, 0.60, 0.20, 0.80, 1.80,
		0.40, 0.00, 1.20, 1.60, 1.80, 0.40, 1.20, 1.20, 0.80, 1.80,
		};
	mat<> z(aa, 10, 10);
	cout.precision(1);
	cout.setf(ios::fixed, ios::floatfield);
	cout << "z" << endl;
	cout << z << endl;
	cout << "e" << endl;
	cout << z / z << endl;
	cout << "inv" << endl;
	cout << z.inv() << endl;
	cout << "cinv" << endl;
	cout << z.cinv() << endl;
	cout << "c" << endl;
	cout << z.cholesky() << endl;
	//another test matrix
	float a[] = {25, 15, -5,
		     15, 18, 0,
		     -5, 0, 11,
		    };
	mat<> t(a, 3, 3);
	cout << "t" << endl;
	cout << t << endl;
	cout << "t^-1" << endl;
	cout << (t^-1) << endl;
	cout << "c" << endl;
	cout << (t^T) << endl;
	cout << "t*t^-1" << endl;
	cout << (t * (t ^ -1)) << endl;
	cout << "c*c^t" << endl;
	cout << t.cholesky() * (t.cholesky()^T) << std::endl;
	return 0;
}
