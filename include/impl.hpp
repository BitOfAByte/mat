#ifndef MAT_UTIL_H
#define MAT_UTIL_H

#include <math.hpp>

///N = cx, M = cy

///Note: to multiply a matrix and a vector, transpose the vector (use mult)

namespace matrix
{
	template<typename Tp, typename U = mathutils>
	class impl
	{
	public:
		///results in BNxAM
		///find the inner product of two matricies (equal to tmultT)
		static Tp* mul(Tp* dest, const Tp* A, unsigned int AN, unsigned int AM, const Tp* B, unsigned int BN, unsigned int BM);
		///results in BNxAN
		///find the inner product of two matricies, with A transposed
		static Tp* tmul(Tp* dest, const Tp* A, unsigned int AN, unsigned int AM, const Tp* B, unsigned int BN, unsigned int BM);
		///results in BMxAM
		///find the inner product of two matricies, with B transposed
		static Tp* mult(Tp* dest, const Tp* A, unsigned int AN, unsigned int AM, const Tp* B, unsigned int BN, unsigned int BM);
		///results in AMxBN
		///find the inner product of two matricies, with the result transposed (equal to tmult())
		static Tp* mulT(Tp* dest, const Tp* A, unsigned int AN, unsigned int AM, const Tp* B, unsigned int BN, unsigned int BM);
		///results in ANxBN
		///find the inner product of two matricies, with the result and A transposed
		static Tp* tmulT(Tp* dest, const Tp* A, unsigned int AN, unsigned int AM, const Tp* B, unsigned int BN, unsigned int BM);
		///results in AMxBM
		///find the inner product of two matricies, with the result and B transposed
		static Tp* multT(Tp* dest, const Tp* A, unsigned int AN, unsigned int AM, const Tp* B, unsigned int BN, unsigned int BM);
		
		///results BNxAM
		///multiply (A)(B^-1)
		static Tp* div(Tp* dest, const Tp* A, unsigned int AN, unsigned int AM, const Tp* B, unsigned int BN);		

		///Multiply matrix by a scalar
		static Tp* smul(Tp* dest, unsigned int N, unsigned int M, const Tp* src, const Tp& C);
		///component-wise add two matricies 
		static Tp* add(Tp* dest, unsigned int N, unsigned int M, const Tp* A, const Tp* B);
		///component-wise subtract two matricies 
		static Tp* sub(Tp* dest, unsigned int N, unsigned int M, const Tp* A, const Tp* B);

		///transpose a matrix
		static Tp* trans(Tp* dest, unsigned int N, unsigned int M, const Tp* src);
		///set a matrix to identity
		static Tp* id(Tp* dest, unsigned int N);
		///set a matrix to zero
		static Tp* zero(Tp* dest, unsigned int N, unsigned int M);
		///creates a rotational matrix in dest
		static Tp* rot(Tp* dest, Tp theta);
		///copy a matrix
		static Tp* cpy(Tp* dest, unsigned int N, unsigned int M, const Tp* src);

		///find inverse matrix using lu decomp.
		static Tp* inv(Tp* dest, unsigned int N, const Tp* src);
		///find inverse matrix by cholesky decomp.
		static Tp* cinv(Tp* dest, unsigned int N, const Tp* src);
		///preform cholskey decomposition on matrix src
		///(dest)*(dest^T) = src
		static Tp* cholesky(Tp*dest, unsigned int N, const Tp* src);

		///find the determinant of a matrix
		static Tp* det(const Tp* src, unsigned int N);

	private:
		static int LUImpl(Tp* A, unsigned int astep, int m, Tp* b, unsigned int bstep, int n);
		static bool CholImpl(Tp* A, unsigned int astep, int m, Tp* b, unsigned int bstep, int n);
	};
}

#include <impl.impl.hpp>
#include <solver.hpp>

#endif
