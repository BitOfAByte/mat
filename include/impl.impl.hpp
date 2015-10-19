#ifndef IMPL_IMPL_H
#define IMPL_IMPL_H

template<typename Tp, typename U>
Tp* matrix::impl<Tp, U>::mul(Tp* dest, const Tp* A, unsigned int ci, unsigned int cy, const Tp* B, unsigned int cx, unsigned int BM)
{
	unsigned int x, y, i;
	if (ci != BM)
		return NULL;
	if (!A || !B)
		return NULL;
	for (i = 0; i < cx*cy; i++)
		dest[i] = 0;
	for (y = 0; y < cy; y++)
		for (x = 0; x < cx; x++)
			for (i = 0; i < ci; i++)
				dest[y * cx + x] += A[y * ci + i] * B[i * ci + x];
	return dest;
}

template<typename Tp, typename U>
Tp* matrix::impl<Tp, U>::tmul(Tp* dest, const Tp* A, unsigned int cy, unsigned int ci, const Tp* B, unsigned int cx, unsigned int BM)
{
	unsigned int x, y, i;
	if (ci != BM)
		return NULL;
	if (!A || !B)
		return NULL;
	for (i = 0; i < cx*cy; i++)
		dest[i] = 0;
	for (y = 0; y < cy; y++)
		for (x = 0; x < cx; x++)
			for (i = 0; i < ci; i++)
				dest[y * cx + x] += A[i * ci + y] * B[i * ci + x];
	return dest;
}

template<typename Tp, typename U>
Tp* matrix::impl<Tp, U>::mult(Tp* dest, const Tp* A, unsigned int ci, unsigned int cy, const Tp* B, unsigned int BN, unsigned int cx)
{
	unsigned int x, y, i;
	if (ci != BN)
		return NULL;
	if (!A || !B)
		return NULL;
	for (i = 0; i < cx*cy; i++)
		dest[i] = 0;
	for (y = 0; y < cy; y++)
		for (x = 0; x < cx; x++)
			for (i = 0; i < ci; i++)
				dest[y * cx + x] += A[y * ci + i] * B[x * ci + i];
	return dest;
}

template<typename Tp, typename U>
Tp* matrix::impl<Tp, U>::mulT(Tp* dest, const Tp* A, unsigned int ci, unsigned int cy, const Tp* B, unsigned int cx, unsigned int BM)
{
	unsigned int x, y, i;
	if (ci != BM)
		return NULL;
	if (!A || !B)
		return NULL;
	for (i = 0; i < cx*cy; i++)
		dest[i] = 0;
	for (x = 0; x < cx; x++)
		for (y = 0; y < cy; y++)
			for (i = 0; i < ci; i++)
				dest[x * cx + y] += A[y * ci + i] * B[i * ci + x];
	return dest;
}

template<typename Tp, typename U>
Tp* matrix::impl<Tp, U>::tmulT(Tp* dest, const Tp* A, unsigned int cy, unsigned int ci, const Tp* B, unsigned int cx, unsigned int BM)
{
	unsigned int x, y, i;
	if (ci != BM)
		return NULL;
	if (!A || !B)
		return NULL;
	for (i = 0; i < cx*cy; i++)
		dest[i] = 0;
	for (x = 0; x < cx; x++)
		for (y = 0; y < cy; y++)
			for (i = 0; i < ci; i++)
				dest[x * cx + y] += A[i * ci + y] * B[i * ci + x];
	return dest;
}

template<typename Tp, typename U>
Tp* matrix::impl<Tp, U>::multT(Tp* dest, const Tp* A, unsigned int ci, unsigned int cy, const Tp* B, unsigned int BN, unsigned int cx)
{
	unsigned int x, y, i;
	if (ci != BN)
		return NULL;
	if (!A || !B)
		return NULL;
	for (i = 0; i < cx*cy; i++)
		dest[i] = 0;
	for (x = 0; x < cx; x++)
		for (y = 0; y < cy; y++)
			for (i = 0; i < ci; i++)
				dest[x * cx + y] += A[y * ci + i] * B[x * ci + i];
	return dest;
}

template<typename Tp, typename U>
Tp* matrix::impl<Tp, U>::div(Tp* dest, const Tp* A, unsigned int AN, unsigned int AM, const Tp* B, unsigned int BN)
{
	Tp* tmp = (Tp*)MAT_ALLOCA(sizeof(Tp) * AM * BN);
	return mul(dest, A, AN, AM, matrix::impl<Tp, U>::inv(tmp, BN, B), BN, BN);
}

template<typename Tp, typename U>
Tp* matrix::impl<Tp, U>::smul(Tp* dest, unsigned int N, unsigned int M, const Tp* src, const Tp& C)
{
	unsigned int x, y;
	if (!src)
		return NULL;
	for (y = 0; y < M; y++)
		for (x = 0; x < N; x++)
			dest[y * N + x] = src[y * N + x] * C;
	return dest;
}

template<typename Tp, typename U>
Tp* matrix::impl<Tp, U>::add(Tp* dest, unsigned int N, unsigned int M, const Tp* A, const Tp* B)
{
	unsigned int i;
	if (!A || !B)
		return NULL;
	for (i = 0; i < M*N; i++)
		dest[i] = A[i] + B[i];
	return dest;
}

template<typename Tp, typename U>
Tp* matrix::impl<Tp, U>::sub(Tp* dest, unsigned int N, unsigned int M, const Tp* A, const Tp* B)
{
	unsigned int i;
	if (!A || !B)
		return NULL;
	for (i = 0; i < M*N; i++)
		dest[i] = A[i] - B[i];
	return dest;
}

template<typename Tp, typename U>
Tp* matrix::impl<Tp, U>::trans(Tp* dest, unsigned int N, unsigned int M, const Tp* src)
{
	unsigned int x, y;
	if (!src)
		return NULL;
	for (y = 0; y < M; y++)
		for (x = 0; x < N; x++)
			dest[x * M + y] = src[y * N + x];
	return dest;
}

template<typename Tp, typename U>
Tp* matrix::impl<Tp, U>::id(Tp* dest, unsigned int N)
{
	unsigned int x, y;
	for (y = 0; y < N; y++)
		for (x = 0; x < N; x++)
			if (x == y)
				dest[y * N + x] = 1;
			else
				dest[y * N + x] = 0;
	return dest;
}

template<typename Tp, typename U>
Tp* matrix::impl<Tp, U>::zero(Tp* dest, unsigned int N, unsigned int M)
{
	unsigned int i;
	for (i = 0; i < N * M; i++)
		dest[i] = 0;
	return dest;
}

template<typename Tp, typename U>
Tp* matrix::impl<Tp, U>::rot(Tp* dest, Tp theta)
{
	dest[0] = U::cos(theta);
	dest[1] = -U::sin(theta);
	dest[2] = U::sin(theta);
	dest[3] = U::cos(theta);
	return dest;
}

template<typename Tp, typename U>
Tp* matrix::impl<Tp, U>::cpy(Tp* dest, unsigned int N, unsigned int M, const Tp* src)
{
	unsigned int i;
	if (!src)
		return NULL;
	for (i = 0; i < N * M; i++)
		dest[i] = src[i];
	return dest;
}

template<typename Tp, typename U>
Tp* matrix::impl<Tp, U>::inv(Tp* dest, unsigned int N, const Tp* src)
{
	Tp* tmp = (Tp*)MAT_ALLOCA(sizeof(Tp) * N * N);
	if (!src)
		return NULL;
	id(dest, N);
	cpy(tmp, N, N, src);
	if (!LUImpl(tmp, N, N, dest, N, N))
		return NULL;
	return dest;
}

template<typename Tp, typename U>
Tp* matrix::impl<Tp, U>::cinv(Tp* dest, unsigned int N, const Tp* src)
{
	Tp* tmp = (Tp*)MAT_ALLOCA(sizeof(Tp) * N * N);
	if (!src)
		return NULL;
	id(dest, N);
	cpy(tmp, N, N, src);
	if (!CholImpl(tmp, N, N, dest, N, N))
		return NULL;
	return dest;
}

template<typename Tp, typename U>
Tp* matrix::impl<Tp, U>::cholesky(Tp* dest, unsigned int N, const Tp* src)
{
	unsigned int x, y, k;
	if (!src)
		return NULL;
	for (y = 0; y < N; y++) {
		for (x = 0; x < y + 1; x++) {
			Tp s = 0;
			for (k = 0; k < x; k++)
				s += dest[y * N + k] * dest[x * N + k];
			if (y == x)
				dest[y * N + x] = U::sqrt(src[y * N + y] - s);
			else
				dest[y * N + x] = 1. / dest[x * N + x] * (src[y * N + x] - s);
			if(dest[y * N + x] != dest[y * N + x])
				return NULL;
		}
	}
	return dest;
}

template<typename Tp, typename U>
Tp* matrix::impl<Tp, U>::det(const Tp* src, unsigned int N)
{
	unsigned int i;
	Tp det;
	Tp* tmp = (Tp*)MAT_ALLOCA(sizeof(Tp) * N * N);
	cpy(tmp, N, N, src);
	det = LUImpl(tmp, N, N, (Tp*)NULL, 0, 0);
	for (i = 0; i < N; i++)
		det *= tmp[i * N + i];
	if (det)
		det = 1. / det;
	return det;
}

#endif
