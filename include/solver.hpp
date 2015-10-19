/*M///////////////////////////////////////////////////////////////////////////////////////
//
//  IMPORTANT: READ BEFORE DOWNLOADING, COPYING, INSTALLING OR USING.
//
//  By downloading, copying, installing or using the software you agree to this license.
//  If you do not agree to this license, do not download, install,
//  copy or use the software.
//
//
//			     License Agreement
//		  For Open Source Computer Vision Library
//
// Copyright (C) 2000-2008, Intel Corporation, all rights reserved.
// Copyright (C) 2009, Willow Garage Inc., all rights reserved.
// Third party copyrights are property of their respective owners.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
//   * Redistribution's of source code must retain the above copyright notice,
//     this list of conditions and the following disclaimer.
//
//   * Redistribution's in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//
//   * The name of the copyright holders may not be used to endorse or promote products
//     derived from this software without specific prior written permission.
//
// This software is provided by the copyright holders and contributors "as is" and
// any express or implied warranties, including, but not limited to, the implied
// warranties of merchantability and fitness for a particular purpose are disclaimed.
// In no event shall the Intel Corporation or contributors be liable for any direct,
// indirect, incidental, special, exemplary, or consequential damages
// (including, but not limited to, procurement of substitute goods or services;
// loss of use, data, or profits; or business interruption) however caused
// and on any theory of liability, whether in contract, strict liability,
// or tort (including negligence or otherwise) arising in any way out of
// the use of this software, even if advised of the possibility of such damage.
//
//M*/

/*LUImpl was taken from OpenCV lapack.cpp and is available under the above
  license*/

#ifndef SOLVER_H
#define SOLVER_H

template<typename Tp, typename M>
int matrix::impl<Tp, M>::LUImpl(Tp* A, unsigned int astep, int m, Tp* b, unsigned int bstep, int n)
{
	int i, j, k, p = 1;

	for (i = 0; i < m; i++) {
		k = i;

		for (j = i + 1; j < m; j++ )
			if (M::abs(A[j * astep + i]) > M::abs(A[k * astep + i]) )
				k = j;

		if (M::abs(A[k * astep + i]) < M::epsilon(A[0]))
			return 0;

		if (k != i) {
			for (j = i; j < m; j++)
				M::swap(A[i * astep + j], A[k * astep + j]);
			if (b)
				for (j = 0; j < n; j++)
					M::swap(b[i * bstep + j], b[k*bstep + j]);
			p = -p;
		}

		Tp d = -1 / A[i * astep + i];

		for (j = i + 1; j < m; j++) {
			Tp alpha = A[j * astep + i] * d;

			for (k = i + 1; k < m; k++)
				A[j * astep + k] += alpha * A[ i * astep + k];

			if (b)
				for (k = 0; k < n; k++)
					b[j * bstep + k] += alpha * b[i * bstep + k];
		}

		A[i * astep + i] = -d;
	}

	if (b) {
		for (i = m - 1; i >= 0; i--)
			for (j = 0; j < n; j++) {
				Tp s = b[i * bstep + j];
				for (k = i + 1; k < m; k++)
					s -= A[ i *astep + k] * b[k * bstep + j];
				b[i * bstep + j] = s * A[i * astep + i];
			}
	}

	return p;
}

template<typename Tp, typename M> bool
matrix::impl<Tp, M>::CholImpl(Tp* A, unsigned int astep, int m, Tp* b, unsigned int bstep, int n)
{
	Tp* L = A;
	int i, j, k;
	double s;
	
	for (i = 0; i < m; i++) {
		for (j = 0; j < i; j++) {
			s = A[i * astep + j];
			for(k = 0; k < j; k++)
				s -= L[i * astep + k] * L[j * astep + k];
			L[i * astep + j] = (Tp)(s * L[j * astep + j]);
		}
		s = A[i * astep + i];
		for (k = 0; k < j; k++) {
			double t = L[i * astep + k];
			s -= t * t;
		}
		if (s < M::epsilon(s))
			return false;
		L[i * astep + i] = (Tp)(1. / M::sqrt(s));
	}
	
	if (!b)
		return true;

	// LLt x = b
	// 1: L y = b
	// 2. Lt x = y
	
	/*
	[ L00		  ]  y0   b0
	[ L10 L11	  ]  y1 = b1
	[ L20 L21 L22	  ]  y2   b2
	[ L30 L31 L32 L33 ]  y3   b3
	
	[ L00 L10 L20 L30 ]  x0   y0
	[     L11 L21 L31 ]  x1 = y1
	[	  L22 L32 ]  x2   y2
	[	      L33 ]  x3   y3
	*/
	
	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++) {
			s = b[i * bstep + j];
			for(k = 0; k < i; k++)
				s -= L[i * astep + k] * b[k * bstep + j];
			b[i * bstep + j] = (Tp)(s * L[i * astep + i]);
		}
	}
	
	for (i = m - 1; i >= 0; i--) {
		for (j = 0; j < n; j++) {
			s = b[i * bstep + j];
			for (k = m - 1; k > i; k--)
				s -= L[k * astep + i] * b[k * bstep + j];
			b[i * bstep + j] = (Tp)(s * L[i * astep + i]);
		}
	}
	
	return true;
}

#endif
