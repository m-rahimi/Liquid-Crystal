/* 
 * File:   Eigen.h
 * Author: amin
 *
 * Created on May 13, 2015, 9:08 AM
 */

#ifndef EIGEN_H
#define	EIGEN_H

/* Symmetric matrix A => eigenvectors in columns of V, corresponding
   eigenvalues in d. */
void eigen_decomposition(float A[3][3], float V[3][3], float d[3]);

#endif	/* EIGEN_H */

