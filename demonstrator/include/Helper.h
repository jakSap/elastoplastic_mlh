//
// Created by Johannes Martin on 26.09.22.
//

#ifndef MESHLESSHYDRO_HELPER_H
#define MESHLESSHYDRO_HELPER_H

#include <cmath>

#include "parameter.h"
#include "Logger.h"

/// Matrix inversion taken from: https://stackoverflow.com/a/3525136/6208997

extern "C" {
    // LU decomoposition of a general matrix
    void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);

    // generate inverse of a matrix given its LU decomposition
    void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);

    // eigenvalues and eigenvectors of a real symmetric matrix
    void dsyev_(char* jobz, char* uplo, int* n, double* a, int* lda,
                double* w, double* work, int* lwork, int* info);
}

class Helper {

public:
    void inverseMatrix(double *A, int N);
    static double dotProduct(double *a, double *b);
    static void crossProduct(double *a, double *b, double *crossProduct);


    /**
     * @brief Matrix-matrix multiplication C = A * B for DIM x DIM matrices
     *
     * @param[in] A left matrix, row-major [DIM*DIM]
     * @param[in] B right matrix, row-major [DIM*DIM]
     * @param[out] C result matrix, row-major [DIM*DIM], must be pre-allocated
     *             C may NOT alias A or B
     */
    static void matMul(double *A, double *B, double *C);

    /**
     * @brief Largest eigenvalue of a symmetric DIM x DIM matrix.
     *
     * 2D uses the closed-form quadratic; 3D uses the analytic trigonometric
     * method for symmetric 3x3 matrices (Smith 1961). Used by the Grady-Kipp
     * damage model to find the maximum tensile principal stress.
     *
     * @param[in] S symmetric matrix, row-major [DIM*DIM]
     * @return the maximum eigenvalue
     */
    static double maxEigenvalueSym(const double *S);

    /**
     * @brief Eigenvalues and eigenvectors of a symmetric DIM x DIM matrix (LAPACK dsyev).
     *
     * Used by the GIZMO elastic stress flux and the tensile-principal damping, which need
     * the full eigenbasis (not just the largest eigenvalue) to project/damp the deviatoric
     * stress in 3D, where there is no closed-form eigenbasis.
     *
     * @param[in]  S    symmetric matrix, row-major [DIM*DIM]
     * @param[out] eval eigenvalues, ascending [DIM]
     * @param[out] evec eigenvectors stored column-wise: eigenvector k is evec[k*DIM + i],
     *                  i = component (LAPACK's column-major output) [DIM*DIM]
     */
    static void eigenDecompositionSym(const double *S, double *eval, double *evec);

    /**
     * good resource for 3D implementation: https://math.stackexchange.com/a/897677
     *
     * @param[in] a normed vector to be aligned with b
     * @param[in] b normed vector to which a shall be rotated by the matrix
     * @param[out] Lambda rotation matrix, must be pre-allocated [DIM*DIM]
     *             indexed lambda_ij = Lambda[j+DIM*i], DIM=2
     */
    static void rotationMatrix2D(double *a, double *b, double *Lambda);
#if DIM==3
    static void rotationMatrix3D(double *a, double *b, double *Lambda);
#endif


private:
    double WORK[DIM*DIM];
    int IPIV[DIM];

};


#endif //MESHLESSHYDRO_HELPER_H
