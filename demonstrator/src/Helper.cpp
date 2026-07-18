//
// Created by Johannes Martin on 26.09.22.
//

#include "../include/Helper.h"



void Helper::inverseMatrix(double *A, int N){
    //int *IPIV = new int[N];
    int LWORK = N*N;
    //double *WORK = new double[LWORK];
    int INFO;

    dgetrf_(&N,&N,A,&N,IPIV,&INFO);
    dgetri_(&N,A,&N,IPIV,WORK,&LWORK,&INFO);

    //delete[] IPIV;
    //delete[] WORK;
}

bool Helper::inverseMatrixChecked(double *A, int N){
    int LWORK = N*N;
    int INFO;

    // Frobenius norm of A captured before dgetrf overwrites it with its LU
    // factors; needed for the condition-number test below.
    double normA = 0.;
    for (int k=0; k<N*N; ++k) normA += A[k]*A[k];

    dgetrf_(&N,&N,A,&N,IPIV,&INFO);
    if (INFO != 0) return false;      // singular: U has a zero pivot
    dgetri_(&N,A,&N,IPIV,WORK,&LWORK,&INFO);
    if (INFO != 0) return false;

    // Reject ill-conditioned (not merely exactly-singular) matrices. The
    // augmented order-2/3 gradient moment matrix goes near-singular on the
    // one-sided stencils of free-surface particles; for order 3 the linear
    // (gradient) and cubic basis terms share odd parity and become collinear
    // there, so the odd block collapses while LAPACK still returns a garbage
    // inverse that blows the recovered gradient up. kappa_F = ||A||_F *
    // ||A^-1||_F is a cheap Frobenius proxy for the condition number; above
    // GRAD_MAT_COND_MAX we report failure so the caller falls back to the
    // well-conditioned first-order weights instead.
    double normInv = 0.;
    for (int k=0; k<N*N; ++k) normInv += A[k]*A[k];
    return std::sqrt(normA*normInv) <= GRAD_MAT_COND_MAX;
}

double Helper::dotProduct(double *a, double *b){
    double res = 0.;
    for (int alpha=0; alpha<DIM; ++alpha){
        res += a[alpha] * b[alpha];
    }
    return res;
}

void Helper::crossProduct(double *a, double *b, double *crossProduct){
#if DIM == 3
    crossProduct[0] = a[1]*b[2] - a[2]*b[1];
    crossProduct[1] = a[2]*b[0] - a[0]*b[2];
    crossProduct[2] = a[0]*b[1] - a[1]*b[0];
#else
    Logger(ERROR) << "Cross product not defined for 2D. - Aborting.";
    exit(9);
#endif
}

void Helper::matMul(double *A, double *B, double *C){
    for (int i=0; i<DIM; ++i){
        for (int j=0; j<DIM; ++j){
            C[j + DIM*i] = 0.;
            for (int k=0; k<DIM; ++k){
                C[j + DIM*i] += A[k + DIM*i] * B[j + DIM*k];
            }
        }
    }
}

double Helper::maxEigenvalueSym(const double *S){
#if DIM == 2
    // S = [[a, b], [b, c]]; eigenvalues = mean +- sqrt(((a-c)/2)^2 + b^2).
    const double a = S[0], b = S[1], c = S[3];
    const double mean = 0.5 * (a + c);
    const double diff = 0.5 * (a - c);
    return mean + std::sqrt(diff*diff + b*b);
#else // DIM == 3
    // Analytic eigenvalues of a symmetric 3x3 matrix (Smith 1961).
    const double a = S[0], b = S[4], c = S[8];   // diagonal
    const double d = S[1], e = S[5], f = S[2];   // S01, S12, S02
    const double p1 = d*d + e*e + f*f;
    if (p1 <= 0.0) {                              // already diagonal
        return std::max(a, std::max(b, c));
    }
    const double q = (a + b + c) / 3.0;
    const double aa = a - q, bb = b - q, cc = c - q;
    const double p2 = aa*aa + bb*bb + cc*cc + 2.0*p1;
    const double p = std::sqrt(p2 / 6.0);
    // B = (1/p) (S - q I); r = det(B)/2
    const double b00 = aa/p, b11 = bb/p, b22 = cc/p;
    const double b01 = d/p, b12 = e/p, b02 = f/p;
    double r = 0.5 * (b00*(b11*b22 - b12*b12)
                      - b01*(b01*b22 - b12*b02)
                      + b02*(b01*b12 - b11*b02));
    if (r < -1.0) r = -1.0; else if (r > 1.0) r = 1.0;
    const double phi = std::acos(r) / 3.0;
    // largest eigenvalue corresponds to the smallest acos argument
    return q + 2.0 * p * std::cos(phi);
#endif
}

void Helper::eigenDecompositionSym(const double *S, double *eval, double *evec){
    // A symmetric matrix has identical row-major and column-major storage, so S can be
    // handed to LAPACK directly. On exit 'a' (column-major) holds the orthonormal
    // eigenvectors in its columns: eigenvector k = a[k*DIM + i].
    int n = DIM, lda = DIM, info;
    int lwork = 3 * DIM;                 // >= max(1, 3n-1) for all DIM >= 1
    double a[DIM * DIM];
    double work[3 * DIM];
    for (int k = 0; k < DIM * DIM; ++k) a[k] = S[k];
    char jobz = 'V', uplo = 'U';
    dsyev_(&jobz, &uplo, &n, a, &lda, eval, work, &lwork, &info);
    for (int k = 0; k < DIM * DIM; ++k) evec[k] = a[k];
}

void Helper::rotationMatrix2D(double *a, double *b, double *Lambda){
    Lambda[0] = a[0]*b[0] + a[1]*b[1];
    Lambda[1] = -(a[0]*b[1] - a[1]*b[0]);
    //Lambda[2] = a[0]*b[1] - a[1]*b[0];
    Lambda[2] = -Lambda[1];
    Lambda[3] = Lambda[0];
}

#if DIM==3
void Helper::rotationMatrix3D(double *a, double *b, double *Lambda){
    double v[DIM];
    crossProduct(a, b, v);
    double cosAB = dotProduct(a, b); // a and b MUST be normed

    // TODO: cosAB == -1 fails !!!

    double n = 1./(1. + cosAB);

    Lambda[0] = 1.-n*(v[2]*v[2]+v[1]*v[1]);
    Lambda[1] = -v[2]+n*v[0]*v[1];
    Lambda[2] = v[1]+n*v[0]*v[2];
    Lambda[3] = v[2]+n*v[0]*v[1];
    Lambda[4] = 1.-n*(v[2]*v[2]+v[0]*v[0]);
    Lambda[5] = -v[0]+n*v[1]*v[2];
    Lambda[6] = -v[1]+n*v[0]*v[2];
    Lambda[7] = v[0]+n*v[1]*v[2];
    Lambda[8] = 1.-n*(v[1]*v[1]+v[0]*v[0]);

    // check for aligned vectors
    //double cosAB = ab/(sqrt(dotProduct(a, a))*sqrt(dotProduct(b, b)));
    //if (abs(cosAB) < 1. + ROT_3D_ALIGN_TOL && abs(cosAB) > 1. - ROT_3D_ALIGN_TOL){
    //    Logger(WARN) << "Very small angle between a and b. - Check Lambda!!";
    //    Logger(DEBUG) << "a = [" << a[0] << ", " << a[1] << ", " << a[2] << "]";
    //    Logger(DEBUG) << "b = [" << b[0] << ", " << b[1] << ", " << b[2] << "]";
    //    Logger(DEBUG) << "Lambda = [" << Lambda[0] << ", " << Lambda[1] << ", " << Lambda[2];
    //    Logger(DEBUG) << "          " << Lambda[3] << ", " << Lambda[4] << ", " << Lambda[5];
    //    Logger(DEBUG) << "          " << Lambda[6] << ", " << Lambda[7] << ", " << Lambda[8] << "]";
    //}
}
#endif
