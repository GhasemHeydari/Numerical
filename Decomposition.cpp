#include <math.h>

// Author: 	GhasemHeydari
// Email: 	ghasem.heydari2@gmail.com

#include "LUD.h"
#include "MATLAB.h"
#include "Matrix.h"
#include "Vector.h"

LUDecomposition::LUDecomposition(const Matrix & A) :
    m(A.rows()),
    n(A.cols()),
    pivsign(1),
    piv(new int[m])
{
    int i, j, k;
    LU = new double*[m];
    for (i = 0; i < m; i++) {
        LU[i] = new double[n];
        for (int j = 0; j < n; j++) {
            LU[i][j] = A(i, j);
        }
    }

    for (i = 0; i < m; i++) {
        piv[i] = i;
    }
    double* LUrowi;
    double* LUcolj = new double[m];

    for (j = 0; j < n; j++) {
        // Make a copy of the j-th column to localize references.

        for (i = 0; i < m; i++) {
            LUcolj[i] = LU[i][j];
        }

        // Apply previous transformations.

        for (i = 0; i < m; i++) {
            LUrowi = LU[i];

            // Most of the time is spent in the following dot product.

            int kmax = (i < j) ? i : j;
            double s = 0.0;
            for (k = 0; k < kmax; k++) {
                s += LUrowi[k] * LUcolj[k];
            }

            LUrowi[j] = LUcolj[i] -= s;
        }

        // Find pivot and exchange if necessary.

        int p = j;
        for (i = j + 1; i < m; i++) {
            if (fabs(LUcolj[i]) > fabs(LUcolj[p])) {
                p = i;
            }
        }
        if (p != j) {
            for (k = 0; k < n; k++) {
                double t = LU[p][k];
                LU[p][k] = LU[j][k];
                LU[j][k] = t;
            }
            k = piv[p];
            piv[p] = piv[j];
            piv[j] = k;
            pivsign = -pivsign;
        }

        // Compute multipliers.

        if (j < m && LU[j][j] != 0.0) {
            for (i = j + 1; i < m; i++) {
                LU[i][j] /= LU[j][j];
            }
        }
    }

    delete LUcolj;
}

LUDecomposition::~LUDecomposition(void) {
    for (int i = 0; i < m; i++) {
        delete[] LU[i];
    }
    delete[] LU;
    delete[] piv;
}

bool LUDecomposition::isNonsingular() {
    for (int j = 0; j < n; j++) {
        if (LU[j][j] == 0)
            return false;
    }
    return true;
}

Matrix LUDecomposition::getL() {
    Matrix X(m, n);
    int i, j;
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            if (i > j) {
                X(i, j) = LU[i][j];
            }
            else if (i == j) {
                X(i, j) = 1.0;
            }
            else {
                X(i, j) = 0.0;
            }
        }
    }
    return X;
}

Matrix LUDecomposition::getL(Matrix X) {
    //Matrix X(m, n);
    int i, j;
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            if (i > j) {
                X(i, j) = LU[i][j];
            }
            else if (i == j) {
                X(i, j) = 1.0;
            }
            else {
                X(i, j) = 0.0;
            }
        }
    }
    return X;
}

Matrix LUDecomposition::getU() {
    Matrix X(n, n);
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (i <= j) {
                X(i, j) = LU[i][j];
            }
            else {
                X(i, j) = 0.0;
            }
        }
    }
    return X;
}

double* LUDecomposition::getdoublePivot() {
    double* vals = new double[m];
    for (int i = 0; i < m; i++) {
        vals[i] = (double)piv[i];
    }
    return vals;
}

double LUDecomposition::det() {
    assert(m == n);

    double d = (double)pivsign;
    for (int j = 0; j < n; j++) {
        d *= LU[j][j];
    }
    return d;
}

Matrix LUDecomposition::solve(Matrix & B) {
    assert(m == B.rows());
    int i, j, k;

    // Copy right hand side with pivoting
    int nx = B.cols();
    Matrix Xmat(m, nx);
    for (i = 0; i < m; i++) {
        Xmat.setRow(i, B.getRow(piv[i]));
    }

    // Solve L*Y = B(piv,:)
    for (k = 0; k < n; k++) {
        for (i = k + 1; i < n; i++) {
            for (j = 0; j < nx; j++) {
                Xmat(i, j) -= Xmat(k, j) * LU[i][k];
            }
        }
    }
    // Solve U*X = Y;
    for (k = n - 1; k >= 0; k--) {
        for (j = 0; j < nx; j++) {
            Xmat(k, j) /= LU[k][k];
        }
        for (i = 0; i < k; i++) {
            for (j = 0; j < nx; j++) {
                Xmat(i, j) -= Xmat(k, j) * LU[i][k];
            }
        }
    }
    return Xmat;
}

Matrix LUDecomposition::inv() {
    assert(m == n);

    Matrix idn = eye(m, m);
    return solve(idn);
}
