#ifndef LUDECOMPOSITION_H
#define LUDECOMPOSITION_H

// Author: 	GhasemHeydari
// Email: 	ghasem.heydari2@gmail.com

class Matrix;

class LUDecomposition {
private:
    double** LU;
    int m, n, pivsign;
    int* piv;
public:
    LUDecomposition(const Matrix &);
    ~LUDecomposition(void);
    bool isNonsingular();
    double det();
    double* getdoublePivot();
    Matrix getL();
    Matrix getL(Matrix X);
    Matrix getU();
    Matrix solve(Matrix & B);
    Matrix inv();
};

#endif /* LUDECOMPOSITION_H */
