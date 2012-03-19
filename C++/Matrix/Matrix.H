#include "Vector.H"

using namespace std;

#ifndef MATRIX_H
#define MATRIX_H


class Matrix
{
  // row vectors of the matrix
    Vector x_, y_, z_;

public:

    // null constructor
    Matrix();

    // construct from components
    Matrix(const Vector& x, const Vector& y, const Vector& z);

    // destructor
    ~Matrix();

    // acess functions
    inline const Vector& x() const { return x_; }
    inline const Vector& y() const { return y_; }
    inline const Vector& z() const { return z_; }

    // determinant
    double det() const;

    // transpose
    Matrix transpose() const;

    // inverse
    Matrix inverse() const;

    // operators
    const Vector operator&(const Vector& b) const;
    const Matrix operator*(const Matrix& m) const;
    const Matrix operator+(const Matrix& m) const;
    const Matrix operator-(const Matrix& m) const;

    friend ostream& operator<<(ostream& os, const Matrix& v)
    {
        os << v.x() << endl << v.y() << endl << v.z();
        return os;
    }

};

#endif
