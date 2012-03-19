#include <iostream>
#include <fstream>

using namespace std;

#ifndef VECTOR_H
#define VECTOR_H

#define SMALL 1.0e-15
typedef double scalar;
typedef int label;

class Vector
{
private:
    double x_;
    double y_;
    double z_;

public:

    //constructors
    // null-constructor
    Vector();

    // construct from components
    Vector
    (
        const scalar& x,
	const scalar& y,
	const scalar& z
    );

    // construct from stream
    Vector(ifstream& is);

    // destructor
    ~Vector();

    // access functions
    inline const double& x() const { return x_; }
    inline const double& y() const { return y_; }
    inline const double& z() const { return z_; }

    //operators

    Vector operator+(const Vector& v) const;
    Vector operator-(const Vector& v) const;
    Vector operator*(const scalar& a) const;
    Vector operator/(const scalar& a) const;

    // inner product
    double operator&(const Vector& v) const;

    // cross product
    const Vector operator^(const Vector& b) const;

    double mag() const;

    friend const double mag(const Vector& v)
    {
        return v.mag();
    }

    // left hand multiplication with double
    friend const Vector operator*(const double& d, const Vector& v)
    {
        return Vector(d*v.x(), d*v.y(), d*v.z());
    }

    friend ostream& operator<<(ostream& os, const Vector& v)
    {
        os << v.x() << " " << v.y() << " " << v.z();
        return os;
    }
};

#endif
