#include <iostream>

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

    //constructor
    Vector();
    Vector
    (
        const scalar& x,
	const scalar& y,
	const scalar& z
    );

    //destructor
    ~Vector();

    //access functions
    inline const scalar& x() const 
    {
        return x_;
    }

    inline const scalar& y() const
    {
        return y_;
    }

    inline const scalar& z() const
    {
        return z_;
    }

    //operators

    double operator&(const Vector& v) const;
    Vector operator+(const Vector& v) const;
    Vector operator-(const Vector& v) const;
    Vector operator*(const scalar& a) const;

    friend Vector operator*(const scalar& a, const Vector& v)
    {
        return Vector(a*v.x(), a*v.y(), a*v.z());
    }

    friend ostream& operator<<(ostream& os, const Vector& v)
    {
        os << v.x() << " " << v.y() << " " << v.z();
        return os;
    }
};

#endif
