#include "Vector.H"
#include <math.h>

Vector::Vector()
:
    x_(0),
    y_(0),
    z_(0)
{}

Vector::Vector
(
    const scalar& x,
    const scalar& y,
    const scalar& z
)
:
    x_(x), y_(y), z_(z)
{}

// construct from stream
Vector::Vector(ifstream& is)
{
    is >> x_;
    is >> y_;
    is >> z_;
}

Vector::~Vector()
{}

double Vector::operator&(const Vector& v) const
{
    return x_*v.x() + y_*v.y() + z_*v.z();
}

Vector Vector::operator+(const Vector& v) const
{
    return Vector(x_+v.x(), y_+v.y(), z_+v.z());
}

Vector Vector::operator-(const Vector& v) const
{
    return Vector(x_-v.x(), y_-v.y(), z_-v.z());
}

Vector Vector::operator*(const scalar& a) const
{
    return Vector(a*x_, a*y_, a*z_);
}

Vector Vector::operator/(const double& d) const
{
    return Vector(x_/d, y_/d, z_/d);
}

double Vector::mag() const
{
  return sqrt( *this & *this );
}

// cross product
const Vector Vector::operator^(const Vector& b) const
{
    double x = y_*b.z() - z_*b.y();
    double y = -(x_*b.z() - z_*b.x());
    double z = x_*b.y() - y_*b.x();
  
    return Vector(x, y, z);
}

