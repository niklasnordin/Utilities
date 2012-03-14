#include "Vector.H"

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

Vector::~Vector()
{}

double Vector::operator&(const Vector& v) const
{
    return x_*v.x() +y_*v.y() + z_*v.z();
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
