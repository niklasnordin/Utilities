#include "Matrix.H"

// null constructor
Matrix::Matrix()
: 
    x_(Vector(0,0,0)),
    y_(Vector(0,0,0)),
    z_(Vector(0,0,0))
{}

// construct from components
Matrix::Matrix
(
    const Vector& x, 
    const Vector& y, 
    const Vector& z
)
:
    x_(x), y_(y), z_(z)
{}

// destructor
Matrix::~Matrix()
{}

double Matrix::det() const
{
    double a = y_.y()*z_.z() - y_.z()*z_.y();
    double b = -(y_.x()*z_.z() - y_.z()*z_.x());
    double c = y_.x()*z_.y() - y_.y()*z_.x();
    return x_.x()*a + x_.y()*b + x_.z()*c;;
}

    // transpose
Matrix Matrix::transpose() const
{
    Vector a(x_.x(), y_.x(), z_.x());
    Vector b(x_.y(), y_.y(), z_.y());
    Vector c(x_.z(), y_.z(), z_.z());
    return Matrix(a, b, c);
}

// inverse
Matrix Matrix::inverse() const
{
    double det = this->det();
    double idet = 1.0/det;
    
    double ax = idet*(y_.y()*z_.z() - z_.y()*y_.z());
    double ay = -idet*(y_.x()*z_.z() - z_.x()*y_.z());
    double az = idet*(y_.x()*z_.y() - z_.x()*y_.y());
    Vector a(ax, ay, az);
    
    double bx = -idet*(x_.y()*z_.z() - z_.y()*x_.z());
    double by = idet*(x_.x()*z_.z() - z_.x()*x_.z());
    double bz = -idet*(x_.x()*z_.y() - z_.x()*x_.y());
    Vector b(bx, by, bz);
    
    double cx = idet*(x_.y()*y_.z() - y_.y()*x_.z());
    double cy = -idet*(x_.x()*y_.z() - y_.x()*x_.z());
    double cz = idet*(x_.x()*y_.y() - y_.x()*x_.y());
    Vector c(cx, cy, cz);
    
    Matrix m(a,b,c);
  
    return m.transpose();
}

// operators
const Vector Matrix::operator&(const Vector& b) const
{
    double x = x_ & b;
    double y = y_ & b;
    double z = z_ & b;
    return Vector(x, y, z);
}

const Matrix Matrix::operator*(const Matrix& m) const
{
    Matrix bT = m.transpose();
    Vector a(x_ & bT.x(), x_ & bT.y(), x_ & bT.z());
    Vector b(y_ & bT.x(), y_ & bT.y(), y_ & bT.z());
    Vector c(z_ & bT.x(), z_ & bT.y(), z_ & bT.z());
    return Matrix(a,b,c);
}

const Matrix Matrix::operator+(const Matrix& m) const
{
    Vector a = x_ + m.x();
    Vector b = y_ + m.y();
    Vector c = z_ + m.z();
    return Matrix(a,b,c);
}

const Matrix Matrix::operator-(const Matrix& m) const
{
    Vector a = x_ - m.x();
    Vector b = y_ - m.y();
    Vector c = z_ - m.z();
    return Matrix(a,b,c);
}

