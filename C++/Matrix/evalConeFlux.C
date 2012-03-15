#include <iostream>
#include <string.h>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <math.h>

using namespace std;

#define SMALL 1.0e-15

class Vector
{
    double x_, y_, z_;

public:

    // null constructor
    Vector()
    : 
        x_(0.0),
	y_(0.0),
	z_(0.0)
    {}

    // construct from components
    Vector(const double& x, const double& y, const double z)
    :
      x_(x), y_(y), z_(z)
    {}

    // construct from stream
    Vector(ifstream& is)
    {
        is >> x_;
	is >> y_;
	is >> z_;
    }

    // destructor
    ~Vector()
    {}

    // access functions
    inline const double& x() const { return x_; }
    inline const double& y() const { return y_; }
    inline const double& z() const { return z_; }

  
    // operators
    double operator&(const Vector& rhs) const
    {
        return x_*rhs.x() + y_*rhs.y() + z_*rhs.z();
    }

    double mag() const
    {
        return sqrt( *this & *this );
    }

    // cross product
    const Vector operator^(const Vector& b) const
    {
        double x = y_*b.z() - z_*b.y();
	double y = -(x_*b.z() - z_*b.x());
	double z = x_*b.y() - y_*b.x();

        return Vector(x, y, z);
    }

    const Vector operator-(const Vector& rhs) const
    {
        return Vector(x_ - rhs.x(), y_ - rhs.y(), z_ - rhs.z());
    }

    const Vector operator+(const Vector& rhs) const
    {
        return Vector(x_ + rhs.x(), y_ + rhs.y(), z_ + rhs.z());
    }
  
    // right hand multiplication with double
    Vector operator*(const double& d) const
    {
        return Vector(d*x_, d*y_, d*z_);
    }

    Vector operator/(const double& d) const
    {
        return Vector(x_/d, y_/d, z_/d);
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

class Matrix
{
  // row vectors of the matrix
    Vector x_, y_, z_;

public:

    // null constructor
    Matrix()
    {}

    // construct from components
    Matrix(const Vector& x, const Vector& y, const Vector& z)
    :
      x_(x), y_(y), z_(z)
    {}

    // destructor
    ~Matrix()
    {}

    // acess functions
    const Vector& x() const { return x_; }
    const Vector& y() const { return y_; }
    const Vector& z() const { return z_; }

    // determinant
    double det() const
    {
        double a = y_.y()*z_.z() - y_.z()*z_.y();
        double b = -(y_.x()*z_.z() - y_.z()*z_.x());
	double c = y_.x()*z_.y() - y_.y()*z_.x();
        return x_.x()*a + x_.y()*b + x_.z()*c;;
    }

    // transpose
    Matrix transpose() const
    {
        Vector a(x_.x(), y_.x(), z_.x());
	Vector b(x_.y(), y_.y(), z_.y());
	Vector c(x_.z(), y_.z(), z_.z());
	return Matrix(a, b, c);
    }

    // inverse
    Matrix inverse() const
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
    const Vector operator&(const Vector& b)
    {
        double x = x_ & b;
	double y = y_ & b;
	double z = z_ & b;
	return Vector(x, y, z);
    }

    const Matrix operator*(const Matrix& m) const
    {
        Matrix bT = m.transpose();
	Vector a(x_ & bT.x(), x_ & bT.y(), x_ & bT.z());
	Vector b(y_ & bT.x(), y_ & bT.y(), y_ & bT.z());
	Vector c(z_ & bT.x(), z_ & bT.y(), z_ & bT.z());
	return Matrix(a,b,c);
    }

    friend ostream& operator<<(ostream& os, const Matrix& v)
    {
        os << v.x() << endl << v.y() << endl << v.z();
        return os;
    }

};


int main(int argc, char* argv[]) 
{

    if (argc < 4)
    {
	cout << "wrong number of arguments!" << endl;
	cout << argv[0] << " <input filename> <nIntervals phi> <nIntervals> <output filename>" << endl;
	return 0;
    }

    // read the command line arguments
    char* inputFileName(argv[1]);
    int nPhi = atoi(argv[2]);
    int nZ = atoi(argv[3]);
    char* outputFileName(argv[4]);
    
    vector<Vector> x1(0), x2(0), x3(0), u(0), xc(0), n(0);
    vector<double> area(0);

    // try to open the input file
    ifstream iFile(inputFileName);
    if (!iFile.is_open())
    {
	cout << "Could not open file " << inputFileName << endl;
	return 0;
    }
    else
    {
	cout << "Reading " << inputFileName << endl;
    }
    
    // read the input file
    int nLines = 0;
    while (iFile.good())
    {
	double tmp;
	Vector a(iFile);
	
	if (iFile.good())
	{
	    nLines++;
	    Vector b(iFile);
	    Vector c(iFile);
	    Vector ut(iFile);
	    Vector xct = (a+b+c)/3.0;
	    Vector normal = (b-a)^(c-a);
	    double dotsqrt = sqrt(normal & normal);
	    
	    normal = normal / (dotsqrt + SMALL);

	    x1.push_back(a);
	    x2.push_back(b);
	    x3.push_back(c);
	    u.push_back(ut);
	    xc.push_back(xct);
	    n.push_back(normal);
	    area.push_back(0.5*dotsqrt);
	}
    }

    cout << "nLines = " << nLines << endl;

    Vector center(0,0,0);
    double totArea = 0;

    // calculate center point
    for(int i=0; i<nLines; i++)
    {
	center = center + (area[i]*xc[i]);
	totArea += area[i];
    }
    center = center / totArea;
    cout << "center = " << center << endl;

    // correct the element normals to point outwards
    for(int i=0; i<nLines; i++)
    {
        double dot = (center - xc[i]) & n[i];
	if (dot > 0)
	{
	    n[i] = -1.0*n[i];
	}
    }
    /*
    Matrix mm(n[1], n[10], n[100]);
    Matrix mInv = mm.inverse();
    Matrix mT = mm.transpose();
    Matrix id = mm * mInv;
    */

    Vector xTop(0, 0, 0);
    int num = 0;
    int skip = 60;
    for (int i=0; i<nLines-1; i += skip)
    {
        if (area[i] > SMALL)
	{
	  //cout << "i = " << i << endl;
	    for (int j=i+1; j<nLines; j += skip)
	    {
	        if (area[j] > SMALL)
		{
		    for (int k=0; k<nLines; k += skip)
		    {
		        if ( (area[k] > SMALL) && ( k != i) && (k != j) )
			{
			    Matrix m(n[i], n[j], n[k]);
			    if (m.det() > 1.0e-4)
			    {
			        Matrix mInv = m.inverse();
				Vector rhs(xc[i] & n[i], xc[j] & n[j], xc[k] & n[k]);
				Vector xt = mInv&rhs;
				Vector diff = xTop - xt;
				xTop = num*xTop + xt;
				num++;
				xTop = xTop/num;
			    }
			}
		    }
		}
	    }
	}
    }

    cout << "xTop = " << xTop << endl;
    cout << "nCombinations = " << num << endl;
    return 0;
}
