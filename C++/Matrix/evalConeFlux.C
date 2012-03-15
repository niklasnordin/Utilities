#include <iostream>
#include <string.h>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <math.h>

using namespace std;

#define SMALL 1.0e-15
#define pi 3.14159265358979323

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

    if (argc < 6)
    {
	cout << "wrong number of arguments!" << endl;
	cout << argv[0] << " <input filename> <delta phi max> <nIntervals> <min dist> <max dist> <output filename>" << endl;
	return 0;
    }

    // read the command line arguments
    char* inputFileName(argv[1]);
    // maximum angle between different bins
    int dPhi = atof(argv[2]);
    int nZ = atoi(argv[3]);
    double dMin = atof(argv[4]); // 0.02
    double dMax = atof(argv[5]); // 0.08
    char* outputFileName(argv[6]);
    
    vector<Vector> x1(0), x2(0), x3(0), u(0), xc(0), n(0);
    vector<double> area(0), phi(0), dist(0), angleBin(0);

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

    // find the top point of the cone
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

    // find the axis of rotation
    Vector z(0, 0, 0);
    for (int i=0; i<nLines; i++)
    {
        Vector di = xc[i] - xTop;
	
        if ( (area[i] > SMALL) && (di.mag() > dMin) && (di.mag() < dMax) )
	{
	    z = z + n[i]*area[i];
	}
    }
    z = z/-z.mag();

    cout << "xTop = " << xTop << endl;
    cout << "z = " << z << endl;
    cout << "nCombinations = " << num << endl;

    // choose two vectors perpendicular to z
    Vector x(0, 1, 0);
    Vector y;
    double xz = x & z;
    if (xz > 1.0e-4)
    {
        // this vector is good enough, remove the z component
        x = x - (xz*z);
	x = x/x.mag();
	y = z ^ x;
    }
    else
    {
        // if (0,1,0) is not good enough, this will be
        x = Vector(0, 0, 1);
        x = x - (xz*z);
	x = x/x.mag();
	y = z ^ x;
    }
    cout << "0 deg is at direction x = ";
    cout << x << endl;
    cout << "90 deg is at direction y = ";
    cout << y << endl;

    // transform the points to the local coordinate system and calculate the 'angle'
    for (int i=0; i<nLines; i++)
    {
        Vector lox = xc[i] - xTop;
	double xLoc = lox & x;
	double yLoc = lox & y;
	dist.push_back(lox & z);
	double alpha = atan(yLoc/xLoc)*180.0/pi;
	if (xLoc < 0) alpha += 180.0;
	if (alpha < 0) alpha += 360.0;
	phi.push_back(alpha);

    }
    
    ofstream stlFile("geometry.stl");
    stlFile << "solid geometry" << endl;
    stlFile << "  color 0.5 0 0.5" << endl;
    for (int i=0; i<nLines; i++)
    {
        if (dist[i] > dMin && dist[i] < dMax)
	{	  
            stlFile << "  facet normal " << n[i] << endl;
	    stlFile << "    outer loop" << endl;
	    stlFile << "      vertex " << x1[i] << endl;
	    stlFile << "      vertex " << x2[i] << endl;
	    stlFile << "      vertex " << x3[i] << endl;
	    stlFile << "    endloop" << endl;
	    stlFile << "  endfacet" << endl;

	    bool needNewBin = true;
	    int n = angleBin.size();
	    for (int j=0; j<n; j++)
	    {
	        double d1 = fabs(phi[i] - angleBin[j]);
	        double d2 = fabs(360.0+phi[i] - angleBin[j]);
		double diff = min(d1, d2);
		if (diff < dPhi)
		{
		    needNewBin = false;
		}
	    }

	    if (needNewBin)
	    {
	        angleBin.push_back(phi[i]);
	    }
	}
    }
    stlFile << "endsolid geometry" << endl;

    cout << "number of slits = " << angleBin.size() << endl;

    double alphaMin = 1.0e+10;
    for (int i=0; i<angleBin.size(); i++)
    {
        alphaMin = min(angleBin[i], alphaMin);
    }

    // find the offset angle
    double aOffset = 0.0;
    int counter = 0;
    for (int i=0; i<nLines; i++)
    {
        if (dist[i] > dMin && dist[i] < dMax)
	{	  
	    double d1 = fabs(phi[i] - alphaMin);
	    double d2 = fabs(360.0+phi[i] - alphaMin);
	    double diff = min(d1, d2);
	    if (diff < dPhi)
	    {
		counter++;
		aOffset += phi[i];
	    }
	}
    }

    // recalculate the angles to have them in ordered fashion
    aOffset /= counter;
    //cout << "alpha offset = " << aOffset << endl;
    for (int i=0; i<angleBin.size(); i++)
    {
        angleBin[i] = i*360.0/angleBin.size() + aOffset;
    }

    double flux[angleBin.size()], slitArea[angleBin.size()];

    // calculate the volume flux thru each slit
    for (int i=0; i<nLines; i++)
    {
        if (dist[i] > dMin && dist[i] < dMax)
	{
	    double fluxi = area[i]*(u[i] & n[i]);
	    
	    for(int j=0; j<angleBin.size(); j++)
	    {
		double d1 = fabs(phi[i] - angleBin[j]);
		double d2 = fabs(360.0+phi[i] - angleBin[j]);
		double diff = min(d1, d2);
		if (diff < dPhi)
		{
		    flux[j] += fluxi;
		    slitArea[j] += fabs(area[i]);
		}
	    }
	}
    }

    ofstream oFile(outputFileName);

    for(int i=0; i<angleBin.size(); i++)
    {
        cout << "a = " << angleBin[i] << ", flux = " << flux[i] 
	     << ", area = " << slitArea[i] 
	     << ", averageFlux = " << flux[i]/slitArea[i]
	     << endl;
	oFile << angleBin[i] << " " << -flux[i] << endl;
    }

    return 0;
}
