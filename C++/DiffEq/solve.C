/*
  This program is intended for two things
  Solve the differential equation using 4-5 th order Runge-Kutta method of:
  1. The spring-mass-damper equation with a random force
  2. Solve the n-body problem.
*/

#include <iostream>
#include <vector>
#include "Vector.H"

typedef double scalar;
typedef int label;

class nBody
{
private:
  label nBodies_;
  vector<Vector> positions_;
  vector<Vector> velocities_;
  vector<scalar> mass_;

public:
  nBody(const label nBodies)
  :
    nBodies_(nBodies),
    positions_(nBodies),
    velocities_(nBodies),
    mass_(nBodies)
  {}

  ~nBody()
  {}

  const vector<Vector>& positions() const
  {
    return positions_;
  }

  vector<Vector>& positions()
  {
    return positions_;
  }

};

class spring
{
};

using namespace std;

int main()
{
  scalar a = 1.0;
  label nBodies = 3;
  nBody system(nBodies);

  cout << "hello world, a = " << a << endl;
  
  label size = system.positions().size();

  cout << "size = " << size << endl;

  return 0;
}
