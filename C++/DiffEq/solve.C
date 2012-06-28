/*
  This program is intended for two things
  Solve the differential equation using 4-5 th order Runge-Kutta method of:
  1. The spring-mass-damper equation with a random force
  2. Solve the n-body problem.
*/

#include <iostream>
#include <vector>
#include <math.h>
#include "Vector.H"

typedef double scalar;
typedef int label;

#define SMALL 1.0e-15

// gravitational constant
//#define G 6.672e-11
#define G 6.672

class nBody
{
private:
  label nBodies_;
  vector<Vector> positions_;
  vector<Vector> velocities_;
  vector<scalar> mass_;

public:

  // constructor
  nBody(const label nBodies)
  :
    nBodies_(nBodies),
    positions_(nBodies),
    velocities_(nBodies),
    mass_(nBodies)
  {}

  // destructor
  ~nBody()
  {}

  // access functions
  const vector<Vector>& positions() const
  {
    return positions_;
  }

  vector<Vector>& positions()
  {
    return positions_;
  }

  const vector<Vector>& velocities() const
  {
    return velocities_;
  }

  vector<Vector>& velocities()
  {
    return velocities_;
  }

  const vector<scalar>& mass() const
  {
    return mass_;
  }
  
  vector<scalar>& mass()
  {
    return mass_;
  }

  vector<Vector> forces
  (
    const vector<Vector>& positions
  )
  {
    Vector zero = Vector(0,0,0);

    vector<Vector> forces(positions.size());
    for(label i=0; i<positions.size(); i++)
    {
	forces[i] = zero;
    }

    for(label i=0; i<positions.size()-1; i++)
    {
	for(label j=i+1; j<positions.size(); j++)
	{
	    Vector rad = positions[i] - positions[j];
	    
	    scalar rad2 = rad & rad;
	    scalar mag = sqrt(rad2);
	    Vector n = rad/(mag + SMALL);

	    Vector fij = G * mass_[i] * mass_[j] * n / ( rad2 + SMALL );
	    forces[i] = forces[i] + fij;
	    forces[j] = forces[j] - fij;
	}
    } 

    return forces;
  }


  void evolve(const scalar dt)
  {
    vector<Vector> f = forces(positions());

    // update velocities
    // m a = F
    for(label i=0; i<nBodies_; i++)
    {
	velocities_[i] = velocities_[i] + dt*f[i]/mass_[i];
    }

    // update positions
    for(label i=0; i<nBodies_; i++)
    {
	positions_[i] = positions_[i] + dt*velocities_[i];
    }
        
  }
  
  void initMasses
  (
    const double minMass,
    const double maxMass
  )
  {
    for(label i=0; i<mass_.size(); i++)
      {
	scalar mass = (maxMass-minMass)*((scalar)rand()/(scalar)RAND_MAX) + minMass;
	mass_[i] = mass;
      }
  }
  
  void randomizeVectors
  (
    vector<Vector>& v, 
    const double& mag
  )
  {
    for(label i=0; i<v.size(); i++)
      {
	scalar x = ((scalar)rand()/(scalar)RAND_MAX) - 0.5;
	scalar y = ((scalar)rand()/(scalar)RAND_MAX) - 0.5;
	scalar z = ((scalar)rand()/(scalar)RAND_MAX) - 0.5;
	z = 0.0;
	Vector p(x, y, z);
	scalar pMag = p.mag() + SMALL;
	v[i] = p*mag/pMag;
      }
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
  scalar minMass = 1.0;
  scalar maxMass = 1000.0;
  scalar maxRadius = 10000.0;
  scalar maxVelocity = 10.0;
  scalar endTime = 100.0;
  scalar deltaT = 1.0e-1;
  nBody system(nBodies);
  system.initMasses(minMass, maxMass);
  system.randomizeVectors(system.positions(), maxRadius);
  system.randomizeVectors(system.velocities(), maxVelocity);

  scalar time = 0.0;
  while (time <= endTime)
    {
      time += deltaT;
      system.evolve(deltaT);
      cout << "t = " << time << ": "
	   << system.positions()[0] << " : "
	   << system.velocities()[0]
	   << endl;
    }
  return 0;
}
