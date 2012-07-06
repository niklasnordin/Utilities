/*
  This program is intended for two things
  Solve the differential equation using 4-5 th order Runge-Kutta method of:
  1. The spring-mass-damper equation with a random force
  2. Solve the n-body problem.
*/

#include <iostream>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include "Vector.H"

typedef double scalar;
typedef int label;

#define SMALL 1.0e-15

// gravitational constant
//#define G 6.672e-11
#define G 6.672e-11

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
        const vector<Vector>& pos
    )
    {
        Vector zero = Vector(0,0,0);

        vector<Vector> forces(pos.size());
        for(label i=0; i<pos.size(); i++)
        {
	    forces[i] = zero;
        }

	scalar rad2Min = 1.0e-1;

        for(label i=0; i<pos.size()-1; i++)
        {
	    for(label j=i+1; j<pos.size(); j++)
	    {
	        Vector rad = pos[i] - pos[j];
	    
                scalar rad2 = rad & rad;
                scalar mag = sqrt(rad2);
                Vector n = rad/(mag + SMALL);
		if (rad2 < rad2Min) rad2 = rad2Min;
                Vector fij = (G * mass_[i] * mass_[j] / rad2 ) * n;
                forces[i] = forces[i] - fij;
                forces[j] = forces[j] + fij;
            }
        } 

        return forces;
    }


    void evolve(const scalar dt)
    {

	Vector center = Vector(0,0,0);
	scalar totMass = 0.0;
	for(label i=0; i<nBodies_; i++)
	{
	    totMass += mass_[i];
	    center = center + mass_[i]*positions_[i];
	}
	center = center/totMass;

        vector<Vector> k1 = forces(positions_);

	vector<Vector> v1 = velocities_;
	vector<Vector> v2(k1.size());
	vector<Vector> v3(k1.size());
	vector<Vector> v4(k1.size());

	vector<Vector> p1(k1.size());
	vector<Vector> p2(k1.size());
	vector<Vector> p3(k1.size());
	vector<Vector> p4(k1.size());

        // update velocities and positions
        // m a = F
        for(label i=0; i<nBodies_; i++)
        {
            v2[i] = velocities_[i] + 0.5*dt*k1[i]/mass_[i];
            p2[i] = positions_[i] + 0.5*dt*v1[i];
        }

        vector<Vector> k2 = forces(p2);

        for(label i=0; i<nBodies_; i++)
        {
            v3[i] = velocities_[i] + 0.5*dt*k2[i]/mass_[i];
            p3[i] = positions_[i] + 0.5*dt*v2[i];
        }

        vector<Vector> k3 = forces(p3);

        for(label i=0; i<nBodies_; i++)
        {
            v4[i] = velocities_[i] + dt*k3[i]/mass_[i];
            p4[i] = positions_[i] + dt*v3[i];
        }

        vector<Vector> k4 = forces(p4);
	/*
        for(label i=0; i<nBodies_; i++)
        {
	    velocities_[i] = velocities_[i] + dt*k1[i]/mass_[i];
	    positions_[i] = positions_[i] + dt*velocities_[i];
        }
	*/
        for(label i=0; i<nBodies_; i++)
        {
	    velocities_[i] = velocities_[i] + dt*(k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i])/(6.0*mass_[i]);
            positions_[i] = positions_[i] + dt*(v1[i] + 2.0*v2[i] + 2.0*v3[i] + v4[i])/6.0;
        }


	// transpose position back to mass center
	Vector newCenter = Vector(0,0,0);
        for(label i=0; i<nBodies_; i++)
        {
	    newCenter = newCenter + mass_[i]*positions_[i];
        }
	newCenter = newCenter/totMass;


	Vector transpose = center - newCenter;
	for(label i=0; i<nBodies_; i++)
	{
	  positions_[i] = positions_[i] + transpose;
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
    label nBodies = 5;
    scalar minMass = 1.0e+10;
    scalar maxMass = 2.0e+10;
    scalar maxRadius = 100.0;
    scalar maxVelocity = 0.01;
    label nSteps = 100000000;
    scalar deltaT = 1.0e-4;
    label nOut = nSteps / 100000;

    nBody system(nBodies);
    system.initMasses(minMass, maxMass);
    system.randomizeVectors(system.positions(), maxRadius);
    system.randomizeVectors(system.velocities(), maxVelocity);

    for(label i=0; i<nBodies; i++)
    {
        cout << "m[" << i << "] = " << system.mass()[i] 
             << ", X = " << system.positions()[i]
             << ", U = " << system.velocities()[i] 
             << endl;
    }

    ofstream file("pos.dat");
    file.precision(10);
    label nCount = 0;
    for(label ii=0; ii<nSteps; ii++)
    {
        nCount++;
        system.evolve(deltaT);

	if (nCount >= nOut)
	{
	    nCount = 0;

            cout << "t = " << time << ": "
                 << system.positions()[0] << " : "
                 << system.velocities()[0].mag()
                 << endl;

	    for(label i=0; i<nBodies; i++)
	    {
	      file << system.positions()[i] << " ";
	    }
	    file << endl;
	}
    }
    file.close();

    cout << "==================" << endl;
    for(label i=0; i<nBodies; i++)
    {
        cout << "m[" << i << "] = " << system.mass()[i] 
             << ", X = " << system.positions()[i] 
             << ", U = " << system.velocities()[i] 
             << endl;
    }
  
    return 0;
}
