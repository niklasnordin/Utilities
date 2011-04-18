/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "TrajectoryCollision.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class CloudType>
Foam::TrajectoryCollision<CloudType>::TrajectoryCollision
(
    const dictionary& dict,
    CloudType& owner
)
:
    CollisionModel<CloudType>(owner),
    coeffsDict_(dict.subDict(typeName + "Coeffs")),
    cSpace_(readScalar(coeffsDict_.lookup("cSpace"))),
    cTime_(readScalar(coeffsDict_.lookup("cTime"))),
    coalescence_(coeffsDict_.lookup("coalescence"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template <class CloudType>
Foam::TrajectoryCollision<CloudType>::~TrajectoryCollision()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::TrajectoryCollision<CloudType>::active() const
{
    return true;
}

template<class CloudType>
bool Foam::TrajectoryCollision<CloudType>::update
(
    const scalar& dt,
    Random& rndGen,
    vector& pos1,
    scalar& m1,
    scalar& d1,
    scalar& N1,
    vector& U1,
    scalar& rho1,
    scalar& T1,
    scalarField& Y1,
    const scalar& sigma1,
    const label celli,
    const scalar voli,
    vector& pos2,
    scalar& m2,
    scalar& d2,
    scalar& N2,
    vector& U2,
    scalar& rho2,
    scalar& T2,
    scalarField& Y2,
    const scalar& sigma2,
    const label cellj,
    const scalar volj
) const
{
    bool coalescense = false;

    vector vRel = U1 - U2;
    //    scalar magVRel = mag(vRel);
    
    vector p = pos2 - pos1;
    scalar dist = mag(p);
    
    scalar vAlign = vRel & (p/(dist+SMALL));

    if (vAlign > 0)
    {
        scalar sumD = d1 + d2;
        
        if (vAlign*dt > dist - 0.5*sumD)
        {
            scalar v1Mag = mag(U1);
            scalar v2Mag = mag(U2);
            vector nv1 = U1/v1Mag;
            vector nv2 = U2/v2Mag;
            
            scalar v1v2 = nv1 & nv2;
            scalar v1p = nv1 & p;
            scalar v2p = nv2 & p;
            
            scalar det = 1.0 - v1v2*v1v2;
            
            scalar alpha = 1.0e+20;
            scalar beta = 1.0e+20;
            
            if (mag(det) > 1.0e-4)
            {
                beta = -(v2p - v1v2*v1p)/det;
                alpha = v1p + v1v2*beta;
            }

            alpha /= v1Mag*dt;
            beta /= v2Mag*dt;
            
            // is collision possible within this timestep
            if ((alpha>0) && (alpha<1.0) && (beta>0) && (beta<1.0)) 
            {
                vector p1c = pos1 + alpha*U1*dt;
                vector p2c = pos2 + beta*U2*dt;
                
                scalar closestDist = mag(p1c-p2c);
                
                scalar collProb = 
                    pow(0.5*sumD/max(0.5*sumD, closestDist), cSpace_)
                  * exp(-cTime_*mag(alpha-beta));

                scalar xx = rndGen.scalar01();

                // collision occur
                if ((xx < collProb) && (m1 > VSMALL) && (m2 > VSMALL))
                {
                    if (d1 > d2)
                    {
                        coalescense = collideSorted
                        (
                            dt,
                            rndGen,
                            pos1,
                            m1,
                            d1,
                            N1,
                            U1,
                            rho1,
                            T1,
                            Y1,
                            sigma1,
                            celli,
                            voli,
                            pos2,
                            m2,
                            d2,
                            N2,
                            U2,
                            rho2,
                            T2,
                            Y2,
                            sigma2,
                            cellj,
                            volj
                        );
                    }
                    else
                    {
                        coalescense = collideSorted
                        (
                            dt,
                            rndGen,
                            pos2,
                            m2,
                            d2,
                            N2,
                            U2,
                            rho2,
                            T2,
                            Y2,
                            sigma2,
                            cellj,
                            volj,
                            pos1,
                            m1,
                            d1,
                            N1,
                            U1,
                            rho1,
                            T1,
                            Y1,
                            sigma1,
                            celli,
                            voli
                        );
                    }
                } // if ( d1 > d2 )
                
            } // if - possible collision (alpha, beta) in timeinterval
            
        } // if - travelled distance is larger distance between parcels

    }

    return coalescense;
}


template<class CloudType>
bool Foam::TrajectoryCollision<CloudType>::collideSorted
(
    const scalar& dt,
    Random& rndGen,
    vector& pos1,
    scalar& m1,
    scalar& d1,
    scalar& N1,
    vector& U1,
    scalar& rho1,
    scalar& T1,
    scalarField& Y1,
    const scalar& sigma1,
    const label celli,
    const scalar voli,
    vector& pos2,
    scalar& m2,
    scalar& d2,
    scalar& N2,
    vector& U2,
    scalar& rho2,
    scalar& T2,
    scalarField& Y2,
    const scalar& sigma2,
    const label cellj,
    const scalar volj
) const
{
    bool coalescence = false;

    vector vRel = U1 - U2;
    scalar magVRel = mag(vRel);

    scalar mdMin = m2/N2;

    scalar mTot = m1 + m2;

    scalar gamma = d1/max(d2, 1.0e-12);
    scalar f = gamma*gamma*gamma + 2.7*gamma - 2.4*gamma*gamma;

    vector momMax = m1*U1;
    vector momMin = m2*U2;

    // use mass-averaged temperature to calculate We number
    scalar Tm = (T1*m1 + T2*m2)/mTot;

    // and mass averaged fractions ...
    //scalarField Yav((m1*Y1 + m2*Y2)/mTot;

    // interpolate the averaged surface tension
    scalar sigma = sigma1 + (sigma2 - sigma1)*(Tm - T1)/(T2 - T1);

    sigma = max(1.0e-6, sigma);
    scalar Vtot = m1/rho1 + m2/rho2;
    scalar rho = mTot/Vtot;

    scalar dMean = sqrt(d1*d2);
    scalar WeColl = max(1.0e-12, 0.5*rho*magVRel*magVRel*dMean/sigma);

    scalar coalesceProb = min(1.0, 2.4*f/WeColl);
    
    scalar prob = rndGen.scalar01();

    // Coalescence
    if ( prob < coalesceProb && coalescence_) 
    {
        coalescense = true;
        // How 'many' of the droplets coalesce
        scalar nProb = prob*N2/N1;

        // Conservation of mass, momentum and energy
        scalar m2Org = m2;
        m2 -= N1*nProb*mdMin;

        scalar newMaxMass = m1 + (m2Org - m2);
        m1 = newMaxMass;
        
        T1 = (Tm*mTot - m2*T2)/m1;

                    pMax().d() = 
                        pow
                        (
                            6.0*newMaxMass/(rhoMax*mathematicalConstant::pi*nMax),
                            1.0/3.0
                        );

                    pMax().U() = 
                        (momMax + (1.0-newMinMass/mMin)*momMin)/newMaxMass;

                    // update the liquid molar fractions
                    scalarField Ymin = spray_.fuels().Y(pMin().X());
                    scalarField Ymax = spray_.fuels().Y(pMax().X());
                    scalarField Ynew = mMax*Ymax + (mMin - newMinMass)*Ymin;
                    scalar Wlinv = 0.0;
                    forAll(Ynew, i)
                    {
                        Wlinv += Ynew[i]/spray_.fuels().properties()[i].W();
                    }
                    forAll(Ynew, i)
                    {
                        pMax().X()[i] = 
                            Ynew[i]/(spray_.fuels().properties()[i].W()*Wlinv);
                    }


                }
                // Grazing collision (no coalescence)
                else
                {
                    scalar gf = sqrt(prob) - sqrt(coalesceProb);
                    scalar denom = 1.0 - sqrt(coalesceProb);
                    if (denom < 1.0e-5) {
                        denom = 1.0;
                    }
                    gf /= denom;

                    // if gf negative, this means that coalescence is turned off
                    // and these parcels should have coalesced
                    gf = max(0.0, gf);

                    scalar rho1 = spray_.fuels().rho(pc, p1().T(), p1().X());
                    scalar rho2 = spray_.fuels().rho(0.0, p2().T(), p2().X());
                    scalar m1 = p1().m();
                    scalar m2 = p2().m();
                    scalar n1 = p1().N(rho1);
                    scalar n2 = p2().N(rho2);

                    // gf -> 1 => v1p -> p1().U() ...
                    // gf -> 0 => v1p -> momentum/(m1+m2)

                    vector mr = m1*v1 + m2*v2;
                    vector v1p = (mr + m2*gf*vRel)/(m1+m2);
                    vector v2p = (mr - m1*gf*vRel)/(m1+m2);

                    if (n1 < n2)
                    {
                        p1().U() = v1p;
                        p2().U() = (n1*v2p + (n2-n1)*v2)/n2;
                    }
                    else
                    {
                        p1().U() = (n2*v1p + (n1-n2)*v1)/n1;
                        p2().U() = v2p;
                    }

                } // if - coalescence or not

            } // if - collision
        */

    return coalescence;
}

// ************************************************************************* //
