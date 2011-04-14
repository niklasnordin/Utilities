vector v1 = p1().U();
vector v2 = p2().U();

vector vRel = v1 - v2;
scalar magVRel = mag(vRel);

scalar sumD = p1().d() + p2().d();
scalar pc = spray_.p()[p1().cell()];

spray::iterator pMin = p1;
spray::iterator pMax = p2;

scalar dMin = pMin().d();
scalar dMax = pMax().d();

if (dMin > dMax) {
    dMin = pMax().d();
    dMax = pMin().d();
    pMin = p2;
    pMax = p1;
}

scalar rhoMax = spray_.fuels().rho(pc, pMax().T(), pMax().X());
scalar rhoMin = spray_.fuels().rho(pc, pMin().T(), pMin().X());
scalar mMax = pMax().m();
scalar mMin = pMin().m();
scalar mTot = mMax + mMin;

scalar nMax = pMax().N(rhoMax);
scalar nMin = pMin().N(rhoMin);

scalar mdMin = mMin/nMin;

scalar nu0 = 0.25*mathematicalConstant::pi*sumD*sumD*magVRel*dt/vols_[cell1];
scalar nu = nMin*nu0;
scalar collProb = exp(-nu);
scalar xx = rndGen_.scalar01();

// collision occur
if (( xx > collProb) && (mMin > VSMALL) && (mMax > VSMALL)) {

    scalar gamma = dMax/max(dMin, 1.0e-12);
    scalar f = gamma*gamma*gamma + 2.7*gamma - 2.4*gamma*gamma;

    vector momMax = mMax*pMax().U();
    vector momMin = mMin*pMin().U();

//------------------------------------AL_____101011--------------------------------------------//
    // use mass-averaged temperature to calculate We number
//    scalar averageTemp = (pMax().T()*mMax + pMin().T()*mMin)/mTot;
    scalar averageTemp = (pMax().T()*dMax*dMax*dMax + pMin().T()*dMin*dMin*dMin)/(dMax*dMax*dMax + dMin*dMin*dMin);
    Info << "average temperature = " << averageTemp << endl;
    // and mass averaged mole fractions ...
//    scalarField Xav((pMax().m()*pMax().X()+pMin().m()*pMin().X())/(pMax().m() + pMin().m()));
//  simple my improvement:
//  scalarField Xav((pMax().X()*dMax*dMax*dMax + pMin().X()*dMin*dMin*dMin)/(dMax*dMax*dMax + dMin*dMin*dMin));
//  more my accurate improvement: 
    scalarField Ymin = spray_.fuels().Y(pMin().X());		//if this line is activated, the same line should be commented in the following
    scalarField Ymax = spray_.fuels().Y(pMax().X());        	//if this line is activated, the same line should be commented in the following
    scalarField Yav = (mMax*Ymax/nMax + mMin*Ymin/nMin)/(mMax/nMax + mMin/nMin);
    scalarField Xav = Yav;
    scalar Wav = 0.0;
    forAll(Yav, i)
        {
            Wav += Yav[i]/spray_.fuels().properties()[i].W();
        }
    forAll(Yav, i)
        {
            Xav[i] = Yav[i]/(spray_.fuels().properties()[i].W()*Wav);
        }
    Info << "Xav = " << Xav << endl;
//---------------------------------------------------------------------------------------------//

    scalar sigma = spray_.fuels().sigma(pc, averageTemp, Xav);
    sigma = max(1.0e-6, sigma);
    Info << "sigma = " << sigma << endl;
//------------------------------------AL_____101011--------------------------------------------//
//    scalar rho = spray_.fuels().rho(pc, averageTemp, Xav);
    scalar rho = (rhoMax*dMax*dMax*dMax + rhoMin*dMin*dMin*dMin)/(dMax*dMax*dMax + dMin*dMin*dMin);
    Info << "rho = " << rho << endl;
//---------------------------------------------------------------------------------------------//

    scalar WeColl = max(1.0e-12, 0.5*rho*magVRel*magVRel*dMin/sigma);
    Info << "WeColl = " << WeColl << endl; 
    scalar coalesceProb = min(1.0, 2.4*f/WeColl);
    scalar prob = rndGen_.scalar01();

    // Coalescence
    if ( prob < coalesceProb && coalescence_) {

        // How 'many' of the droplets coalesce
        // NN. This is the kiva way ... which actually works best

        scalar zz = collProb;
        scalar vnu = nu*collProb;
        label n=2;

        // xx > collProb=zz
        while ((zz < xx) && (n<1000)) {
            zz += vnu;
            vnu *= nu/n;
	    n++;
        }
        Info<< "vnu = " << vnu << ", n = " << n << endl;
        scalar nProb = n - 1;
        // All droplets coalesce
        if (nProb*nMax > nMin) {
            nProb = nMin/nMax;
        }

        // Conservation of mass, momentum and energy
        pMin().m() -= nProb*nMax*mdMin;

        scalar newMinMass = pMin().m();
        scalar newMaxMass = mMax + (mMin - newMinMass);
        pMax().m() = newMaxMass;

//------------------------------------AL_____101011--------------------------------------------//
//        pMax().T() = (averageTemp*mTot - newMinMass*pMin().T())/newMaxMass;
//        rhoMax = spray_.fuels().rho(pc, pMax().T(), pMax().X());
//        scalar d3 = pow(dMax, 3) + nProb*pow(dMin,3);
//        pMax().d() = cbrt(d3);
//        pMax().U() = (momMax + (1.0-newMinMass/mMin)*momMin)/newMaxMass;
	scalar factorCollAv = 1./(dMax*dMax*dMax + nProb*dMin*dMin*dMin);
	pMax().T() = (pMax().T()*dMax*dMax*dMax + pMin().T()*dMin*dMin*dMin*nProb)*factorCollAv;
	pMax().U() = (pMax().U()*dMax*dMax*dMax + pMin().U()*dMin*dMin*dMin*nProb)*factorCollAv; 
        rhoMax = spray_.fuels().rho(pc, pMax().T(), pMax().X());
        scalar d3 = pow(dMax, 3) + nProb*pow(dMin,3);
        pMax().d() = cbrt(d3);
        Info << "pMax().T() = " << pMax().T() << ", pMax().U()" << pMax().U() << ", rhoMax = " << rhoMax << ", pMax().d() = " << pMax().d() << endl;
//---------------------------------------END---------------------------------------------------//        
        // update the liquid molar fractions
//      scalarField Ymin = spray_.fuels().Y(pMin().X());
//      scalarField Ymax = spray_.fuels().Y(pMax().X());
// calculate the mass fraction Ynew in the big droplet, which is the new droplet.
        scalarField Ynew = mMax*Ymax + (mMin - newMinMass)*Ymin;
        Info << "(mMax + mMin - newMinMass)= " << mMax + mMin - newMinMass << endl;
//------------------------------------AL_____101011--------------------------------------------//
	Ynew = Ynew/(mMax + mMin - newMinMass + SMALL);
        Info << "Ynew = " << Ynew << endl;
//---------------------------------------END---------------------------------------------------//
        scalar Wlinv = 0.0;
        forAll(Ynew, i)
        {
            Wlinv += Ynew[i]/spray_.fuels().properties()[i].W();
        Info << "Wlinv = " << Wlinv << endl;
        }
        forAll(Ynew, i)
        {
            pMax().X()[i] = Ynew[i]/(spray_.fuels().properties()[i].W()*Wlinv);
        }
        Info << "pMax().X() = " << pMax().X() << endl;
    }
    // Grazing collision (no coalescence)
    else {

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
        scalar rho2 = spray_.fuels().rho(pc, p2().T(), p2().X());
        scalar m1 = p1().m();
        scalar m2 = p2().m();
        scalar n1 = p1().N(rho1);
        scalar n2 = p2().N(rho2);

        // gf -> 1 => v1p -> p1().U() ...
        // gf -> 0 => v1p -> momentum/(m1+m2)
//------------------------------------AL_____101011--------------------------------------------//
//      vector mr = m1*v1 + m2*v2;
//      vector v1p = (mr + m2*gf*vRel)/(m1+m2);
//      vector v2p = (mr - m1*gf*vRel)/(m1+m2);
	scalar vol1 = m1/(rho1*n1);
	scalar vol2 = m2/(rho2*n2);
	vector v1p = (v1*vol1 + (v2+gf*vRel)*vol2)/(vol1 + vol2);
        vector v2p = ((v1-gf*vRel)*vol1 + v2*vol2)/(vol1 + vol2);
//---------------------------------------------------------------------------------------------//

        if (n1 < n2) {
            p1().U() = v1p;
            p2().U() = (n1*v2p + (n2-n1)*v2)/n2;
        }
        else {
            p1().U() = (n2*v1p + (n1-n2)*v1)/n1;
            p2().U() = v2p;
        }

    } // if - coalescence or not

} // if - collision

