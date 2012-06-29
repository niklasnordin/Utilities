#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <cstdlib>
#include <list>
#include <vector>
#include <math.h>

using namespace std;

#include "function.H"

int main(int argc, char *argv[])
{
  //int nFunction = 2;

    const double dx = 1.0e-5;
    //    const double urlx = 1.0e-3;
    //const int maxIter = 100000;
    const int output = 0;
    const double errMax = 1.0e-10;

    int nEntries;
    vector<double> func;

    cout.precision(5);

    if (argc < 8)
    {
        cout << "missing input of parameter and result file. Exiting" << endl;
	cout << "Usage:\n\t " << argv[0] << " <data file> <nParameters> <nFunctions> <targetFunction> <urlx> <itMax> <errDiffMax> <coeffs>" << endl;
	exit(1);
    }

    string fileName(argv[1]);  
    //string resFileName(argv[2]);

    int nParameters = atoi(argv[2]);
    int nResponses = atoi(argv[3]);
    int nFunction = atoi(argv[4]) - 1;
    double urlx = atof(argv[5]);
    int maxIter = atoi(argv[6]);
    double errDiff = atof(argv[7]);

    // use this if you are including quadratic terms
    int nCoeffs = 1 + 2*nParameters + (nParameters-1)*nParameters/2;
    //int nCoeffs = 1 + nParameters + (nParameters-1)*nParameters/2;
    double A[nCoeffs];
    double delta[nCoeffs];
    double xin[nParameters];
    vector<double> x[nParameters];
    vector<string> paraNames;
    vector<double> y[nResponses];
    double xMin[nParameters];
    double xMax[nParameters];
    string names[nCoeffs];
    vector<int> deletedEntries;

    for(int i=0; i<nCoeffs; i++)
    {
        A[i] = 0.00;
    }
    if (argc > 8)
    {
        string coeffsFileName(argv[8]);
#       include "readCoeffs.H"
    }

#   include "readData.H"

    setNames(names, paraNames, nParameters);
    nEntries = x[0].size();

    for(int i=0; i<nParameters; i++)
    {
        xMin[i] = 1.0e+10;
	xMax[i] = -1.0e+10;
    }

    for(int j=0; j<nEntries; j++)
    {
	for(int i=0; i<nParameters; i++)
	{
	    xMin[i] = fmin(xMin[i], x[i][j]);
	    xMax[i] = fmax(xMax[i], x[i][j]);
	}
    }

    for(int i=0; i<nEntries; i++)
    {
        double yi = y[nFunction][i];
	//cout << "i = " << i << ", yi = " << yi << endl;
        func.push_back(yi);
    }
    vector<double> errs;
    errs.resize(func.size());

    int nCounter = 10;
    int iter = -1;
    double eMax = 1.0;

    while ( (iter<maxIter) && (eMax > errMax) && (urlx > 1.0e-30))
    {
        iter++;
	double err0 = err(x, A, nParameters, func, errs, xMin, xMax);

	nCounter++;
	int iMax = 0;
	for(int i=0; i<nCoeffs; i++)
	{
	    double A2[nCoeffs];
	    for(int j=0; j<nCoeffs; j++)
	    {
	        A2[j] = A[j];
	    }

	    A2[i] = A[i] + dx ;
	    double erp = err(x, A2, nParameters, func, errs, xMin, xMax);

	    A2[i] = A[i] - dx;
	    double erm = err(x, A2, nParameters, func, errs, xMin, xMax);

	    delta[i] = (erp - erm)/dx;

	    A[i] -= urlx*delta[i];
	    if (fabs(delta[i]) > fabs(delta[iMax]) )
	    {
		  iMax = i;
	    }
	}

	if (nCounter > output)
	{
	    nCounter = 1;

	    cout << iter << ", e = " << setw(7) << err0;
	    cout << ", d[" << setw(2) << iMax << "]=" << setw(7) << fabs(delta[iMax]);
	    for(int n=0; n<=nParameters; n++)
	    {
	      cout << ", " << setw(7) << A[n];
	    }
	    cout << ", " << setw(7) << A[iMax];
	    cout << ", " << urlx;
	    cout << ", " << errs.size();

	    cout << endl;
	    
	}
	//#       include "checkErrs.H"

    }
    
#   include "writeCoeffs.H"
#   include "orderCoeffs.H"

    for(int i=0; i<deletedEntries.size(); i++)
      {
	cout << deletedEntries[i]+1 << endl;
      }
    return 0;
}
