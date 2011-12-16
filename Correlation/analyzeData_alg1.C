#include <iostream>
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

    const double dx = 1.0e-4;
    const double urlx = 1.0e-4;
    const int maxIter = 100000;
    const int output = 0;
    const double errMax = 1.0e-10;

    int nEntries;
    vector<double> func;

    cout.precision(5);
    if (argc < 6)
    {
        cout << "missing input of parameter and result file. Exiting" << endl;
	cout << "Usage:\n\t " << argv[0] << " <params> <results> <nParameters> <nFunctions> <targetFunction> <coeffs>" << endl;
	exit(1);
    }

    string fileName(argv[1]);  
    string resFileName(argv[2]);

    int nParameters = atoi(argv[3]);
    int nResponses = atoi(argv[4]);
    int nFunction = atoi(argv[5]) - 1;

    int nCoeffs = 1 + 2*nParameters + (nParameters-1)*nParameters/2;
    double A[nCoeffs];
    double delta[nCoeffs];
    double xin[nParameters];
    vector<double> x[nParameters];
    vector<string> paraNames;
    vector<double> y[nResponses];
    double xMin[nParameters];
    double xMax[nParameters];
    string names[nCoeffs];

    for(int i=0; i<nCoeffs; i++)
    {
        A[i] = 0.00;
    }
    if (argc > 6)
    {
        string coeffsFileName(argv[6]);
#       include "readCoeffs.H"
    }

#   include "readGraph.H"
#   include "readResult.H"
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
	cout << "i = " << i << ", yi = " << yi << endl;
        func.push_back(yi);
    }

    int nCounter = 0;
    int iter = 0;
    double eMax = 1.0;
    //for(int iter=0; iter<maxIter; iter++)
    while ( (iter<maxIter) && (eMax > errMax) )
    {
        iter++;
	double err0 = err(x, A, nParameters, func, xMin, xMax);

	nCounter++;

	for(int i=0; i<nCoeffs; i++)
	{
	    double A2[nCoeffs];
	    for(int j=0; j<nCoeffs; j++)
	    {
	        A2[j] = A[j];
	    }
	    //double scale = fabs(A[i]);

	    A2[i] = A[i] + dx;
	    double erp = err(x, A2, nParameters, func, xMin, xMax);

	    A2[i] = A[i] - dx;
	    double erm = err(x, A2, nParameters, func, xMin, xMax);

	    delta[i] = (erp - erm)/dx;
	    A[i] -= urlx*delta[i];
	}


	if (nCounter > output)
	{
	    nCounter = 0;

	    double deltaMax = 0.0;
	    for(int i=0; i<nCoeffs; i++)
	    {
		deltaMax = fmax(deltaMax, fabs(delta[i]));
	    }
	    eMax = deltaMax;

	    cout << "error = " << err0 << ", " 
                << A[0] << ", " << A[1] << ", "
                << A[2] << ", " << A[3]
		 << ", " << A[nCoeffs-3]
		 << ", " << A[nCoeffs-2]
		 <<  ", " << A[nCoeffs-1]
		 << ", " << deltaMax
                << endl;
	    
	}
    }
    
#   include "writeCoeffs.H"
#   include "orderCoeffs.H"

    return 0;
}
