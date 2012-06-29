#include <iostream>
#include <math.h>
#include <string>
#include <fstream>
#include <cstdlib>
#include <list>
#include <vector>

using namespace std;

double value
(
    const double *x,
    const double *A,
    const int np
)
{
    double y = A[0];
    int index = 1;
    
    for(int i=0; i<np; i++)
    {
	y += A[index++]*x[i];
    }
    
    for(int i=0; i<np; i++)
    {
	for(int j=0; j<=i; j++)
	{
	    y += A[index++]*x[i]*x[j];
	}
    }

  return y;
}

int main(int argc, char *argv[])
{
    if (argc < 3)
    {
        cout << "missing input of parameter and result file. Exiting" << endl;
	cout << "Usage:\n\t " << argv[0] << " <nParameters> <nData>" << endl;
	exit(1);
    }

    int nParameters = atoi(argv[1]);
    int nData = atoi(argv[2]);
    int nCoeffs = 1 + 2*nParameters + (nParameters-1)*nParameters/2;

    cout << "nCoeffs = " << nCoeffs << endl;
    
    ifstream cFile("myCoeffs.dat");

    if (!cFile)
    {
        cout << "could not open myCoeffs.dat file " << endl;
        exit(1);
    }

    double A[nCoeffs];
    for(int i=0; i<nCoeffs; i++)
    {
      cFile >> A[i];
    }
    cFile.close();

    double x[nData][nParameters];

    ofstream paraFile("myData.dat");
    //ofstream resFile("myRes.dat");
    srand((unsigned)time(0)); 
    paraFile.precision(10);

    paraFile << "id";
    for(int i=0; i<nParameters; i++)
    {
        paraFile << "\t" << "var" << i;
    }
    
    paraFile << "\t" << "res" << endl;

    for(int i=0; i<nData; i++)
    {
        paraFile << i+1;
	double xin[nParameters];
        for(int j=0; j<nParameters; j++)
	{
	    x[i][j] = rand()/(double)RAND_MAX;
	    xin[j] = x[i][j];
	    paraFile << "\t" << x[i][j];
	}
	//paraFile << "\t0" << endl;


	double func = value(xin, A, nParameters);
	paraFile << i+1
		<< "\t" << func << endl;
    }
}

