double value
(
    const double *x,
    const double *A,
    const int np,
    const double *xMin,
    const double *xMax
)
{
    double y = A[0];
    double xin[np];
    int index = 1;
    
    for(int j=0; j<np; j++)
    {
	xin[j] = (x[j] - xMin[j])/(xMax[j]-xMin[j]);
	//xin[j] = x[j];
    }

    for(int i=0; i<np; i++)
    {
	y += A[index++]*xin[i];
    }
    
    for(int i=0; i<np; i++)
    {
        for(int j=0; j<i; j++) // include = for quadratic terms
	{
	    y += A[index++]*xin[i]*xin[j];
	}
    }

    return y;
}

double err
(
   const vector<double> *x,
   const double *A,
   const int nParam,
   const vector<double>& y,
   vector<double>& errs,
   const double *xMin,
   const double *xMax
)
{
    double e=0.0;
    double xin[nParam];
    
    int N = y.size();

    for(int i=0; i<N; i++)
    {
	for(int j=0; j<nParam; j++)
	{
	    xin[j] = x[j][i];
	}
	double f = value(xin, A, nParam, xMin, xMax);
	double e0 = f - y[i];
	double e02 = e0*e0;
	errs[i] = e02;
	e += e02;
    }
    //    double eAv = e/N;

    return sqrt(e);
}

void setNames
(
 string *names,
 const vector<string> paraNames,
 const int np
)
{
    int index = 1;

    for(int i=0; i<np; i++)
    {
	names[index] = paraNames[i];
	index++;
    }
    
    for(int i=0; i<np; i++)
    {
        for(int j=0; j<i; j++) // include = for quadratic terms
	{
	    names[index] = paraNames[i] + " * " + paraNames[j];
	    index++;
	}
    }
}
