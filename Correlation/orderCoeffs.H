int pointer[nCoeffs];
double oA[nCoeffs];

for(int i=0; i<nCoeffs; i++)
{
  oA[i] = A[i];
  pointer[i] = i;
}
for(int i=1; i<nCoeffs; i++)
{

    for(int j=1; j<nCoeffs; j++)
    {
      if (i != j)
	{
	  double magA = fabs(oA[j]);
	  if (magA < fabs(oA[i]))
	    {
	      int id = pointer[i];
	      pointer[i] = pointer[j];
	      pointer[j] = id;
	      double t = oA[i];
	      oA[i] = oA[j];
	      oA[j] = t;
	    }
	}
    }
}

ofstream file("sensParams.dat");
file.precision(10);

for(int i=1; i<nCoeffs; i++)
{
  file << "i = " << pointer[i] << ", A = " << oA[i] << ", " << names[pointer[i]] << endl;
}
