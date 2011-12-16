{

    ofstream file("Acoeffs.dat");
    file.precision(10);
    for(int i=0; i<nCoeffs; i++)
    {
      file << A[i] << endl;
    }
    for(int i=0; i<nParameters; i++)
    {
      file << xMin[i] << endl;
    }

    for(int i=0; i<nParameters; i++)
    {
      file << xMax[i] << endl;
    }

    file.close();
}
