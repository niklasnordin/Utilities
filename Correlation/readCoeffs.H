{

    ifstream file(coeffsFileName.c_str());

    if (!file)
    {
        cout << "could not open file " << coeffsFileName << endl;
        exit(1);
    }
    else
    {
        cout << "reading file " << coeffsFileName << "... ";
    }

    double numInput;
    for(int i=0; i<nCoeffs; i++)
    {
      file >> A[i];
    }
    file.close();
    cout << "Done!" << endl;
}
