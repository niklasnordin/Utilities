{

    ifstream file(resFileName.c_str());

    if (!file)
    {
        cout << "could not open file " << resFileName << endl;
        exit(1);
    }
    else
    {
        cout << "reading file " << resFileName << "... ";
    }

    string input;
    double numInput;
    int nLine = 0;

    while(!file.eof())
    {
	nLine++;

        file >> input; // dummy, first entry is not used (id);
	if (!file.eof())
	{
	    // parameters
	    for(int i=0; i<nResponses; i++)
	    {
	        file >> numInput;
		y[i].push_back(numInput);
	    }
	}

    }
    file.close();
    cout << "Done!" << endl;
}
