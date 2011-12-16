{

    ifstream file(fileName.c_str());

    if (!file)
    {
        cout << "could not open file " << fileName << endl;
        exit(1);
    }
    else
    {
        cout << "reading file " << fileName << "... ";
    }

    string input;
    double numInput;
    int nLine = 0;

    int index = 0;
    while(!file.eof())
    {
	nLine++;

        file >> input; // dummy, first entry is not used (id);

	if (!file.eof())
	{
	    // parameters
	    for(int i=0; i<nParameters; i++)
	    {
		if (nLine > 1)
		{
		    file >> numInput;
		    x[i].push_back(numInput);
		    //cout << numInput << " ";
		}
		else
		{
		    file >> input;
		    paraNames.push_back(input);
		}
	    }

	    file >> input; // dummy object funtion
	}
	//cout << endl;
    }
    file.close();
    cout << "Done!" << endl;
}
