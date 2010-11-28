    Info<< "Reading controFoamDict properties.\n" << endl;

    IOdictionary controlFoamDict
    (
        IOobject
        (
            "controlFoamDict",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    List<word> methods(controlFoamDict.lookup("methods"));
    string caseDir(controlFoamDict.lookup("case"));

    PtrList<dictionary> methodsDictPtr(methods.size());
    PtrList<word> applicationsPtr(methods.size());
    PtrList<dictControl> actionsDictPtr(methods.size());

    forAll(methods, i)
    {
        methodsDictPtr.set
	(
	    i,
	    new dictionary(controlFoamDict.subDict(methods[i] + "Dict"))
        );
        actionsDictPtr.set(i, new dictControl(controlFoamDict, methods[i]));
    }

    forAll(actionsDictPtr, i)
    {
      //        applicationsPtr.set(i, new word(methodsDictPtr[i].lookup("application")));
        Info << "name = " << actionsDictPtr[i].name() << endl;
	Info << "options " << actionsDictPtr[i].options() << endl;
	Info << actionsDictPtr[i].dict() << endl;
    }
