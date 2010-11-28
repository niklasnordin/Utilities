IOdictionary snapEdgeDict
(
    IOobject
    (
        "snapEdgeDict",
         runTime.constant(),
         runTime,
         IOobject::MUST_READ,
         IOobject::NO_WRITE
    )
);

wordList snapPatches(snapEdgeDict.lookup("snapPatches"));
wordList snapZones(snapEdgeDict.lookup("snapZones"));

List<word> stlFileNames(snapEdgeDict.lookup("stlFileNames"));

scalar tolerance(readScalar(snapEdgeDict.lookup("tolerance")));
scalar relaxation(readScalar(snapEdgeDict.lookup("relaxation")));

scalar featureAngle(readScalar(snapEdgeDict.lookup("featureAngle")));
scalar excludeEdgeAngle(readScalar(snapEdgeDict.lookup("excludeEdgeAngle")));
scalar parallelAngle(readScalar(snapEdgeDict.lookup("parallelAngle")));

bool includeInterior(readBool(snapEdgeDict.lookup("includeInterior")));
label nIterations(readLabel(snapEdgeDict.lookup("nIterations")));

scalar overlapFactor(readScalar(snapEdgeDict.lookup("overlapFactor")));
scalar fitFactor(readScalar(snapEdgeDict.lookup("fitFactor")));

label nSmootherIterations(readLabel(snapEdgeDict.lookup("nSmootherIterations")));
scalar straightEdgePerturbation(readScalar(snapEdgeDict.lookup("straightEdgePerturbation")));

//scalar smootherRelaxation(readScalar(snapEdgeDict.lookup("smootherRelaxation")));
//scalar mapTolerance(readScalar(snapEdgeDict.lookup("mapTolerance")));
