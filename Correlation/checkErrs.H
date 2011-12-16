int eMax = 0;
int N = errs.size();
double errAv = 0.0;
for(int n=0; n<N; n++)
  {
    errAv += errs[n];
    if (errs[n] > errs[eMax])
      {
	eMax = n;
      }
  }
errAv /= N;

cout << "eMax = " << errs[eMax] << ", i = " << eMax << ", av = " << errAv << endl;
if (errs[eMax] > errDiff*errAv)
  {
    int lid = N-1;
    func[eMax] = func[lid];
    func.pop_back();
    errs.pop_back();
    for(int i=0; i<nParameters; i++)
      {
	x[i][eMax] = x[i][lid];
	x[i].pop_back();
      }

    for(int i=0; i<nResponses; i++)
      {
	y[i][eMax] = y[i][lid];
	y[i].pop_back();
      }
    deletedEntries.push_back(eMax);
  }

