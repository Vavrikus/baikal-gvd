#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

struct RunInfo
{
	int m_seasonID, m_clusterID, m_runID;
	double m_runTime; //days
	int m_Nentries, m_NFil, m_SixThreeFil, m_QFilChi2, m_TFil, m_TFilChi2, m_LikelihoodFit;
};

//operator overloading for printing
ostream& operator<<(ostream& stream, const RunInfo& rinfo)
{
	stream << "\nRUN INFO:\n";
	stream << "  seasonID: " << rinfo.m_seasonID << " clusterID: " << rinfo.m_clusterID;
	stream << " runID: " << rinfo.m_runID << "\n";
	stream << "  runTime: " << rinfo.m_runTime << " days";
	stream << "\n  Nentries:               " << rinfo.m_Nentries;
	stream << "\n  After NFilter:          " << rinfo.m_NFil;
	stream << "\n  After SixThreeFilter:   " << rinfo.m_SixThreeFil;
	stream << "\n  After QFilterChi2:      " << rinfo.m_QFilChi2;
	stream << "\n  After TFilter:          " << rinfo.m_TFil;
	stream << "\n  After TFilterChi2:      " << rinfo.m_TFilChi2;
	stream << "\n  After LikelihoodFitter: " << rinfo.m_LikelihoodFit << "\n";

	return stream;
}

void InputShort(ifstream& inf, string& variable)
{
	string line;

	inf >> line;
	inf >> variable;
}

void InputLong(ifstream& inf, string& variable)
{
	string line;

	getline(inf,line);
	inf >> line;
	inf >> line;
	inf >> variable;
}

//load info about runs into vector from log files
vector<RunInfo> parseRuns(string path)
{
	vector<RunInfo> dataOut;
	ifstream inf{path};

	while(inf)
	{
		string line;
		getline(inf,line);

		if(line == "RunInfo (Number of entries, RunTime [hours], runTime [days])")
		{
			//skip to line containing IDs
			getline(inf,line);

			string seasonID, clusterID, runID, Nentries, runTime;

			InputShort(inf,seasonID);
			InputShort(inf,clusterID);
			InputShort(inf,runID);
			InputShort(inf,Nentries);
			InputShort(inf,runTime);

			getline(inf,line);
			getline(inf,line);
			getline(inf,line);

			while(line != "*********************************************************************************")
				getline(inf,line);

			string AfterNFilter, AfterSixThreeFilter, AfterQFilterChi2, AfterTFilter,
				   AfterTFilterChi2, AfterLikelihoodFitter;			

			InputLong(inf,AfterNFilter);
			InputLong(inf,AfterSixThreeFilter);
			InputLong(inf,AfterQFilterChi2);
			InputLong(inf,AfterTFilter);
			InputLong(inf,AfterTFilterChi2);
			InputLong(inf,AfterLikelihoodFitter);

			int sID    = stoi(seasonID);
			int cID    = stoi(clusterID);
			int rID    = stoi(runID);
			double rT  = stod(runTime);
			int N      = stoi(Nentries);
			int ANF    = stoi(AfterNFilter);
			int ASTF   = stoi(AfterSixThreeFilter);
			int AQF2   = stoi(AfterQFilterChi2);
			int ATF    = stoi(AfterTFilter);
			int ATF2   = stoi(AfterTFilterChi2);
			int ALF    = stoi(AfterLikelihoodFitter);

			dataOut.emplace_back(RunInfo{sID,cID,rID,rT,N,ANF,ASTF,AQF2,ATF,ATF2,ALF});

			// cout << seasonID << " " << clusterID << " " << runID << " " << Nentries << " ";
			// cout << runTime << " " << AfterNFilter << " " << AfterSixThreeFilter;
			// cout << " " << AfterQFilterChi2 << " " << AfterTFilter << " " << AfterTFilterChi2;
			// cout << " " << AfterLikelihoodFitter << endl;
		}
	}

	for(RunInfo& ri : dataOut) cout << ri;

	return dataOut;
}

int readlogs()
{
	parseRuns("/home/vavrik/work/Baikal-GVD/cascade_flux/logs/programOutput_19_0.log");
	return 0;
}

int main()
{
	return readlogs();
}