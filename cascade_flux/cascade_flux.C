#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THStack.h"
#include "TMath.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTimeStamp.h"
#include "TVector3.h"

#include <algorithm> //std::min
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <queue>

#define PROFILLING 0
#include "../Instrumentor.h"
#include "../transformations.h"

#define DEBUG() cout << "Current Line: " << __LINE__ << endl;

using namespace std;

// std::vector<double> stringXPositions = {5.22,52.13,57.54,25.17,-29.84,-53.6,-42.32,0};
// std::vector<double> stringYPositions = {62.32,37.15,-13.92,-52.01,-52.36,-7.49,42.74,0};

constexpr double xPos[40] = {-13.76,32.14,45.06,5.13,-45.03,-76.21,-59.85,-14.47,-195.19,-164.79,-180.08,-227.51,-276.24,-279.59,-248.17,-222.70,-270.25,-228.58,-220.89,-261.89,-309.86,-337.48,-319.74,-282.27,65.85,108.73,113.87,74.19,25.1,-2.48,16.08,58.37,-163.91,-119.26,-113.90,-152.28,-202.59,-230.83,-213.25,-170.30};
constexpr double yPos[40] = {-211.35,-235.88,-285.45,-325.83,-319.82,-281.63,-231.37,-270.17,-340.62,-384.09,-435.13,-450.13,-424.31,-372.59,-337.03,-391.09,-37.36,-65.26,-117.78,-153.57,-146.26,-101.43,-55.24,-96.82,-435.47,-462.39,-514.68,-549.90,-544.25,-500.53,-453,-491.97,-628.26,-656.49,-707.52,-744.24,-738.58,-694.13,-645.06,-685.35};

//{{{}},{{},{}},{{},{},{}},{{},{},{},{},{}},{{},{},{},{},{},{},{}}};//
const vector<vector<vector<int>>> ledMatrixRuns = {{{2,3,4,5,6,7,8,9,10,11,118,119,177,193,194,200,201,228,229,230,231,232,233,234,235,236,237,560,598}},{{},{}},{{7,117,412,429,443,459,474,490,505,520,548,564,579,595},{1,2,3,6,7,37,134,340,428,450,464,480,495,510,527,540,568,584,599,615,631,647,668},{35,36,117,120,131,151,412,429,443,459,474,489,504,519,520,547,575,591,607,623,644}},{{17,18,37,38,39,40,44,61,77,93,97,111,126,142,158,174,190,203,218,232,247,264,277,292,362,377,392,407,422,437,452,467,484,536,551,566,583,596,611,628,644,661,676,677,693},{8,41,54,56,60,61,77,92,107,123,138,154,169,184,201,215,231,245,260,276,306,375,391,406,421,436,451,466,481,498,553,571,586,603,616,631,648,663,679,694,709},{8,9,10,24,80,93,109,124,139,155,170,185,201,216,233,247,262,276,291,329,330,331,337,406,422,437,453,468,483,498,513,530,594,595,596,597,611,612,629,642,657,674,689,705,720,735},{13,23,36,51,67,82,100,116,131,146,162,179,193,208,222,237,251,268,283,350,367,384},{13,23,34,50,67,82,86,88,89,90,91,92,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,112,113,114,116,117,118,120,121,122,123,124,129,130,132,137,147,163,180,193,208,222,237,238,253,265,279,363,379}},{{3,19,32,42,51,52,62,71,82,92,102,122,145,156,165,180},{12,14,24,33,35,42,51,60,69,83,90,111,134,145,146,147,155,156,157,158,159,160,162,164},{9,13,14,17,132,143,153,164,165,167,169,172},{1,15,17,21,26,36,46,58,67,76,86,94,103,112,114},{2,12,17,19,23,24,26,36,44,55,63,73,82,89,98,106,117,131,143,151,160,166,168,175},{18,20,25,31,41,51,62,71,90,110,118,128,130,143,145,154,157,163,166,173,178,185,195,220,232,241,250,260,282,296,301,312,326,336,346,356,367,384,394,404,414,425,434,442,451},{7,10,12,16,17,22,30,40,49,58,67,76,84,93,102,105,113,115,129,131,143,144,149,152,159,165,169,174,177,207,219,228,237,245,254,264,277,281,290,301,312,322,332,342,359,369,380,389,398,407,417,426}}};

constexpr int unix1995 = 788918400;

struct Event
{
	int m_seasonID, m_clusterID, m_runID, m_eventID, m_nHits, m_nHitsAfterCaus;
	int m_nHitsAfterTFilter, m_nStringsAfterCaus, m_nStringsAfterTFilter, m_nTrackHits;
	double m_energy, m_theta, m_phi, m_mcEnergy, m_mcTheta, m_mcPhi;
	double m_energySigma, m_thetaSigma, m_phiSigma, m_directionSigma;
	double m_chi2AfterCaus, m_chi2AfterTFilter, m_cascTime, m_likelihood, m_likelihoodHitOnly;
	double m_qTotal;
	double m_rightAscension, m_declination;
	TVector3* m_position = new TVector3();
	TVector3* m_mcPosition = new TVector3();
	TTimeStamp* m_eventTime = new TTimeStamp();
	int m_coincidenceID = -1; //-1 if not in any coincidence

	//returns angular distance between this event and event in argument
	double angDist(const Event& ev) const
	{
	    TVector3 v1(0,0,1);
	    v1.SetTheta(TMath::Pi()/2.0+this->m_declination);
	    v1.SetPhi(this->m_rightAscension);

	    TVector3 v2(0,0,1);
	    v2.SetTheta(TMath::Pi()/2.0+ev.m_declination);
	    v2.SetPhi(ev.m_rightAscension);

	    return 180.0*v1.Angle(v2)/TMath::Pi();
	}

	//returns angular distance from given point on sky (input and output in degrees)
	double angDist(double ra, double dec)
	{
	    TVector3 v1(0,0,1);
	    v1.SetTheta(TMath::Pi()/2.0+this->m_declination);
	    v1.SetPhi(this->m_rightAscension);

	    if(ra < 0) ra += 360;

	    TVector3 v2(0,0,1);
	    v2.SetTheta(TMath::Pi()/2.0+degToRad(dec));
	    v2.SetPhi(degToRad(ra));

	    return radToDeg(v1.Angle(v2));
	}

	void LowTimeWarning()
	{
		//event before 01/01/2016 warning
		if(m_eventTime->GetSec() < 1451606400)
		{
			cout << "Event " << m_eventID << " has low eventTime: ";
			m_eventTime->Print();
			cout << " seasonID: " << m_seasonID << " clusterID: " << m_clusterID << " runID: " << m_runID << "\n";
		}
	}

	double Dist(const Event& ev) const
	{
		double dx2 = pow(((*(this->m_position))(0) - (*(ev.m_position))(0)),2);
		double dy2 = pow(((*(this->m_position))(1) - (*(ev.m_position))(1)),2);
		double dz2 = pow(((*(this->m_position))(2) - (*(ev.m_position))(2)),2);

		return TMath::Sqrt(dx2+dy2+dz2);
	}

	bool IsContained(double distFromCluster = 0) const
	{
		if (TMath::Sqrt(TMath::Power(m_position->X(),2)+TMath::Power(m_position->Y(),2)) < 60+distFromCluster && TMath::Abs(m_position->Z()) < 265+distFromCluster)
			return true;
		else
			return false;
	}

	bool IsUncontained(double near, double far) const
	{
		double horizontalDist = TMath::Sqrt(TMath::Power(m_position->X(),2)+TMath::Power(m_position->Y(),2));
		double verticalDist = TMath::Abs(m_position->Z());
		if ((horizontalDist < far && horizontalDist > near && verticalDist < 263) || (horizontalDist < far && verticalDist < 263+(far-60) && verticalDist > 263+(near-60)))
			return true;
		else
			return false;
	}

	bool IsLEDMatrixRun() const
	{
		bool isLEDMatrixRun = false;
		for (int i = 0; i < ledMatrixRuns[m_seasonID-2016][m_clusterID].size(); ++i)
		{
			if (m_runID == ledMatrixRuns[m_seasonID-2016][m_clusterID][i])
			{
				isLEDMatrixRun = true;
				break;
			}
		}
		return isLEDMatrixRun;
}

	static bool IsEarlier(const Event& ev1, const Event& ev2)
	{
		if(ev1.m_eventTime->GetSec() < ev2.m_eventTime->GetSec()) return true;
		else return false;
	}
};

// struct Coincidence
// {
// 	int m_id;
// 	vector<int> m_indexes; //indexes of events in sortedEvents

// 	Coincidence()
// 	{
// 		this->m_id = -1;
// 	}

// 	Coincidence(const Coincidence&) = delete;
// 	Coincidence& operator=(Coincidence other) = delete;

// 	// ~Coincidence() {}

// 	//returns angular distance of events number i and j in degrees
// 	double angDist(int i, int j) const
// 	{
// 	    return sortedEvents[m_indexes[i]].angDist(sortedEvents[m_indexes[j]]);
// 	}

// 	//adds event from sortedEvents
// 	void AddEvent(int index)
// 	{
// 		sortedEvents[index].m_coincidenceID = this->m_id;
// 		this->m_indexes.push_back(index);
// 	}
// };

struct RunInfo
{
	int m_seasonID, m_clusterID, m_runID;
	double m_runTime; //days
	int m_Nentries, m_NFil, m_SixThreeFil, m_QFilChi2, m_TFil, m_TFilChi2, m_LikelihoodFit;
	int m_CustomFil;

	static void InputShort(ifstream& inf, string& variable);
	static void InputLong(ifstream& inf, string& variable);
	static void parseRuns(vector<RunInfo>& dataOut, const string& path); 
};

void RunInfo::InputShort(ifstream& inf, string& variable)
{
	string line;

	inf >> line;
	inf >> variable;
}

void RunInfo::InputLong(ifstream& inf, string& variable)
{
	string line;

	getline(inf,line);
	inf >> line;
	inf >> line;
	inf >> variable;
}

//load info about runs into vector from log files
void RunInfo::parseRuns(vector<RunInfo>& dataOut, const string& path)
{
	ifstream inf{path};

	if(!inf) cout << "File: " << path << " was not found!" << endl;
	else     cout << "Parsing file " << path << endl;

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

			bool IsLEDrun = false;

			while(line != "*********************************************************************************")
			{
				if(line.find("Probably LED matrix run. Processing terminated!") != string::npos) 
					IsLEDrun = true;
				getline(inf,line);
			}

			string AfterNFilter, AfterSixThreeFilter, AfterQFilterChi2, AfterTFilter,
				   AfterTFilterChi2, AfterLikelihoodFitter;			

			InputLong(inf,AfterNFilter);
			InputLong(inf,AfterSixThreeFilter);
			InputLong(inf,AfterQFilterChi2);
			InputLong(inf,AfterTFilter);
			InputLong(inf,AfterTFilterChi2);
			InputLong(inf,AfterLikelihoodFitter);

			// cout << seasonID << " " << clusterID << " " << runID << " " << Nentries << " ";
			// cout << runTime << " " << AfterNFilter << " " << AfterSixThreeFilter;
			// cout << " " << AfterQFilterChi2 << " " << AfterTFilter << " " << AfterTFilterChi2;
			// cout << " " << AfterLikelihoodFitter << endl;

            if((AfterNFilter != "to") and !IsLEDrun)
			{
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

				dataOut.emplace_back(RunInfo{sID,cID,rID,rT,N,ANF,ASTF,AQF2,ATF,ATF2,ALF,0});
			}
		}
	}

	// for(const RunInfo& ri : dataOut) cout << ri;
}

class IDrawable
{
	typedef std::function<void(const Event&)> FillFn;
	typedef std::function<void()> DrawFn;

protected:
	FillFn fillfunc;
	DrawFn drawfunc;

public:
	virtual ~IDrawable() {}
	virtual void Fill(const Event&) = 0;
	virtual void Draw() = 0;
	virtual void Save() = 0;
	void SetFillFunc(FillFn f) {fillfunc = f;}
	void SetDrawFunc(DrawFn f) {drawfunc = f;}
};

template<typename DrawType>
class DrawMap : public IDrawable
{
public:
	map<TString,TCanvas*> canvasmap;
	map<TString,DrawType*> drawmap;

public:
	void Fill(const Event& e) override {fillfunc(e);}
	void Draw() override {drawfunc();}
	void Save() override {for(auto x : canvasmap) x.second->Write();}
};


template<typename DrawType>
class DrawSingle : public IDrawable
{
public:
	TCanvas* canvas;
	DrawType* drawsingle;

public:
	DrawSingle() {drawsingle = new DrawType();}

	void Fill(const Event& e) override {fillfunc(e);}
	void Draw() override {drawfunc();}
	void Save() override {canvas->Write();}
};

class EventLoop
{
	typedef std::function<bool(const Event&)> FilterFn;

private:
	int year,cluster,startID,endID,startSeason,endSeason;

	vector<FilterFn> filters;
	vector<IDrawable*> drawables;
	Event current_ev;

	map<TString,tuple<TGraph*,TGraph*>> flux_graphs;
	map<TString,TCanvas*> flux_canv;

public:
	vector<Event>   sortedEvents;
	vector<RunInfo> runs;
	vector<vector<const RunInfo*>> plotruns;
	TChain reconstructedCascades;
	TTree* filteredCascades;

private:
	void PrintProgress(int, int);
	bool CheckFilters();

public:
	EventLoop(int year, int cluster)
	{
		this->year    = year;
		this->cluster = cluster;

		startID = cluster!=-1?cluster:0;
		endID   = cluster!=-1?cluster+1:10;

		startSeason = year!=-1?year:16;
		endSeason   = year!=-1?year+1:20+1;
	}
	EventLoop(const EventLoop&) = delete;
	EventLoop& operator=(EventLoop other) = delete;

	void LoadReco(const char* env_p);
	void LoadRunLogs(string dir);
	void SetUpTTrees();
	void UseLEDfilter();
	void UseContainedFilter(double dist);
	void UseLikelihoodFilter(double max);
	void UseEnergyFilter(double min);
	void AddFilter(FilterFn f)     {filters.push_back(f);}
	void AddDrawable(IDrawable* d) {drawables.push_back(d);}
	void FillDrawables(const Event& e) {for(IDrawable* d : drawables) d->Fill(e);}
	void RunLoop();
	void DrawAll() 				 	   {for(IDrawable* d : drawables) d->Draw();}
	void SaveAll();
	void DrawFluxGraphs();
	int FindRunInfo(int seasonID, int clusterID, int runID);
	int FindRunInfo(Event* e) {return FindRunInfo(e->m_seasonID,e->m_clusterID,e->m_runID);}
	int FindRunInfo() {return FindRunInfo(&current_ev);}
};

void EventLoop::PrintProgress(int done, int all)
{
	if(floor(10*done/all) != floor(10*(done-1)/all)) 
		cout << "Filtering cascades progress: " << floor((100.0*done)/all) << "%" << endl;
	if(done == all-1) cout << "Filtering cascades progress: 100%" << endl;
}

bool EventLoop::CheckFilters()
{
	for(auto f : filters) if(!f(current_ev)) return false;
	return true;
}

void EventLoop::LoadReco(const char* env_p)
{
	reconstructedCascades.SetName("Tree/t_RecCasc");

	TString filesDir;

	for (int j = startSeason; j < endSeason; j++)
	{
		for (int i = startID; i < endID; ++i)
		{
			filesDir = Form("%s/exp20%d/cluster%d/",env_p,j,i);
			cout << filesDir << endl;

			auto dir = gSystem->OpenDirectory(filesDir.Data());
			while (auto f = gSystem->GetDirEntry(dir))
			{
			  	if (!strcmp(f,".") || !strcmp(f,"..")) continue;
			  	TString fullFilePath = filesDir + f + "/recCascResults.root";
			  	if (!gSystem->AccessPathName(fullFilePath))
			  	{
			  		// cout << f << endl;
			  		reconstructedCascades.Add(TString(filesDir) + f + "/recCascResults.root");
			  	}
			}
			
			gSystem->FreeDirectory(dir);
		}
	}
}

void EventLoop::LoadRunLogs(string dir)
{
	cout << "\n";

	for (int j = startSeason; j < endSeason; ++j)
	{
		for (int i = startID; i < endID; ++i)
		{
			string path = dir;

			path += to_string(j);
			path += "_";
			path += to_string(i);
			path += ".log";

			RunInfo::parseRuns(runs,path);
		}
	}
}

void EventLoop::SetUpTTrees()
{
	filteredCascades = new TTree("filteredCascades","Filtered Cascades");

	reconstructedCascades.SetBranchAddress("seasonID", 			   &current_ev.m_seasonID);
	reconstructedCascades.SetBranchAddress("clusterID", 		   &current_ev.m_clusterID);
	reconstructedCascades.SetBranchAddress("runID", 			   &current_ev.m_runID);
	reconstructedCascades.SetBranchAddress("eventID", 			   &current_ev.m_eventID);
	reconstructedCascades.SetBranchAddress("nHits", 			   &current_ev.m_nHits);
	reconstructedCascades.SetBranchAddress("nHitsAfterCaus", 	   &current_ev.m_nHitsAfterCaus);
	reconstructedCascades.SetBranchAddress("nStringsAfterCaus",    &current_ev.m_nStringsAfterCaus);
	reconstructedCascades.SetBranchAddress("chi2AfterCaus", 	   &current_ev.m_chi2AfterCaus);
	reconstructedCascades.SetBranchAddress("nHitsAfterTFilter",    &current_ev.m_nHitsAfterTFilter);
	reconstructedCascades.SetBranchAddress("nStringsAfterTFilter", &current_ev.m_nStringsAfterTFilter);
	reconstructedCascades.SetBranchAddress("chi2AfterTFilter", 	   &current_ev.m_chi2AfterTFilter);
	reconstructedCascades.SetBranchAddress("energy", 			   &current_ev.m_energy);
	reconstructedCascades.SetBranchAddress("energySigma", 		   &current_ev.m_energySigma);
	reconstructedCascades.SetBranchAddress("theta",		 		   &current_ev.m_theta);
	reconstructedCascades.SetBranchAddress("thetaSigma",		   &current_ev.m_thetaSigma);
	reconstructedCascades.SetBranchAddress("phi", 				   &current_ev.m_phi);
	reconstructedCascades.SetBranchAddress("phiSigma", 			   &current_ev.m_phiSigma);
	reconstructedCascades.SetBranchAddress("directionSigma", 	   &current_ev.m_directionSigma);
	reconstructedCascades.SetBranchAddress("declination",		   &current_ev.m_declination);
	reconstructedCascades.SetBranchAddress("rightAscension",	   &current_ev.m_rightAscension);
	reconstructedCascades.SetBranchAddress("position", 			   &current_ev.m_position);
	reconstructedCascades.SetBranchAddress("eventTime",			   &current_ev.m_eventTime);
	reconstructedCascades.SetBranchAddress("time", 				   &current_ev.m_cascTime);
	reconstructedCascades.SetBranchAddress("mcEnergy", 			   &current_ev.m_mcEnergy);
	reconstructedCascades.SetBranchAddress("mcTheta", 			   &current_ev.m_mcTheta);
	reconstructedCascades.SetBranchAddress("mcPhi", 			   &current_ev.m_mcPhi);
	reconstructedCascades.SetBranchAddress("mcPosition", 		   &current_ev.m_mcPosition);
	reconstructedCascades.SetBranchAddress("likelihood", 		   &current_ev.m_likelihood);
	reconstructedCascades.SetBranchAddress("likelihoodHitOnly",    &current_ev.m_likelihoodHitOnly);
	reconstructedCascades.SetBranchAddress("qTotal", 			   &current_ev.m_qTotal);
	reconstructedCascades.SetBranchAddress("nTrackHits", 		   &current_ev.m_nTrackHits);

	filteredCascades->Branch("seasonID", 				&current_ev.m_seasonID);
	filteredCascades->Branch("clusterID", 				&current_ev.m_clusterID);
	filteredCascades->Branch("runID", 					&current_ev.m_runID);
	filteredCascades->Branch("eventID", 				&current_ev.m_eventID);
	filteredCascades->Branch("nHits", 					&current_ev.m_nHits);
	filteredCascades->Branch("nHitsAfterCaus", 			&current_ev.m_nHitsAfterCaus);
	filteredCascades->Branch("nStringsAfterCaus", 		&current_ev.m_nStringsAfterCaus);
	filteredCascades->Branch("chi2AfterCaus", 			&current_ev.m_chi2AfterCaus);
	filteredCascades->Branch("nHitsAfterTFilter", 		&current_ev.m_nHitsAfterTFilter);
	filteredCascades->Branch("nStringsAfterTFilter",	&current_ev.m_nStringsAfterTFilter);
	filteredCascades->Branch("chi2AfterTFilter", 		&current_ev.m_chi2AfterTFilter);
	filteredCascades->Branch("energy", 					&current_ev.m_energy);
	filteredCascades->Branch("energySigma", 			&current_ev.m_energySigma);
	filteredCascades->Branch("theta", 					&current_ev.m_theta);
	filteredCascades->Branch("thetaSigma", 				&current_ev.m_thetaSigma);
	filteredCascades->Branch("phi", 					&current_ev.m_phi);
	filteredCascades->Branch("phiSigma", 				&current_ev.m_phiSigma);
	filteredCascades->Branch("directionSigma", 			&current_ev.m_directionSigma);
	filteredCascades->Branch("declination",				&current_ev.m_declination);
	filteredCascades->Branch("rightAscension",			&current_ev.m_rightAscension);
	filteredCascades->Branch("position", 				&current_ev.m_position);
	filteredCascades->Branch("eventTime","TTimeStamp",	&current_ev.m_eventTime);
	filteredCascades->Branch("time", 					&current_ev.m_cascTime);
	filteredCascades->Branch("mcEnergy", 				&current_ev.m_mcEnergy);
	filteredCascades->Branch("mcTheta", 				&current_ev.m_mcTheta);
	filteredCascades->Branch("mcPhi", 					&current_ev.m_mcPhi);
	filteredCascades->Branch("mcPosition", 				&current_ev.m_mcPosition);
	filteredCascades->Branch("likelihood", 				&current_ev.m_likelihood);
	filteredCascades->Branch("likelihoodHitOnly", 		&current_ev.m_likelihoodHitOnly);
	filteredCascades->Branch("qTotal", 					&current_ev.m_qTotal);
	filteredCascades->Branch("nTrackHits", 				&current_ev.m_nTrackHits);

	int nRecCasc = reconstructedCascades.GetEntries();

	cout << "\nnRecCasc: " << reconstructedCascades.GetEntries() << endl;

	sortedEvents.reserve(nRecCasc);
}

void EventLoop::UseLEDfilter()
{
	FilterFn LEDfilter = [](const Event& e){if(e.IsLEDMatrixRun()) return false; else return true;};
	filters.push_back(LEDfilter);
}

void EventLoop::UseContainedFilter(double dist)
{
	FilterFn Contained40Filter = [dist](const Event& e){if(e.IsContained(dist)) return true; else return false;};
	filters.push_back(Contained40Filter);
}

void EventLoop::UseLikelihoodFilter(double max)
{
	FilterFn LikelihoodFilter = [max](const Event& e){if(e.m_likelihoodHitOnly > max) return false; else return true;};
	filters.push_back(LikelihoodFilter);
}

void EventLoop::UseEnergyFilter(double min)
{
	FilterFn EnergyFilter = [min](const Event& e){if(e.m_energy < min) return false; else return true;};
	filters.push_back(EnergyFilter);
}

void EventLoop::RunLoop()
{
	for (int i = 0; i < reconstructedCascades.GetEntries(); ++i)
	{
		PrintProgress(i,reconstructedCascades.GetEntries());
		reconstructedCascades.GetEntry(i);
		if(!CheckFilters()) continue;

		filteredCascades->Fill();
		sortedEvents.push_back(current_ev);

		int run_index = FindRunInfo();
		if(run_index != -1) runs[run_index].m_CustomFil++;
		//else cout << "No info available for run " << runID << " from cluster " << clusterID << " from season " << seasonID << "." << endl;

		current_ev.LowTimeWarning();

		FillDrawables(current_ev);
	}

	cout << "\nnFilCasc: " << filteredCascades->GetEntries() << endl;
	sort(sortedEvents.begin(),sortedEvents.end(),Event::IsEarlier);
}

void EventLoop::SaveAll()
{
	TString outputFileName = Form("cascFlux_y%dc%d.root",year,cluster);
	TFile* outputFile = new TFile(outputFileName,"RECREATE");

	for(IDrawable* d : drawables) d->Save();
	for(auto x : flux_canv) x.second->Write();
}

//making graphs from logs (and CustomCut events)
void EventLoop::DrawFluxGraphs()
{
	for(const RunInfo& rinfo : runs)
	{
		//if(rinfo.m_LikelihoodFit/rinfo.m_runTime > 500) cout << rinfo;

		TString graph_key = to_string(rinfo.m_clusterID) + to_string(rinfo.m_seasonID);

		if(flux_graphs.find(graph_key) == flux_graphs.end())
		{
			TGraph* runs = new TGraph();
			TGraph* saverage = new TGraph();

			flux_graphs[graph_key] = make_tuple(runs,saverage);

			TString graph_title = Form("Cascade flux year %d cluster %d;RunID;Cascades per day",rinfo.m_seasonID,rinfo.m_clusterID);
			runs->SetTitle(graph_title);
			runs->SetName(Form("g_cascFlux_y%dc%d",rinfo.m_seasonID-2000,rinfo.m_clusterID));

			saverage->SetName(Form("g2_cascFlux_y%dc%d",rinfo.m_seasonID-2000,rinfo.m_clusterID));
			saverage->SetLineColor(kRed);
			saverage->SetLineWidth(4);

			plotruns.push_back(vector<const RunInfo*>());
		}

		get<0>(flux_graphs[graph_key])->SetPoint(get<0>(flux_graphs[graph_key])->GetN(),
										rinfo.m_runID,rinfo.m_CustomFil/rinfo.m_runTime);
		plotruns.back().push_back(&rinfo);
	}


	for(auto vec : plotruns) 
	{
		deque<tuple<double,double,double>> avg; //queue for sliding average (x,custom,runtime)

		for(auto rinfo : vec)
		{
			// cout << *rinfo;
			TString graph_key = to_string(rinfo->m_clusterID) + to_string(rinfo->m_seasonID);

			int n_avg = 10;

			avg.push_back(make_tuple(rinfo->m_runID,rinfo->m_CustomFil,rinfo->m_runTime));
			
			//cout << "almost there, deque size: " << avg.size() <<"\n";
			if(avg.size() == n_avg+1) 
			{
				avg.pop_front();

				double newX       = 0;
				double newCustom  = 0;
				double newRunTime = 0;

				for(auto x : avg)
				{
					// cout << "x: " << get<0>(x) << endl;
					newX 	   += get<0>(x);
					newCustom  += get<1>(x);
					newRunTime += get<2>(x);
				}

				newX = newX/n_avg;

				// cout << "adding point\n";
				// cout << "X: " << newX << "\n";
				get<1>(flux_graphs[graph_key])->SetPoint(get<1>(flux_graphs[graph_key])->GetN(),
												newX,newCustom/newRunTime);
			}
		}
	}

	for(auto const& x : flux_graphs)
	{
		TString season  = x.first(1,4);
		TString cluster = x.first(0,1);

		TString canvas_title = "Cascade flux year " + season + " cluster " + cluster;
		TString canvas_name  = "c_fluxGraph_y" + season(2,2) + "c" + cluster;

		flux_canv[x.first] = new TCanvas(canvas_name,canvas_title,800,600);
		get<0>(x.second)->Draw();
		get<1>(x.second)->Draw("same");
	}
}

//returns index of given run in vector of RunInfo
int EventLoop::FindRunInfo(int seasonID, int clusterID, int runID)
{
	for(int i = 0; i < runs.size(); i++)
	{
		if((runs[i].m_seasonID == seasonID) and (runs[i].m_clusterID == clusterID) and (runs[i].m_runID == runID))
			return i;
	}

	return -1;
}

//operator overloading for printing
std::ostream& operator<<(std::ostream& stream, const Event& ev)
{
	stream << "     =====================================  EVENT INFO  =====================================\n";
	stream << "     Time Stamp:\n          ";
	ev.m_eventTime->Print();

	stream << "\n     IDs:\n     ";
    stream << "     eventID: " << ev.m_eventID << ", coincidenceID: " << ev.m_coincidenceID;
    stream << ", seasonID: " << ev.m_seasonID << ", clusterID: " << ev.m_clusterID;
    stream << ", runID: " << ev.m_runID << "\n\n";

    stream << "     Reconstructed variables:\n     ";
    stream << "     energy = " << ev.m_energy << " TeV, sigma = " << ev.m_energySigma << " TeV\n     ";
    stream << "     theta = " << ev.m_theta/TMath::Pi()*180 << ", sigma = " << ev.m_thetaSigma/TMath::Pi()*180;
    stream << ", phi = " << ev.m_phi/TMath::Pi()*180 << ", sigma = " << ev.m_phiSigma/TMath::Pi()*180;
    stream << ", direction sigma = " << ev.m_directionSigma << "\n     ";
    stream << "     Position (XYZ): " << ev.m_position->X();
    stream << " " << ev.m_position->Y() << " " << ev.m_position->Z() << "\n     ";
    stream << "     right ascension = " << ev.m_rightAscension/TMath::Pi()*180;
    stream << ", declination = " << ev.m_declination/TMath::Pi()*180 << "\n\n";

    stream << "     Reconstruction parameters:\n     ";
    stream << "     nHits = " << ev.m_nHits << ", nTrackHits = " << ev.m_nTrackHits << ", qTotal = ";
    stream << ev.m_qTotal << "\n     ";
    stream << "     likelihood = " << ev.m_likelihood << ", likelihoodHitOnly = " << ev.m_likelihoodHitOnly;
    stream << "\n          cascTime = " << ev.m_cascTime << "\n   ";

    stream << "     After TFilter:\n     ";
    stream << "     nHits = " << ev.m_nHitsAfterTFilter <<  ", nStrings = " << ev.m_nStringsAfterTFilter;
    stream << ", chi2 = " << ev.m_chi2AfterTFilter << "\n   ";

    stream << "     After Caus:\n     ";
    stream << "     nHits = " << ev.m_nHitsAfterCaus << ", nStrings = " << ev.m_nStringsAfterCaus;
    stream << ", chi2 = " << ev.m_chi2AfterCaus << "\n\n";

    stream << "     MC data:\n     ";
    stream << "     energy = " << ev.m_mcEnergy << " TeV, theta = " << ev.m_mcTheta/TMath::Pi()*180; 
    stream << ", phi = " << ev.m_mcPhi/TMath::Pi()*180 << "\n     ";
    stream << "     Position (XYZ): " << ev.m_mcPosition->X();
    stream << " " << ev.m_mcPosition->Y() << " " << ev.m_mcPosition->Z() << endl;

    return stream;
}

// std::ostream& operator<<(std::ostream& stream, const Coincidence& c)
// {
// 	stream << "#######################################  COINCIDENCE  #######################################\n";
// 	stream << "ID: " << c.m_id << "\n";
// 	stream << "Number of events: " << c.m_indexes.size() << "\n";

// 	stream << "Time differences: ";

// 	for(int i = 1; i < c.m_indexes.size(); i++)
// 	{
// 		if(i!=1) cout << ", ";
// 		stream << sortedEvents[c.m_indexes[i]].m_eventTime->GetSec()-sortedEvents[c.m_indexes[i-1]].m_eventTime->GetSec(); 
// 	}

// 	stream << "\nMinimal angular distance: ";
// 	double minAngDist = c.angDist(0,1);

// 	for (int i = 0; i < c.m_indexes.size(); i++)
// 	{
// 		for(int j = i+1; j < c.m_indexes.size(); j++)
// 		{
// 			if(c.angDist(i,j) < minAngDist)
// 			{
// 				minAngDist = c.angDist(i,j);
// 			}
// 		}
// 	}

// 	stream << minAngDist << "\n\n";

// 	for(int ev_index : c.m_indexes) stream << sortedEvents[ev_index] << "\n" << endl;

//     return stream;
// }

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
	stream << "\n  After LikelihoodFitter: " << rinfo.m_LikelihoodFit;
	stream << "\n  After Custom Filter:    " << rinfo.m_CustomFil << "\n";

	return stream;
}

//returns number of leap years since 2016 (copied from transformations.h UTCtoUnix function)
int GetLeapYears(int season)
{
    int leapYears = std::floor((season - 2016) / 4);
    leapYears -= std::floor((season - 2000) / 100);
    leapYears += std::floor((season - 2000) / 400);

    return leapYears;
}

//returns unix (1995) time of 01/04/YYYY 00:00:00
int GetStartTime(int season)
{
	return 1459468800+((season-2016)*365+GetLeapYears(season))*86400-unix1995;
}

//returns unix (1995) time of 31/03/YYYY+1 23:59:59
int GetEndTime(int season)
{
	return 1491004799+((season-2016)*365+GetLeapYears(season))*86400-unix1995;
}

// //attempts to find a coincidence with given id and returns its index (-1 if unsuccessful)
// int FindCID(const int& id) 
// {
// 	int index = -1;
// 	for (int i = 0; i < coincidences.size(); ++i)
// 	{
// 		if(coincidences[i]->m_id == id) {index = i; break;}
// 	}

// 	if(index == -1) 
// 	{
// 		cout << "FindCID warning: Unable to find coincidence with id " << id << "." << endl;
// 		cout << "coincidences.size(): " << coincidences.size();
// 		cout << ", coincidences.back()->m_id: " << coincidences.back()->m_id << endl;
// 	}

// 	return index;
// }

// //outputs coincidence counts into console
// void WriteCStats()
// {
// 	vector<int> counts; //how many coincidences with given number of events, 0th index = 2

// 	for(shared_ptr<Coincidence> c : coincidences)
// 	{
// 		int cEvents = c->m_indexes.size();

// 		if(cEvents-1 > counts.size())
// 			for (int i = counts.size(); i < cEvents-1; ++i) counts.push_back(0);

// 		counts[cEvents-2]++;
// 	}

// 	cout << "#NoE   #NoC\n";
// 	cout << "=============\n";
// 	for (int i = 0; i < counts.size(); ++i)
// 	{
// 		random_coincidences->Fill(counts[i],i+2);
// 		if(i < 8)       cout << i+2 << "      " << counts[i] << "\n";
// 		else if(i < 98) cout << i+2 << "     " << counts[i] << "\n";
// 		else            cout << i+2 << "    " << counts[i] << "\n";
// 	}
// 	cout << "\nTotal coincidences: " << coincidences.size() << "\n\n";
// }

// //returns if events with indexes i,j are part coincidence with given time and angle differences
// //and if both events have greater energy (TeV) than given 
// bool IsTAEC(int i, int j, long int maxTimeDiff, double maxAngDist = 360, double minEnergy = 0)
// {
// 	bool timeOK   = sortedEvents[j].m_eventTime->GetSec() - sortedEvents[i].m_eventTime->GetSec() <= maxTimeDiff;
// 	bool angleOK  = sortedEvents[i].angDist(sortedEvents[j]) <= maxAngDist;
// 	bool energyOK = (sortedEvents[i].m_energy >= minEnergy) and (sortedEvents[j].m_energy >= minEnergy);
	
// 	return timeOK and angleOK and energyOK;
// }

// //writes coincidences with stats into console
// void WriteTAEC(long int maxTimeDiff, double maxAngDist = 360, double minEnergy = 0)
// {	
// 	cout << "\n\nCoincidences with maximal time difference " << maxTimeDiff;
// 	cout << " seconds, maximal distance " << maxAngDist << " degrees\n";
// 	cout << "and minimal energy " << minEnergy << " TeV:\n" << endl;
// 	for(auto const& c : coincidences) cout << *c;

// 	WriteCStats();
// }

// //writes warning if two or more cascades are separated by smaller than selected amount of time
// //and smaler than selected angle, saves coincidences into a vector
// template<typename... args>
// void FindCoincidences(bool(*IsCoin)(int, int, long int, args...), uint random_offset, long int maxdt, args... a)
// {
// 	int numOfCoincidences = 0; //number of created coincidences, may be bigger than actual count

// 	coincidences.clear();

// 	//reset coincidence IDs
// 	for(int i = 0; i < sortedEvents.size(); ++i) sortedEvents[i].m_coincidenceID = -1;

// 	if(sortedEvents.size() > 0)
// 	{
// 		for(int i = 0; i < sortedEvents.size()-1; ++i)
// 		{
// 			shared_ptr<Coincidence> c = make_shared<Coincidence>();

// 			bool IsCoincidence = false;
// 			int currentID = sortedEvents[i].m_coincidenceID;
// 			if(currentID != -1) continue;

// 			int startEvent;

// 			//if random time shift is applied, find first later event and set it as start 
// 			if(random_offset != 0)
// 			{
// 				//add random shift
// 				//cout << "old " << i << " " << sortedEvents[i].m_eventTime.GetSec() << "\n";
// 				sortedEvents[i].m_eventTime->SetSec(sortedEvents[i].m_eventTime->GetSec()+random_offset);
// 				//cout << "new " << i << " " << sortedEvents[i].m_eventTime.GetSec() << "\n";

// 				int lower_bound  = i;
// 				int higher_bound = sortedEvents.size()-1;

// 				while(lower_bound != higher_bound-1)
// 				{
// 					int mid = floor((higher_bound+lower_bound)/2.0);
// 					//cout << lower_bound << " " << mid << " " << higher_bound << endl;
// 					if(sortedEvents[i].m_eventTime->GetSec() > sortedEvents[mid].m_eventTime->GetSec())
// 						lower_bound = mid;
// 					else higher_bound = mid;
// 				}

// 				//if should not happen for good shift
// 				if(higher_bound == i) 
// 				{
// 					startEvent = i+1;
// 					cout << "WARNING: Random shift insufficient.";
// 				}

// 				else startEvent = higher_bound;
// 				//cout << "startEvent " << startEvent << " " << sortedEvents[startEvent].m_eventTime.GetSec() << "\n";
// 			}

// 			else startEvent = i+1;

// 			//searching events coinciding with event i
// 			for (int j = startEvent; j < sortedEvents.size(); ++j)
// 			{
// 				if(sortedEvents[j].m_eventTime->GetSec()-maxdt > sortedEvents[i].m_eventTime->	GetSec())
// 					break;

// 				if(IsCoin(i,j,maxdt,a...) and (sortedEvents[j].m_coincidenceID == -1))
// 				{	
// 					if(!IsCoincidence)
// 					{
// 						IsCoincidence = true;

// 						c->m_indexes.clear();
// 						c->m_id = numOfCoincidences;
// 						c->AddEvent(i);

// 						numOfCoincidences++;
// 					}

// 					c->AddEvent(j);					
// 				}
// 			}

// 			//adding new coincidence
// 			if(IsCoincidence and (currentID == -1)) coincidences.push_back(c);

// 			if(random_offset != 0) //remove random shift
// 			{
// 				//cout << "old2 " << i << " " << sortedEvents[i].m_eventTime.GetSec() << "\n";				
// 				sortedEvents[i].m_eventTime->SetSec(sortedEvents[i].m_eventTime->GetSec()-random_offset);
// 				//cout << "new2 " << i << " " << sortedEvents[i].m_eventTime.GetSec() << "\n";
// 			}
// 		}
// 	}


// 	for(shared_ptr<Coincidence> c : coincidences) sort(c->m_indexes.begin(),c->m_indexes.end());

// 	if((void*)IsCoin == (void*)IsTAEC) WriteTAEC(maxdt,a...);
// }

// //returns if events with indexes i,j satisfy time, run and position criterions
// bool IsTRPC(int i, int j, long int maxTimeDiff = 3600, double maxDist = 5)
// {
// 	bool timeOK = sortedEvents[j].m_eventTime->GetSec() - sortedEvents[i].m_eventTime->GetSec() <= maxTimeDiff;
// 	bool posOK  = sortedEvents[i].Dist(sortedEvents[j]) <= maxDist;
// 	bool yearOK	= sortedEvents[i].m_seasonID  == sortedEvents[j].m_seasonID;
// 	bool clusOK = sortedEvents[i].m_clusterID == sortedEvents[j].m_clusterID;
// 	bool runOK  = sortedEvents[i].m_runID	  == sortedEvents[j].m_runID;

// 	return timeOK and posOK and yearOK and clusOK and runOK;
// }

// void WarnLEDMatrixRun(int minCoinSize = 3, long int maxTimeDiff = 3600, double maxDist = 5)
// {
// 	PROFILE_FUNCTION();

// 	cout << "\nPossible LED matrix runs detected:\n  coincidences with ";
// 	cout << minCoinSize << " or more events\n  position difference max ";
// 	cout << maxDist <<" meters\n  time difference max " << maxTimeDiff << " seconds\n";
// 	cout << "======================================\n";

// 	bool noLEDRunsDetected = true;

// 	FindCoincidences(IsTRPC,0,maxTimeDiff,maxDist);

// 	for(shared_ptr<Coincidence> c : coincidences)
// 	{
// 		for(int season = 2016; season < 2021; season++)
// 		{
// 			for(int cluster = 0; cluster < 10; cluster++)
// 			{
// 				if(sortedEvents[c->m_indexes[0]].m_seasonID == season && sortedEvents[c->m_indexes[0]].m_clusterID == cluster && c->m_indexes.size() >= minCoinSize)
// 				{
// 					noLEDRunsDetected = false;
// 					cout << "seasonID: "   << sortedEvents[c->m_indexes[0]].m_seasonID;
// 					cout << " clusterID: " << sortedEvents[c->m_indexes[0]].m_clusterID;
// 					cout << " runID: "     << sortedEvents[c->m_indexes[0]].m_runID << "\n";
// 				}

// 			}
// 		}
// 	}

// 	if(noLEDRunsDetected) cout << "No runs detected." << endl;
// 	cout << endl;
// }

// //returns index of given run in vector of RunInfo
// int FindRunInfo(int seasonID, int clusterID, int runID)
// {
// 	for(int i = 0; i < runs.size(); i++)
// 	{
// 		if((runs[i].m_seasonID == seasonID) and (runs[i].m_clusterID == clusterID) and (runs[i].m_runID == runID))
// 			return i;
// 	}

// 	return -1;
// }

// void RandomCoincidences(long int maxTimeDiff,double maxAngDist = 20)
// {
// 	int iterations = 10000;

// 	for (int i = 0; i < iterations; ++i)
// 	{
// 		FindCoincidences(IsTAEC,(i+1)*20,maxTimeDiff,maxAngDist,20.0);
// 	}

// 	random_coincidences->Fill(0.0,2.1,iterations-random_coincidences->GetEntries());
// 	random_coincidences->GetXaxis()->SetTitle("#NoC");
// 	random_coincidences->GetYaxis()->SetTitle("#NoE");
// 	random_coincidences->Draw("Lego2");
// }

int cascade_flux(int val = 0, int year = -1, int cluster = -1)
{
#if PROFILLING
	Instrumentor::Get().BeginSession("Session Name");
#endif
{
	PROFILE_SCOPE("MAIN");

	int LCut;
	long int maxTimeDiff;
	double energyCut = 200;

	cout << "Contained + 40 cut is applied.\n";
	cout << "Apply likelihood cut? [1/0]\n";
	cin >> LCut;

	bool DoLikelihoodCut = LCut == 1;
	if(DoLikelihoodCut) cout << "Likelihood (Lho <= 1.5) cut is applied.\n";

	cout << "\nMaximal time difference for coincidences in seconds: \n";
	cin >> maxTimeDiff;
	cout << "Maximal time difference set to " << maxTimeDiff << " seconds.\n";

	typedef std::function<void(const Event&)> FillFn;
	typedef std::function<void()> DrawFn;

	//make aitoff map for energyCut TeV+
	DrawSingle<TGraph>* aitoff = new DrawSingle<TGraph>();

	FillFn fillAitoff = [&energyCut,&aitoff](const Event& e)
	{
		if(e.m_energy >= energyCut)
		{
			double ra = radToDeg(e.m_rightAscension);
			if(ra > 180) ra -= 360;
			XY pos = toAitoff(ra,radToDeg(e.m_declination));
			// cout << "ra: " << radToDeg(rightAscension) << " dec: " << radToDeg(declination) << endl;
			// cout << "aitoff x: " << pos.x << " aitoff y: " << pos.y << endl;
			aitoff->drawsingle->SetPoint(aitoff->drawsingle->GetN(),pos.x,pos.y);
		}
	};
	aitoff->SetFillFunc(fillAitoff);

	DrawFn drawAitoff = [&energyCut,&aitoff]()
	{
		aitoff->canvas = new TCanvas("c_aitoff", "", 1000, 500);
		drawmap(Form("Cascades with energy > %f TeV;Right ascension;Declination", energyCut));
		aitoff->drawsingle->SetMarkerColor(kBlue);
		aitoff->drawsingle->SetMarkerSize(1);
		aitoff->drawsingle->SetMarkerStyle(4);
		aitoff->drawsingle->Draw("*");
		drawLabels();
	};
	aitoff->SetDrawFunc(drawAitoff);

	//making flux histograms
	DrawMap<TH1F>* flux_hist = new DrawMap<TH1F>();

	FillFn fillFlux = [&flux_hist](const Event& e)
	{
		//key for histogram identification
		TString hist_key = to_string(e.m_clusterID)+to_string(e.m_seasonID);

		//if histogram with given key does not exist, create one, fill histogram
		if(flux_hist->drawmap.find(hist_key)==flux_hist->drawmap.end())
		{
			TString hist_name = Form("h_cascFlux_y%dc%d",e.m_seasonID-2000,e.m_clusterID);
			flux_hist->drawmap[hist_key] = new TH1F(hist_name,hist_name(11,5)+"; Month; NoE [#]",50,GetStartTime(2016),GetEndTime(2016));

			flux_hist->drawmap[hist_key]->GetXaxis()->SetTimeDisplay(1);
			flux_hist->drawmap[hist_key]->GetXaxis()->SetTimeFormat("%m");//("%m/%Y");

			int color = e.m_seasonID-2014+e.m_clusterID;
			if(color > 9) color += 30;

			flux_hist->drawmap[hist_key]->SetLineColor(color);	
		}

		flux_hist->drawmap[hist_key]->Fill(e.m_eventTime->GetSec()-unix1995-GetStartTime(e.m_seasonID)+GetStartTime(2016)); //1970 unix to 1995 unix
	};
	flux_hist->SetFillFunc(fillFlux);

	DrawFn drawFlux = [&flux_hist]()
	{
		map<TString,THStack*> flux_stack;
		for(auto const& x : flux_hist->drawmap)
		{
			// flux_canv[x.first] = new TCanvas(x.first,"CascadeFlux",800,600);
			// x.second->Draw();

			TString season = x.first(1,4);
			TString cluster = x.first(0,1);

			//if THStack with given key does not exist, create one, fill THStack
			if(flux_stack.find(season)!=flux_stack.end()) flux_stack[season]->Add(x.second);
			else
			{
				TString stack_name = "hs_cascFlux_y"+season(2,2);
				flux_stack[season] = new THStack(stack_name,"Cascade flux over time; Month; NoE [#] per 7.3 days");

				flux_stack[season]->Add(x.second);
			}

			if(flux_stack.find(cluster)!=flux_stack.end()) flux_stack[cluster]->Add(x.second);
			else
			{
				TString stack_name = "hs_cascFlux_c"+cluster;
				flux_stack[cluster] = new THStack(stack_name,"Cascade flux over time; Month; NoE [#] per 7.3 days");

				flux_stack[cluster]->Add(x.second);
			}		
		}

		for(auto const& x : flux_stack)
		{
			flux_hist->canvasmap[x.first] = new TCanvas("c_cascFlux_" + x.first,"CascadeFlux",800,600);
			x.second->Draw("nostack");

			x.second->GetXaxis()->SetTimeDisplay(1);
			x.second->GetXaxis()->SetTimeFormat("%m");
			x.second->Draw("nostack");
	  		gPad->BuildLegend(0.75,0.75,0.95,0.95,"");	
		}
	};
	flux_hist->SetDrawFunc(drawFlux);

	EventLoop* eloop = new EventLoop(year,cluster);

	const char* data_path = val==1?"/home/vavrik/bajkal/recoCascades/v1.2":"/home/vavrik/work/data";//"/media/vavrik/Alpha/BaikalData/dataGidra_v1.3";//
	eloop->LoadReco(data_path);

	string logs_path;
	if(val == 0) logs_path = "/home/vavrik/work/Baikal-GVD/cascade_flux/logs/programOutput_";
	if(val == 1) logs_path = "/home/vavrik/storage/casc_flux/logs/programOutput_";
	eloop->LoadRunLogs(logs_path);
	eloop->SetUpTTrees();

	eloop->UseLEDfilter();
	eloop->UseContainedFilter(40);
	if(DoLikelihoodCut) eloop->UseLikelihoodFilter(1.5);
	eloop->UseEnergyFilter(energyCut);

	eloop->AddDrawable(aitoff);
	eloop->AddDrawable(flux_hist);
	eloop->RunLoop();
	eloop->DrawAll();
	eloop->DrawFluxGraphs();
	eloop->SaveAll();

// 	// FindCoincidences(IsTAEC,0,maxTimeDiff,360.0,20.0);
// 	// FindCoincidences(IsTAEC,0,maxTimeDiff,20.0,20.0);
// 	// FindCoincidences(IsTAEC,0,maxTimeDiff,10.0,20.0);

// 	//RandomCoincidences(maxTimeDiff);

// 	WarnLEDMatrixRun();

// 	gStyle->SetOptStat(111111);

// 	//DrawResults(val);
// 	SaveResults(year,cluster);

// 	//fixed window events (turn on energyCut)
// 	double window = 6; //angle window in degrees
// 	int step = 3;

// 	TString title = Form("Events E > %f TeV per %f degree window;Right ascension;Declination",energyCut,window);
// 	TH2F* h_density = new TH2F("h_density",title,360/step,-180,180,180/step+1,-90,91);
// 	TH1F* h_density_stats = new TH1F("h_density_stats","Bins with given number of events;#NoE;#NoB",16,0,16);

// 	for(int ra = -180; ra < 181; ra += step)
// 	{
// 		for(int dec = -90; dec < 91; dec += step)
// 		{
// 			int noe = 0;

// 			for(auto ev : sortedEvents)
// 			{
// 				if(ev.angDist(ra,dec) < window)
// 				{
// 					h_density->Fill(ra,dec);
// 					noe++;
// 				}
// 			}
// 			h_density_stats->Fill(noe);
// 		}
// 	}

// 	TCanvas* c2 = new TCanvas("c2", "", 1000, 500);
// 	h_density->Draw("Lego2");

// 	TCanvas* c3 = new TCanvas("c3", "", 1000, 500);
// 	h_density_stats->Draw();


// 	cout << "nProcessedEvents: " << nProcessedEvents << endl;
// 	TString outputFileName = Form("filteredCascades_y%dc%d.root",year,cluster);
// 	TFile *newFile = new TFile(outputFileName,"recreate");
// 	filteredCascades->Write();
// 	newFile->Close();

// 	//for(auto ev : sortedEvents) cout << ev << "\n";
}

// #if PROFILLING
// 	Instrumentor::Get().EndSession();
// #endif

	return 0;
}

int main(int argc, char** argv) 
{
	int val, year, cluster;

	if(argc < 4) cluster = -1;
	else cluster = stoi(argv[3]);

	if(argc < 3) year = -1;
	else year = stoi(argv[2]);
	
	if(argc < 2) val = 0;
	else val = stoi(argv[1]);

	return cascade_flux(val,year,cluster);
}