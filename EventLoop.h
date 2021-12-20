#pragma once

#include <functional>
#include <iostream>
#include <map>
#include <vector>

#include "TCanvas.h"
#include "TChain.h"
#include "TTimeStamp.h"
#include "TVector3.h"

#ifdef NEW_CASC_STRUCTURE
#include "BRecoCascade.h"
#include "BJointHeader.h"
#include "BRunInfo.h"
#endif

#include "transformations.h"

using namespace std;

// std::vector<double> stringXPositions = {5.22,52.13,57.54,25.17,-29.84,-53.6,-42.32,0};
// std::vector<double> stringYPositions = {62.32,37.15,-13.92,-52.01,-52.36,-7.49,42.74,0};

constexpr double xPos[40] = {-13.76,32.14,45.06,5.13,-45.03,-76.21,-59.85,-14.47,-195.19,-164.79,-180.08,-227.51,-276.24,-279.59,-248.17,-222.70,-270.25,-228.58,-220.89,-261.89,-309.86,-337.48,-319.74,-282.27,65.85,108.73,113.87,74.19,25.1,-2.48,16.08,58.37,-163.91,-119.26,-113.90,-152.28,-202.59,-230.83,-213.25,-170.30};
constexpr double yPos[40] = {-211.35,-235.88,-285.45,-325.83,-319.82,-281.63,-231.37,-270.17,-340.62,-384.09,-435.13,-450.13,-424.31,-372.59,-337.03,-391.09,-37.36,-65.26,-117.78,-153.57,-146.26,-101.43,-55.24,-96.82,-435.47,-462.39,-514.68,-549.90,-544.25,-500.53,-453,-491.97,-628.26,-656.49,-707.52,-744.24,-738.58,-694.13,-645.06,-685.35};

//{{{}},{{},{}},{{},{},{}},{{},{},{},{},{}},{{},{},{},{},{},{},{}},{{},{},{},{},{},{},{},{},{}}};//
const vector<vector<vector<int>>> ledMatrixRuns = {{{2,3,4,5,6,7,8,9,10,11,118,119,177,193,194,200,201,228,229,230,231,232,233,234,235,236,237,560,598}},{{},{}},{{7,117,412,429,443,459,474,490,505,520,548,564,579,595},{1,2,3,6,7,37,134,340,428,450,464,480,495,510,527,540,568,584,599,615,631,647,668},{35,36,117,120,131,151,412,429,443,459,474,489,504,519,520,547,575,591,607,623,644}},{{17,18,37,38,39,40,44,61,77,93,97,111,126,142,158,174,190,203,218,232,247,264,277,292,362,377,392,407,422,437,452,467,484,536,551,566,583,596,611,628,644,661,676,677,693},{8,41,54,56,60,61,77,92,107,123,138,154,169,184,201,215,231,245,260,276,306,375,391,406,421,436,451,466,481,498,553,571,586,603,616,631,648,663,679,694,709},{8,9,10,24,80,93,109,124,139,155,170,185,201,216,233,247,262,276,291,329,330,331,337,406,422,437,453,468,483,498,513,530,594,595,596,597,611,612,629,642,657,674,689,705,720,735},{13,23,36,51,67,82,100,116,131,146,162,179,193,208,222,237,251,268,283,350,367,384},{13,23,34,50,67,82,86,88,89,90,91,92,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,112,113,114,116,117,118,120,121,122,123,124,129,130,132,137,147,163,180,193,208,222,237,238,253,265,279,363,379}},{{3,19,32,42,51,52,62,71,82,92,102,122,145,156,165,180},{12,14,24,33,35,42,51,60,69,83,90,111,134,145,146,147,155,156,157,158,159,160,162,164},{9,13,14,17,132,143,153,164,165,167,169,172},{1,15,17,21,26,36,46,58,67,76,86,94,103,112,114},{2,12,17,19,23,24,26,36,44,55,63,73,82,89,98,106,117,131,143,151,160,166,168,175},{18,20,25,31,41,51,62,71,90,110,118,128,130,143,145,154,157,163,166,173,178,185,195,220,232,241,250,260,282,296,301,312,326,336,346,356,367,384,394,404,414,425,434,442,451},{7,10,12,16,17,22,30,40,49,58,67,76,84,93,102,105,113,115,129,131,143,144,149,152,159,165,169,174,177,207,219,228,237,245,254,264,277,281,290,301,312,322,332,342,359,369,380,389,398,407,417,426}},{{167,178,182,204,207,252,259,267,275,281,284,292,300},{158,160,168,171,177,195,244,255,263,271,280,288},{111,118,120,127,129,138,140,147,149,156,158,166,168,191,193,236,240,247},{110,118,120,127,129,137,139,146,148,156,157,165,172,191,193,233,238,247,255,263,269,272,280,288},{161,162,177,193,196,236,248,256,264,270,273,281},{170,171,179,186,203,206,247,252,259,268,276,282,285,293},{167,168,176,184,201,203,246,250,258,266,275,281,284,292,300},{164,165,172,175,182,200,202,245,252,259,269,279,285,288,296,304},{}}};

constexpr int unix1995 = 788918400;

struct BasicEvent
{
	double m_energy, m_theta, m_phi;
	double m_rightAscension, m_declination;
	TTimeStamp m_eventTime;

#ifndef NEW_CASC_STRUCTURE
	//returns angular distance between this event and event in argument
	double angDist(const BasicEvent& ev) const
	{
	    TVector3 v1(0,0,1);
	    v1.SetTheta(TMath::Pi()/2.0+this->m_declination);
	    v1.SetPhi(this->m_rightAscension);

	    TVector3 v2(0,0,1);
	    v2.SetTheta(TMath::Pi()/2.0+ev.m_declination);
	    v2.SetPhi(ev.m_rightAscension);

	    return radToDeg(v1.Angle(v2));
	}

	//returns angular distance from given point on sky (input and output in degrees)
	double angDist(double ra, double dec) const
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
#else
	//returns angular distance between this event and event in argument
	double angDist(const BasicEvent& ev) const
	{
	    TVector3 v1(0,0,1);
	    v1.SetTheta(TMath::Pi()/2.0+degToRad(this->m_declination));
	    v1.SetPhi(degToRad(this->m_rightAscension));

	    TVector3 v2(0,0,1);
	    v2.SetTheta(TMath::Pi()/2.0+degToRad(ev.m_declination));
	    v2.SetPhi(degToRad(ev.m_rightAscension));

	    return radToDeg(v1.Angle(v2));
	}

	//returns angular distance from given point on sky (input and output in degrees)
	double angDist(double ra, double dec) const
	{
	    TVector3 v1(0,0,1);
	    v1.SetTheta(TMath::Pi()/2.0+degToRad(this->m_declination));
	    v1.SetPhi(degToRad(this->m_rightAscension));

	    if(ra < 0) ra += 360;

	    TVector3 v2(0,0,1);
	    v2.SetTheta(TMath::Pi()/2.0+degToRad(dec));
	    v2.SetPhi(degToRad(ra));

	    return radToDeg(v1.Angle(v2));
	}
#endif //NEW_CASCADE_STRUCTURE

	//trasform from horizontal to equatorial coordinates
	void computeRaDec(const double& latDet = latB, const double& lonDet = lonB)
	{
		double alt = -TMath::Pi()/2+m_theta;
		double azm = 3*TMath::Pi()/2-m_phi;

	    //transform using equations from http://star-www.st-and.ac.uk/~fv/webnotes/chapter7.htm
	    double sinDec = sin(alt) * sin(latDet) + cos(alt) * cos(azm) * cos(latDet);

	    //making sure asin does not return NaN
	    if(sinDec > 1) sinDec = 1;
	    else if (sinDec < -1) sinDec = -1;

	    m_declination = asin(sinDec);

	    double cosLHA = (sin(alt) - sin(m_declination) * sin(latDet)) / (cos(m_declination) * cos(latDet));

	    //making sure acos does not return NaN
	    if(cosLHA > 1) cosLHA = 1;
	    else if (cosLHA < -1) cosLHA = -1;

	    double localHourAngle = acos(cosLHA);

	    //expressing localHourAngle in degrees
	    localHourAngle = radToDeg(localHourAngle);

	    if (sin(azm) > 0) localHourAngle = 360.0 - localHourAngle;

	    double localSiderealTime = LST(m_eventTime, lonDet);

	    //calculating right ascension
	    double rightAscension = localSiderealTime - localHourAngle;
	    if (rightAscension < 0) rightAscension += 360.0;

	    m_rightAscension = degToRad(rightAscension);
	}

	static bool IsEarlier(const BasicEvent& ev1, const BasicEvent& ev2)
	{
		if(ev1.m_eventTime.GetSec() < ev2.m_eventTime.GetSec()) return true;
		else return false;
	}
};

struct Event : public BasicEvent
{
	int m_seasonID, m_clusterID, m_runID, m_eventID, m_nHits, m_nHitsAfterCaus;
	int m_nHitsAfterTFilter, m_nStringsAfterCaus, m_nStringsAfterTFilter, m_nTrackHits;
	double m_mcEnergy, m_mcTheta, m_mcPhi;
	double m_energySigma, m_thetaSigma, m_phiSigma, m_directionSigma;
	double m_chi2AfterCaus, m_chi2AfterTFilter, m_cascTime, m_likelihood, m_likelihoodHitOnly;
	double m_qTotal;
	TVector3 m_position;
	TVector3 m_mcPosition;
	int m_coincidenceID = -1; //-1 if not in any coincidence
	double m_distanceCS;

	Event() = default;
#ifdef NEW_CASC_STRUCTURE	
	Event(BRecoCascade*,BJointHeader*);
#endif

	void LowTimeWarning()
	{
		//event before 01/01/2016 warning
		if(m_eventTime.GetSec() < 1451606400)
		{
			cout << "Event " << m_eventID << " has low eventTime: ";
			m_eventTime.Print();
			cout << " seasonID: " << m_seasonID << " clusterID: " << m_clusterID << " runID: " << m_runID << "\n";
		}
	}

	double Dist(const Event& ev) const
	{
		double dx2 = pow(((this->m_position)(0) - (ev.m_position)(0)),2);
		double dy2 = pow(((this->m_position)(1) - (ev.m_position)(1)),2);
		double dz2 = pow(((this->m_position)(2) - (ev.m_position)(2)),2);

		return TMath::Sqrt(dx2+dy2+dz2);
	}

#ifndef NEW_CASC_STRUCTURE
	bool IsContained(double distFromCluster = 0) const
	{
		if (TMath::Sqrt(TMath::Power(m_position.X(),2)+TMath::Power(m_position.Y(),2)) < 60+distFromCluster && TMath::Abs(m_position.Z()) < 265+distFromCluster)
			return true;
		else
			return false;
	}
#else
	bool IsContained(double distFromCluster = 0) const
	{
		if (m_distanceCS < 60+distFromCluster && TMath::Abs(m_position.Z()) < 625+distFromCluster && TMath::Abs(m_position.Z()) > 100-distFromCluster)
			return true;
		else
			return false;
	}
#endif

	bool IsUncontained(double near, double far) const
	{
		double horizontalDist = TMath::Sqrt(TMath::Power(m_position.X(),2)+TMath::Power(m_position.Y(),2));
		double verticalDist = TMath::Abs(m_position.Z());
		if ((horizontalDist < far && horizontalDist > near && verticalDist < 263) || (horizontalDist < far && verticalDist < 263+(far-60) && verticalDist > 263+(near-60)))
			return true;
		else
			return false;
	}

#ifndef NEW_CASC_STRUCTURE
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
#else
	bool IsLEDMatrixRun() const
	{
		bool isLEDMatrixRun = false;
		for (int i = 0; i < ledMatrixRuns[m_seasonID-2016][m_clusterID-1].size(); ++i)
		{
			if (m_runID == ledMatrixRuns[m_seasonID-2016][m_clusterID-1][i])
			{
				isLEDMatrixRun = true;
				break;
			}
		}
		return isLEDMatrixRun;
	}
#endif //NEW_CASCADE_STRUCTURE
};

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
	template<typename... args>
	DrawSingle(args... a) {drawsingle = new DrawType(a...);}

	void Fill(const Event& e) override {fillfunc(e);}
	void Draw() override {drawfunc();}
	void Save() override {canvas->Write();}
};

class CoincidenceFinder
{
public:
	struct Coincidence
	{
		int m_id;
		vector<int> m_indexes; //indexes of events in sortedEvents
		CoincidenceFinder* cfinder;

		Coincidence(CoincidenceFinder* cfinder)
		{
			this->m_id = -1;
			this->cfinder = cfinder;
		}

		Coincidence(const Coincidence&) = delete;
		Coincidence& operator=(Coincidence other) = delete;

		//~Coincidence() {}

		//returns angular distance of events number i and j in degrees
		double angDist(int i, int j) const
		{
		    return cfinder->sortedEvents[m_indexes[i]].angDist(cfinder->sortedEvents[m_indexes[j]]);
		}

		//adds event from sortedEvents
		void AddEvent(int index)
		{
			cfinder->sortedEvents[index].m_coincidenceID = this->m_id;
			this->m_indexes.push_back(index);
		}
	};

public:
	vector<Event> sortedEvents;
	vector<Coincidence*> coincidences;
	TH2F* random_coincidences = new TH2F("ran_coin","Random coincidences",5,0,5,2,2,4);

private:
	bool IsTAEC(int i, int j, long int maxT, double maxAD, double minE);
	bool IsTRPC(int i, int j, long int maxT, double maxD);

public:
	CoincidenceFinder(vector<Event> sortedEvents)
	{
		this->sortedEvents = sortedEvents;
	}

	CoincidenceFinder(const CoincidenceFinder&) = delete;
	CoincidenceFinder& operator=(CoincidenceFinder other) = delete;

	void WriteCStats();
	void WriteTAEC(long int maxT, double maxAD, double minE);

	template<typename... args>
	void FindCoincidences(bool(CoincidenceFinder::*IsCoin)(int, int, long int, args...), uint random_offset, long int maxdt, args... a);
	void FindTAEC(long int maxTimeDiff, double maxAngDist = 360, double minEnergy = 0) {FindCoincidences(&CoincidenceFinder::IsTAEC,0,maxTimeDiff,maxAngDist,minEnergy);}
	void WarnLEDMatrixRun(int minCSize = 3, long int maxT = 3600, double maxD = 5);
	void RandomCoincidences(long int maxT,double maxAD = 20);
};

class EventLoop
{
	typedef std::function<bool(const Event&)> FilterFn;

private:
	int year,cluster,startID,endID,startSeason,endSeason;

	vector<FilterFn> filters;
	vector<IDrawable*> drawables;
	Event current_ev;
#ifndef NEW_CASC_STRUCTURE
	TVector3* position = new TVector3();
	TVector3* mcPosition = new TVector3();
	TTimeStamp* eventTime = new TTimeStamp();
#else
	BRecoCascade* myCascade = NULL;
	BJointHeader* myHeader = NULL;
#endif //NEW_CASC_STRUCTURE

	map<TString,tuple<TGraph*,TGraph*>> flux_graphs;
	map<TString,TCanvas*> flux_canv;

public:
	vector<Event>   sortedEvents;
	vector<RunInfo> runs;
	vector<vector<const RunInfo*>> plotruns;
	TChain reconstructedCascades;
	TTree* filteredCascades;
	CoincidenceFinder* cfinder;

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

	void LoadReco(const char* env_p, bool v379);
	void LoadRunLogs(string dir);
	void SetUpTTrees();
	void UseLEDfilter();
	void UseContainedFilter(double dist);
	void UseLikelihoodFilter(double max);
	void UseEnergyFilter(double min);
	void UseUpGoingFilter();
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

//returns number of leap years since 2016 (copied from transformations.h UTCtoUnix function)
int GetLeapYears(int season);
//returns unix (1995) time of 01/04/YYYY 00:00:00
int GetStartTime(int season);
//returns unix (1995) time of 31/03/YYYY+1 23:59:59
int GetEndTime(int season);

#include "EventLoop.cpp"