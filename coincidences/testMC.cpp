#include "MCevents.h"
#include "../threading.h"
#include "../transformations.h"

#define PROFILLING 0
#include "../Instrumentor.h"

#include <fstream>

#include "TF1.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TVirtualFitter.h"

static std::vector<MCEvent> fitData;
static TVirtualFitter* gFitter;

//number of simulated fits
static double numOfSimulations = 10000;

//generated events
static int backCount;
static int signalCount;

//signal distribution
static double sigTimeMean  = 58000;
static double sigTimeSigma;
static double sigRa 	   = 0;
static double sigDec	   = 0;
static double sigPosSigma;

//chosen time window
static double timeWindow;

//energy window
static double minEnergy    = 1;
static double maxEnergy	   = 10000000;

//energy distribution normalization parameter
static double eNormSig     = powerLawArea(-2,minEnergy,maxEnergy);
static double eNormBack    = powerLawArea(-3.7,minEnergy,maxEnergy);

//operator overloading for appending std::vector
template<typename T>
void operator+=(std::vector<T>& v1, const std::vector<T>& v2)
{
    v1.insert(v1.end(), v2.begin(), v2.end());
}

//timeWindow in days and rate in events per day
//MJDstart is modified julian date
//either generate number of events from rate or take them as parameter
std::vector<MCEvent>* getBackround(double MJDstart, double timeWindow, int numEvents, double timeMean,
                                   double timeSigma, double ra, double dec, double posSigma, double rate = 0)
{
	PROFILE_FUNCTION();

	TRandom3 rnd(0);

	int bEvents = numEvents;

	if(rate != 0) bEvents = rnd.Poisson(timeWindow*rate);
	double uTime = MJDtoUnix(MJDstart);

	std::vector<MCEvent>* output = new std::vector<MCEvent>;
	output->reserve(bEvents);

	//energy limits
	double minEnergy = 1;
	double maxEnergy = 10000000;

	#if PROFILING
		InstrumentationTimer t("Setting up energy distribution");
	#endif
	//energy distrubution for atmospheric neutrinos
	TF1* energyDist = new TF1("energyDist", "[0]*std::pow(x,[1])", minEnergy, maxEnergy);

	double power = -3.7; //-3.7 for atmospheric neutrinos
	energyDist->SetParameter(1,power);

	#if PROFILING
		t.Stop();
	#endif

	#if PROFILING
		InstrumentationTimer t2("Setting up altitude distribution");
	#endif
	//altitude distribution for uniform position distribution
	TF1* altDist = new TF1("altDist", "cos(x*3.14159265358979323/180)", -90, 90);

	//normalization (unnecesary, only for plotting)
	double area  = powerLawArea(power, minEnergy, maxEnergy);
	energyDist->SetParameter(0,1/area);
	#if PROFILING
		t2.Stop();
	#endif

	//generating random values
	for (int i = 0; i < bEvents; ++i)
	{
		PROFILE_SCOPE("For loop");

		double time = UnixToMJD(uTime + 86400*timeWindow*rnd.Rndm());
		double energy = energyDist->GetRandom();
		double altitude = altDist->GetRandom();
		double azimuth  = 360*rnd.Rndm();

		eqCoor pos 		= horToEq(horCoor(altitude,azimuth,uTime+time));

		output->emplace_back(MCEvent{energy, pos.rAsc, pos.dec, time, 'B'});
	}

	delete energyDist;
	delete altDist;

	return output;
}

//signal events with gaussian time distribution
//position in degrees, time in days (modified julian date)
std::vector<MCEvent>* getSignalGaus(int eventCount, double ra, double dec, double posSigma,
									double timeMean, double timeSigma)
{
	PROFILE_FUNCTION();

	TRandom3 rnd(0);

	std::vector<MCEvent>* output = new std::vector<MCEvent>;
	output->reserve(eventCount);

	//energy limits
	double minEnergy = 1;
	double maxEnergy = 10000000;

	//energy distrubution for astrophysical neutrinos
	TF1* energyDist = new TF1("energyDist", "[0]*std::pow(x,[1])", minEnergy, maxEnergy);

	double power = -2; //-2 for astrophysical neutrinos
	energyDist->SetParameter(1,power);

	//normalization (unnecesary, only for plotting)
	double area  = powerLawArea(power, minEnergy, maxEnergy);
	energyDist->SetParameter(0,1/area);

	//converting right ascension and declination to radians
	double raRad  = degToRad(ra);
	double decRad = degToRad(dec);

	//generating random values
	for(int i = 0; i < eventCount; i++)
	{
		double energy = energyDist->GetRandom();
		double time   = rnd.Gaus(timeMean, timeSigma);

		//angular distance and angle in radians
		double dev 	  = abs(rnd.Gaus(0,degToRad(posSigma))); //only works for small enough sigma (dev<PI)
		double angle  = 2*PI*rnd.Rndm();

		//transformation using horToEq with zero local sidereal time
		eqCoor coorOut = shiftSpherTrans(radToDeg(dev),radToDeg(angle),raRad,decRad);

		output->emplace_back(MCEvent{energy,coorOut.rAsc,coorOut.dec,time,'S'});
	}

	delete energyDist;

	return output;
}

//get vector of MC events sorted by time
std::vector<MCEvent>* getMixed(int numBackround, int numSignal, double MJDstart, double timeWindow, 
							   double signalMean, double signalSigma, double posSigma, double ra, double dec)
{
	PROFILE_FUNCTION();

	std::vector<MCEvent>* background  = getBackround(MJDstart, timeWindow, numBackround, signalMean, signalSigma,
												ra, dec, posSigma);
	std::vector<MCEvent>* signal 	  = getSignalGaus(numSignal, ra, dec, posSigma, signalMean, signalSigma);

	std::vector<MCEvent>* mixed = new std::vector<MCEvent>();
	mixed->reserve(numBackround+numSignal);

	for (int i = 0; i < numBackround+numSignal; ++i)
	{
		int simulID = i/(backCount+signalCount);
		int eventID = i%(backCount+signalCount);

		if(eventID < backCount)	mixed->push_back((*background)[simulID*backCount + eventID]);
		else mixed->push_back((*signal)[simulID*signalCount - backCount + eventID]);
	}
	
	delete background;
	delete signal;

	return mixed;
}

double signalProbability(const MCEvent& ev)
{
	double timeProb     = TMath::Gaus(ev.mjd,sigTimeMean,sigTimeSigma);

	//folded gaussian distribution
	double positionProb = TMath::Gaus(angularDistance(ev.ra, ev.dec, sigRa, sigDec), 0, sigPosSigma)*2;
	double energyProb	= std::pow(ev.energy,-2)/eNormSig;

	return timeProb*positionProb*energyProb;
}

double backgroundProbability(const MCEvent& ev)
{
	double timeProb     = 1/timeWindow;//(14*sigTimeSigma);
	double positionProb = 1/(4*PI);//(capAngle(7*sigPosSigma));
	double energyProb	= std::pow(ev.energy,-3.7)/eNormBack;

	return timeProb*positionProb*energyProb;
}

double probability(const MCEvent& ev, double nSignal)
{
	double partSignal = nSignal/fitData.size();
	return partSignal*signalProbability(ev) + (1-partSignal)*backgroundProbability(ev);
}

double testStatistic(double bestFit)
{
	PROFILE_FUNCTION();

	double lBest = 0;

	for(const MCEvent& ev : fitData)
	{
		lBest += std::log(probability(ev,bestFit));
	}

	double lNone = 0;

	for(const MCEvent& ev : fitData)
	{
		lNone += std::log(backgroundProbability(ev));
	}

	return -2*(lNone-lBest);
}

//logLikelihood with output parameter outL and input parameters array par
void logLikelihood(int& npar, double* gin, double& outL, double* par, int iflag)
{
	double& nSignal = par[0]; //creating an alias for par[0]

	outL = 0;

	for(const MCEvent& ev : fitData)
	{
		outL -= std::log(probability(ev,nSignal));
	}
}

void SetFitter(int parameters, bool print = true)
{
	PROFILE_FUNCTION();
	
	gFitter = TVirtualFitter::Fitter(nullptr,parameters); // the second number is number of parameters
	gFitter->SetFCN(logLikelihood);

	if(!print)
	{
		double arg = -1;
		gFitter->ExecuteCommand("SET PRINTOUT",&arg,1);
		gFitter->ExecuteCommand("SET NOW", &arg ,1);
	}
}

void fit(double& nSignal, double& nSignalSigma)
{
	PROFILE_FUNCTION();

	gFitter->SetParameter(0,"nSignal",nSignal,0.1,0,backCount+signalCount);

	double arglist[10] = {500,0.1}; //max iterations, step size
	gFitter->ExecuteCommand("MIGRAD",arglist,2); //last one num of prints
	nSignal      = gFitter->GetParameter(0);
	nSignalSigma = gFitter->GetParError(0);
}

// ./testMC outpath backCount signalCount timeWindow timeSigma posSigma (numOfSimulations = 10000)
int main(int argc, char** argv)
{
#if PROFILLING
	Instrumentor::Get().BeginSession("Session Name");
#endif
	{
		PROFILE_SCOPE("MAIN");

		std::string outpath;

		if(argc < 6) std::cerr << "Not enough arguments!\n";
		else
		{
			outpath        = argv[1];
			backCount      = std::stoi(argv[2]);
			signalCount    = std::stoi(argv[3]);
			timeWindow 	   = std::stod(argv[4]);
			sigTimeSigma   = std::stod(argv[5]);
			sigPosSigma    = std::stod(argv[6]);

			//if(argc == 7) numOfSimulations = std::stoi(argv[7]);
		}
		
		std::ofstream outf{outpath, std::ios::app};

		double nSignal;
		double nSignalSigma;

		SetFitter(1,false);

		std::vector<MCEvent> allData = *getMixed(backCount*numOfSimulations,signalCount*numOfSimulations,sigTimeMean-timeWindow/2,
						   						 timeWindow,sigTimeMean,sigTimeSigma,sigPosSigma,0,0);

		for(int i = 0; i < numOfSimulations; i++)
		{
			nSignal = 0;

			std::vector<MCEvent>::const_iterator first = allData.begin() + i*(backCount+signalCount);
			std::vector<MCEvent>::const_iterator last  = allData.begin() + (i+1)*(backCount+signalCount);

			fitData = std::vector<MCEvent>(first, last);

			fit(nSignal,nSignalSigma);

			outf << testStatistic(nSignal) << '\n';
		}
	}
#if PROFILLING
	Instrumentor::Get().EndSession();
#endif

	return 0;
}