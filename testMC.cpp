#include "MCevents.h"
#include "threading.h"
#include "transformations.h"

#include <fstream>

#include "TF1.h"
#include "TFile.h"
#include "TRandom2.h"
#include "TTree.h"
#include "TVirtualFitter.h"

#define PROFILING 0
#if PROFILING
	#include "Instrumentor.h"
	#define PROFILE_SCOPE(name) InstrumentationTimer timer##__LINE__(name)
	#define PROFILE_FUNCTION() PROFILE_SCOPE(__FUNCSIG__) 
#else
	#define PROFILE_SCOPE(name)
#endif

static std::vector<MCEvent> MCdata;
static TVirtualFitter* gFitter;

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

//sigma cut
static double numOfSigma   = 7;

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
                                   double timeSigma, double ra, double dec, double posSigma, double numOfSigma,
                                   double rate = 0)
{
	PROFILE_FUNCTION();

	TRandom2 rnd(0);

	int bEvents = numEvents;

	if(rate != 0) bEvents = rnd.Poisson(timeWindow*rate);
	double uTime = MJDtoUnix(MJDstart);

	std::vector<MCEvent>* output = new std::vector<MCEvent>;
	output->reserve(bEvents);

	//energy limits
	double minEnergy = 1;
	double maxEnergy = 10000000;

	//energy distrubution for atmospheric neutrinos
	TF1* energyDist = new TF1("energyDist", "[0]*std::pow(x,[1])", minEnergy, maxEnergy);

	//altitude distribution for uniform position distribution
	TF1* altDist = new TF1("altDist", "cos(x*3.14159265358979323/180)", -90, 90);

	double power = -3.7; //-3.7 for atmospheric neutrinos
	energyDist->SetParameter(1,power);

	//normalization (unnecesary, only for plotting)
	double area  = powerLawArea(power, minEnergy, maxEnergy);
	energyDist->SetParameter(0,1/area);

	//generating random values
	for (int i = 0; i < bEvents; ++i)
	{
		double time = UnixToMJD(uTime + 86400*timeWindow*rnd.Rndm());

		if(true)//abs(time - timeMean) < numOfSigma*timeSigma)
		{
			double altitude = altDist->GetRandom();
			double azimuth  = 360*rnd.Rndm();
	
			eqCoor pos 		= horToEq(horCoor(altitude,azimuth,uTime+time));
	
			//saving only events close to signal
			if(true)//angularDistance(ra,dec,pos.rAsc,pos.dec)<numOfSigma*posSigma)
			{
				double energy = energyDist->GetRandom();
				output->emplace_back(MCEvent{energy, pos.rAsc, pos.dec, time, 'B'});
			}
		}
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

	TRandom2 rnd(0);

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

	std::vector<MCEvent>* mixed  = getBackround(MJDstart, timeWindow, numBackround, signalMean, signalSigma,
												ra, dec, posSigma, 7);
	std::vector<MCEvent>* signal = getSignalGaus(numSignal, ra, dec, posSigma, signalMean, signalSigma);

	*mixed += *signal;
	//std::sort(mixed->begin(), mixed->end(), isSooner);

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
	PROFILE_FUNCTION();

	double partSignal = nSignal/MCdata.size();
	return partSignal*signalProbability(ev) + (1-partSignal)*backgroundProbability(ev);
}

double testStatistic(double bestFit)
{
	PROFILE_FUNCTION();

	double lBest = 0;

	for(const MCEvent& ev : MCdata)
	{
		lBest += std::log(probability(ev,bestFit));
	}

	double lNone = 0;

	for(const MCEvent& ev : MCdata)
	{
		lNone += std::log(backgroundProbability(ev));
	}

	return -2*(lNone-lBest);
}

//logLikelihood with output parameter outL and input parameters array par
void logLikelihood(int& npar, double* gin, double& outL, double* par, int iflag)
{
	PROFILE_FUNCTION();

	double& nSignal = par[0]; //creating an alias for par[0]

	outL = 0;

	for(const MCEvent& ev : MCdata)
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

void fit(double& nSignal, double& nSignalSigma, bool print = true)
{
	PROFILE_FUNCTION();

	SetFitter(1,print);
	gFitter->SetParameter(0,"nSignal",nSignal,0.01,0,MCdata.size());

	double arglist[10] = {500,0.01}; //max iterations, step size
	gFitter->ExecuteCommand("MIGRAD",arglist,2); //last one num of prints
	nSignal      = gFitter->GetParameter(0);
	nSignalSigma = gFitter->GetParError(0);
}

// ./testMC outpath backCount signalCount timeWindow timeSigma posSigma (numOfSigma = 7)
int main(int argc, char** argv)
{
#if PROFILING
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

			if(argc == 7) numOfSigma = std::stoi(argv[7]);
		}

		MCdata = *getMixed(backCount,signalCount,sigTimeMean-timeWindow/2,timeWindow,sigTimeMean,sigTimeSigma,sigPosSigma,0,0);

		double nSignal = 0;
		double nSignalSigma;

		fit(nSignal,nSignalSigma,false);

		std::ofstream outf{outpath, std::ios::app};
		outf << testStatistic(nSignal) << '\n';
	}
#if PROFILING
	Instrumentor::Get().EndSession();
#endif

	return 0;
}