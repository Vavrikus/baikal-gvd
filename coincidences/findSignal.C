#include "../transformations.h"
#include "MCevents.h"
#include "../threading.h"

#include <fstream>

#include "TVirtualFitter.h"
#include "TFile.h"
#include "TTree.h"

#define ASYNC 0

static std::vector<MCEvent> MCdata;
static TVirtualFitter* gFitter;

//signal distribution
static double sigTimeMean  = 58000;
static double sigTimeSigma = 0.1;
static double sigRa 	   = 0;
static double sigDec	   = 0;
static double sigPosSigma  = 5;

//chosen time window
static double timeWindow   = 365.25;

//energy window
static double minEnergy    = 1;
static double maxEnergy	   = 10000000;

//energy distribution normalization parameter
static double eNormSig     = powerLawArea(-2,minEnergy,maxEnergy);
static double eNormBack    = powerLawArea(-3.7,minEnergy,maxEnergy);

//parse MC data from root file to vector of MCEvents
std::vector<MCEvent>* readMC(const char* data_path)
{
    //access tree with data inside root file
    TFile file(data_path);
	TTree* events;

	file.GetObject("MC_events",events);
    
    std::vector<MCEvent>* dataOut = new std::vector<MCEvent>;
    dataOut->reserve(events->GetEntries());

    //variables for data
    char type;
    double mjd,ra,dec,energy; 
    
    events->SetBranchAddress("type",   &type);
    events->SetBranchAddress("mjd",    &mjd);
    events->SetBranchAddress("ra",     &ra);
    events->SetBranchAddress("dec",    &dec);
    events->SetBranchAddress("energy", &energy);

    //parsing data to MCEvent objects
    for (int i = 0; i < events->GetEntries(); ++i)
    {
        events->GetEntry(i);
        dataOut->emplace_back(MCEvent{energy, ra, dec, mjd, type});
    }

    return dataOut;
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
	double timeProb     = 1/(14*sigTimeSigma);
	double positionProb = 1/(capAngle(7*sigPosSigma));
	double energyProb	= std::pow(ev.energy,-3.7)/eNormBack;

	return timeProb*positionProb*energyProb;
}

double probability(const MCEvent& ev, double nSignal)
{
	double partSignal = nSignal/MCdata.size();
	return partSignal*signalProbability(ev) + (1-partSignal)*backgroundProbability(ev);
}

#if ASYNC
double func(int i, double nSignal)
{
	return -std::log(probability(MCdata[i],nSignal));
}
#endif

double testStatistic(double bestFit)
{
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
	//Timer t("logLike");
	double& nSignal = par[0]; //creating an alias for par[0]

#if ASYNC
	outL = threadProbabilityFor(8, MCdata.size(), func, nSignal);
#else
	outL = 0;

	for(const MCEvent& ev : MCdata)
	{
		outL -= std::log(probability(ev,nSignal));
	}
#endif
	//std::cout << std::setprecision(16) << "outL: " << outL << '\n';
}

void SetFitter(int parameters, bool print = true)
{
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
	SetFitter(1,print);
	gFitter->SetParameter(0,"nSignal",nSignal,0.01,0,MCdata.size());

	double arglist[10] = {500,0.01}; //max iterations, step size
	gFitter->ExecuteCommand("MIGRAD",arglist,2); //last one num of prints
	nSignal      = gFitter->GetParameter(0);
	nSignalSigma = gFitter->GetParError(0);
}

void findSignal(const char* inpath = "experiments_data/MCevents.root", const char* outpath = "")
{
	MCdata = *readMC(inpath);

	double nSignal = MCdata.size()/2;
	double nSignalSigma;

	fit(nSignal,nSignalSigma,false);

	std::ofstream outf{outpath, std::ios::app};
	outf << testStatistic(nSignal) << '\n';

	//std::cout << "Signal events: " << nSignal << " +- " << nSignalSigma << "\n";
	//std::cout << "Test statistic: " << testStatistic(nSignal) << "\n";

	//for(auto ev : MCdata) std::cout << "Type: " << ev.type << " Signal probability: " << signalProbability(ev) << " Background probability: " << backgroundProbability(ev) << '\n';
}

// ./findSignal inpath outpath timeWindow timeSigma posSigma
int main(int argc, char** argv)
{
	if(argc > 1)
	{
		timeWindow = std::stod(argv[3]);
		sigTimeSigma = std::stod(argv[4]);
		sigPosSigma = std::stod(argv[5]);

		findSignal(argv[1], argv[2]);
	}

	else findSignal();

	return 0;
}