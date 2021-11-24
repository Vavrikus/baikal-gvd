#include "../EventLoop.h"
#include "../threading.h"
#include "../transformations.h"

#define PROFILLING 0
#include "../Instrumentor.h"

#include <fstream>
#include <mutex>
#include <thread>
#include <algorithm>

#include "TF1.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TVirtualFitter.h"

static std::vector<BasicEvent> simulatedEvents;
static TVirtualFitter* gFitter;

//data cosine of theta distribution
static TF1* costheta;

//signal distribution degrees
static double sigRa 	   = 0;
static double sigDec	   = 0;
static double sigPosSigma  = 5;

//energy window
static double minEnergy    = 1;
static double maxEnergy	   = 10000000;

//energy coefficients
const static double s_gamma = -2;
const static double bckg_gamma = -3.7;

time_t tStart = GetStartTime(2019)+unix1995;
time_t tEnd   = GetEndTime(2020)+unix1995;

//return definite integral of power law
double powerLawArea(double power, double xmin, double xmax)
{
	PROFILE_FUNCTION();
	return std::pow(xmax, power+1)/(power+1) - std::pow(xmin, power+1)/(power+1);
}

//energy distribution normalization parameter
static double eNormSig     = powerLawArea(s_gamma,minEnergy,maxEnergy);
static double eNormBack    = powerLawArea(bckg_gamma,minEnergy,maxEnergy);

TF1* bckg_energyDist = new TF1("energyDist", "std::pow(x,[0])",minEnergy,maxEnergy);

std::function<double(double*,double*)> thetaDistEval = [](double* x, double* par)
{
	double xx = x[0];

	return max(0.0,abs(sin(xx))*costheta->Eval(cos(xx)));
};

TF1* bckg_thetaDist  = new TF1("thetaDist",thetaDistEval,0,TMath::Pi(),0);
// TH1F* enhist = new TH1F("enhist","enhist",100,1,1000);

double signalProbability(const BasicEvent& ev)
{
	PROFILE_FUNCTION();
	//folded gaussian distribution
	double positionProb = TMath::Gaus(ev.angDist(sigRa,sigDec), 0, sigPosSigma)*2;
	double energyProb	= std::pow(ev.m_energy,s_gamma)/eNormSig;

	return positionProb*energyProb;
}

double backgroundProbability(const BasicEvent& ev)
{
	PROFILE_FUNCTION();
	double positionProb = abs(sin(ev.m_theta)*costheta->Eval(cos(ev.m_theta)));
	double energyProb	= std::pow(ev.m_energy,bckg_gamma)/eNormBack;

	return positionProb*energyProb;
}

double probability(const BasicEvent& ev, double nSignal)
{
	PROFILE_FUNCTION();
	return nSignal*signalProbability(ev) + simulatedEvents.size()*backgroundProbability(ev);
}

double testStatistic(double bestFit)
{
	PROFILE_FUNCTION();

	double lBest = 0;

	for(const BasicEvent& ev : simulatedEvents)
	{
		lBest += std::log(probability(ev,bestFit));
	}

	lBest -= bestFit;

	double lNone = 0;

	for(const BasicEvent& ev : simulatedEvents)
	{
		lNone += std::log(probability(ev,0));
	}

	return (lBest-lNone);
}

//logLikelihood with output parameter outL and input parameters array par
void logLikelihood(int& npar, double* gin, double& outL, double* par, int iflag)
{
	PROFILE_FUNCTION();
	double& nSignal = par[0]; //creating an alias for par[0]

	outL = 0;

	for(const BasicEvent& ev : simulatedEvents)
	{
		outL -= std::log(probability(ev,nSignal));
	}

	outL += nSignal;
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

	gFitter->SetParameter(0,"nSignal",nSignal,0.1,0,simulatedEvents.size());

	double arglist[10] = {500,0.1}; //max iterations, step size
	gFitter->ExecuteCommand("MIGRAD",arglist,2); //last one num of prints
	nSignal      = gFitter->GetParameter(0);
	nSignalSigma = gFitter->GetParError(0);
}

void generate_background(int events)
{
	PROFILE_FUNCTION();
	simulatedEvents.clear();

	for (int i = 0; i < events; ++i)
	{
		BasicEvent bckg_ev;
		
		bckg_ev.m_energy = bckg_energyDist->GetRandom();

		bckg_ev.m_theta  = bckg_thetaDist->GetRandom();
		bckg_ev.m_phi    = 2*TMath::Pi()*gRandom->Rndm();

		bckg_ev.m_eventTime = TTimeStamp((tEnd-tStart)*gRandom->Rndm()+tStart,0);

		bckg_ev.computeRaDec();

		simulatedEvents.push_back(bckg_ev);

		// cout << "m_energy: " << bckg_ev.m_energy << "\n";
		// cout << "m_theta: " << bckg_ev.m_theta << "\n";
		// cout << "m_phi: " << bckg_ev.m_phi << "\n";
		// cout << "m_eventTime: " << bckg_ev.m_eventTime << "\n";
		// cout << "m_declination: " << bckg_ev.m_declination << "\n";
		// cout << "m_rightAscension: " << bckg_ev.m_rightAscension << "\n\n";

		// enhist->Fill(bckg_ev.m_energy);
	}
}

int pseudo_exp(double input_dec, int id, int nSimulations = 10000)
{
	gErrorIgnoreLevel = 6001; //no ROOT errors please

	sigDec = input_dec;

	string outpath  = "./data/data_nSign_dec_";
	string outpath2 = "./data/data_tStat_dec_";

	outpath  += to_string(input_dec) + "_" + id + ".txt";
	outpath2 += to_string(input_dec) + "_" + id + ".txt";

	std::ofstream outf{outpath, std::ios::app};
	std::ofstream outf2{outpath2, std::ios::app};

	TFile* pdf_input = TFile::Open("cos_theta.root","READ");

	pdf_input->GetObject("f_spline",costheta);
	// costheta->SetNormalized(true);
	// cout << "costheta(-0.5): " << (*costheta)(-0.5) << endl;

	double nSignal;
	double nSignalSigma;

	SetFitter(1,false);
	gRandom = new TRandom3(0);
	bckg_energyDist->SetParameter(0,bckg_gamma);
	bckg_thetaDist->SetNormalized(true);
	// bckg_thetaDist->Draw(); //for drawing costheta must not be deleted!!!

	for (int i = 0; i < nSimulations; ++i)
	{
		generate_background(3365);

		nSignal = 1;

		fit(nSignal,nSignalSigma);

		outf  << nSignal << "\n";
		outf2 << testStatistic(nSignal) << "\n";
	}

	delete pdf_input;
	delete costheta;
	delete gFitter;
	// enhist->Draw();

	return 0;
}

int main(int argc, char** argv)
{
	double input_dec;
	int id;

	if(argc < 3) cerr << "Not enough arguments!\n";
	else
	{
		input_dec = stod(argv[1]);
		id 		  = stoi(argv[2]);
	}

	return pseudo_exp(input_dec,id);
}