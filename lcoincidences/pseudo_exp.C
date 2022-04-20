#define PROFILLING 0
#include "../EventLoop.h"
// #include "../threading.h"
#include "../transformations.h"

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
const static double sig_gamma = -2;
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
static double eNormSig     = powerLawArea(sig_gamma,minEnergy,maxEnergy);
static double eNormBack    = powerLawArea(bckg_gamma,minEnergy,maxEnergy);

TF1* bckg_energyDist = new TF1("energyDist", "std::pow(x,[0])",minEnergy,maxEnergy);
TF1* sig_energyDist  = new TF1("energyDist", "std::pow(x,[0])",minEnergy,maxEnergy);

std::function<double(double*,double*)> thetaDistEval = [](double* x, double* par)
{
	double xx = x[0];

	return max(0.0,abs(sin(xx))*costheta->Eval(cos(xx)));
};

TF1* bckg_thetaDist  = new TF1("thetaDist",thetaDistEval,0,TMath::Pi(),0);
// TH1F* enhist = new TH1F("enhist","enhist",100,1,1000);

static vector<double> sigprobs;
static vector<double> bkgprobs;
static int nSimulEvents;
static double lastFitResult;

void GetProbs()
{
	PROFILE_FUNCTION();

	for(const BasicEvent& ev : simulatedEvents)
	{			
		//folded gaussian distribution
		double positionProb = TMath::Gaus(ev.angDist(sigRa,sigDec), 0, sigPosSigma)*2;
		double energyProb	= std::pow(ev.m_energy,sig_gamma)/eNormSig;
		sigprobs.push_back(positionProb*energyProb);

		positionProb = abs(sin(ev.m_theta)*costheta->Eval(cos(ev.m_theta)));
		energyProb	= std::pow(ev.m_energy,bckg_gamma)/eNormBack;
		bkgprobs.push_back(positionProb*energyProb);
	}
}

void GetLastProb()
{
	PROFILE_FUNCTION();

	const BasicEvent& ev = simulatedEvents.back();

	//folded gaussian distribution
	double positionProb = TMath::Gaus(ev.angDist(sigRa,sigDec), 0, sigPosSigma)*2;
	double energyProb	= std::pow(ev.m_energy,sig_gamma)/eNormSig;
	sigprobs.push_back(positionProb*energyProb);

	positionProb = abs(sin(ev.m_theta)*costheta->Eval(cos(ev.m_theta)));
	energyProb	= std::pow(ev.m_energy,bckg_gamma)/eNormBack;
	bkgprobs.push_back(positionProb*energyProb);
}

double testStatistic(double bestFit)
{
	PROFILE_FUNCTION();

	double lNone = 0;

	for(int i = 0; i < nSimulEvents; i++)
	{
		lNone += std::log(simulatedEvents.size()*bkgprobs[i]);
	}

	return (lastFitResult-lNone);
}

//logLikelihood with output parameter outL and input parameters array par
void logLikelihood(int& npar, double* gin, double& outL, double* par, int iflag)
{
	PROFILE_FUNCTION();
	double& nSignal = par[0]; //creating an alias for par[0]

	outL = 0;

	for(int i = 0; i < nSimulEvents; i++)
	{
		outL -= std::log(nSignal*sigprobs[i] + simulatedEvents.size()*bkgprobs[i]);
	}

	outL += nSignal;
	lastFitResult = -outL;
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

void generate_signal()
{
	PROFILE_FUNCTION();

	BasicEvent sig_ev;
	
	sig_ev.m_energy = sig_energyDist->GetRandom();

	//converting right ascension and declination to radians
	double raRad  = degToRad(sigRa);
	double decRad = degToRad(sigDec);

	//angular distance and angle in radians
	double dev 	  = abs(gRandom->Gaus(0,degToRad(sigPosSigma))); //only works for small enough sigma (dev<PI)
	double angle  = 2*PI*gRandom->Rndm();

	//transformation using horToEq with zero local sidereal time
	eqCoor coorOut = shiftSpherTrans(radToDeg(dev),radToDeg(angle),raRad,decRad);
	sig_ev.m_rightAscension = coorOut.rAsc;
	sig_ev.m_declination    = coorOut.dec;

	sig_ev.m_eventTime = TTimeStamp((tEnd-tStart)*gRandom->Rndm()+tStart,0);

	sig_ev.computeThetaPhi();

	simulatedEvents.push_back(sig_ev);
}

void FitAndOutput(int signal_events, double& nSignal, double& nSignalSigma, ofstream& outf, ofstream& outf2)
{
	GetProbs();

	for(int sig = 0; sig <= signal_events; sig++)
	{
		if(sig != 0) {generate_signal(); GetLastProb();}

		nSignal = 1;

		fit(nSignal,nSignalSigma);

		PROFILE_SCOPE("Writing to files.");
		outf  << nSignal << "\n";
		outf2 << testStatistic(nSignal) << "\n";					
	}

	simulatedEvents.resize(simulatedEvents.size() - signal_events);

	sigprobs.clear();
	bkgprobs.clear();
}

void RunSimulation(int signal_events, double input_dec, double end_dec, double step_dec, double ra_step, int id, int nSimulations)
{
	double nSignal;
	double nSignalSigma;

	SetFitter(1,false);
	gRandom = new TRandom3(0);
	bckg_energyDist->SetParameter(0,bckg_gamma);
	sig_energyDist->SetParameter(0,sig_gamma);
	bckg_thetaDist->SetNormalized(true);
	// bckg_thetaDist->Draw(); //for drawing costheta must not be deleted!!!

	sigprobs.reserve(nSimulEvents+signal_events);
	bkgprobs.reserve(nSimulEvents+signal_events);

	if (step_dec == 0)
	{
		cout << "IF" << endl;
		sigDec = input_dec;

		if (ra_step == 0) //dont iterate ra or dec
		{
			string outpath  = "./data/data_nSign_dec_";
			string outpath2 = "./data/data_tStat_dec_";

			outpath  += to_string(sigDec) + "_" + id + ".txt";
			outpath2 += to_string(sigDec) + "_" + id + ".txt";

			std::ofstream outf{outpath, std::ios::app};
			std::ofstream outf2{outpath2, std::ios::app};

			for (int i = 0; i < nSimulations; ++i)
			{
				generate_background(nSimulEvents);
				FitAndOutput(signal_events, nSignal, nSignalSigma, outf, outf2);
			}
		}
		else //iterate ra but not dec
		{
			for (int i = 0; i < nSimulations; ++i)
			{
				generate_background(nSimulEvents);

				for (sigRa = -180; sigRa < 180; sigRa += ra_step)
				{
					string outpath  = "./data/data_nSign_dec_";
					string outpath2 = "./data/data_tStat_dec_";

					outpath  += to_string(sigDec) + "_" + to_string(sigRa) + "_" + id + ".txt";
					outpath2 += to_string(sigDec) + "_" + to_string(sigRa) + "_" + id + ".txt";

					std::ofstream outf{outpath, std::ios::app};
					std::ofstream outf2{outpath2, std::ios::app};

					FitAndOutput(signal_events, nSignal, nSignalSigma, outf, outf2);

					outf.close();
					outf2.close();
				}
			}
		}
	}

	else //iterate dec but not ra
	{
		for (int i = 0; i < nSimulations; ++i)
		{
			generate_background(nSimulEvents);

			for (sigDec = input_dec; sigDec <= end_dec; sigDec += step_dec)
			{
				string outpath  = "./data/data_nSign_dec_";
				string outpath2 = "./data/data_tStat_dec_";

				outpath  += to_string(sigDec) + "_" + id + ".txt";
				outpath2 += to_string(sigDec) + "_" + id + ".txt";

				std::ofstream outf{outpath, std::ios::app};
				std::ofstream outf2{outpath2, std::ios::app};

				FitAndOutput(signal_events, nSignal, nSignalSigma, outf, outf2);	

				outf.close();
				outf2.close();			
			}
		}
	}
}

int pseudo_exp(int signal_events, double input_dec, int id, double end_dec = 0, double step_dec = 0, int iterate_ra = 0, int nSimulations = 10000)
{
	PROFILLING_START_UNIQUE("pseudo_exp");
	PTIMER_START("MAIN",MAIN);

	gErrorIgnoreLevel = 6001; //no ROOT errors please
	nSimulEvents = 3365;

	TFile* pdf_input = TFile::Open("cos_theta.root","READ");

	pdf_input->GetObject("f_spline",costheta);
	// costheta->SetNormalized(true);
	// cout << "costheta(-0.5): " << (*costheta)(-0.5) << endl;

	if (iterate_ra == 1) RunSimulation(signal_events, input_dec, 0, 0, 10, id, nSimulations);
	else RunSimulation(signal_events, input_dec, end_dec, step_dec, 0, id, nSimulations);

	delete pdf_input;
	delete costheta;
	delete gFitter;
	// enhist->Draw();

	PTIMER_STOP(MAIN);
	PROFILLING_END();

	return 0;
}

int main(int argc, char** argv)
{
	double input_dec;
	int id, signal_events;

	if(argc < 3) cerr << "Not enough arguments!\n";
	else
	{
		signal_events = stoi(argv[1]);
		input_dec     = stod(argv[2]);
		id 		      = stoi(argv[3]);
	}

	double end_dec = 90;
	if(argc > 3) end_dec = stod(argv[4]);

	double step_dec = 0;
	if(argc > 4) step_dec = stod(argv[5]);

	int iterate_ra = 0;
	if(argc > 5) iterate_ra = stoi(argv[6]);

	return pseudo_exp(signal_events,input_dec,id,end_dec,step_dec,iterate_ra);
}