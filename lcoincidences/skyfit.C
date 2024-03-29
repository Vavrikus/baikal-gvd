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

//1e-5 probability bounds
const vector<double> fit_bounds = {10.9983,10.888,10.6235,10.4007,10.3206,10.3551,10.3797,10.5331,10.6906,10.8716,11.1583,11.6488,12.2753,13.1082,14.1282,15.2548,16.4592,17.9392,19.7688,22.1002,24.8501,28.3045,31.8804,35.5551,39.5829,44.0158,49.259,54.5009,58.8381,61.8463,63.9316,65.1841,66.6971,68.6378,70.9946,72.7542,73.3589};

static std::vector<Event> sortedEvents;
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

//return definite integral of power law
double powerLawArea(double power, double xmin, double xmax)
{
	return std::pow(xmax, power+1)/(power+1) - std::pow(xmin, power+1)/(power+1);
}

//energy distribution normalization parameter
static double eNormSig     = powerLawArea(-2,minEnergy,maxEnergy);
static double eNormBack    = powerLawArea(-3.7,minEnergy,maxEnergy);

double signalProbability(const Event& ev)
{
	//folded gaussian distribution
	double positionProb = TMath::Gaus(ev.angDist(sigRa,sigDec), 0, sigPosSigma)*2;
	double energyProb	= std::pow(ev.m_energy,-2)/eNormSig;

	return positionProb*energyProb;
}

double backgroundProbability(const Event& ev)
{
	double positionProb = max(0.0,abs(sin(ev.m_theta))*costheta->Eval(cos(ev.m_theta)));
	double energyProb	= std::pow(ev.m_energy,-3.7)/eNormBack;

	return positionProb*energyProb;
}

double probability(const Event& ev, double nSignal)
{
	return nSignal*signalProbability(ev) + sortedEvents.size()*backgroundProbability(ev);
}

double testStatistic(double bestFit)
{
	PROFILE_FUNCTION();

	double lBest = 0;

	for(const Event& ev : sortedEvents)
	{
		lBest += std::log(probability(ev,bestFit));
	}

	lBest -= bestFit;

	double lNone = 0;

	for(const Event& ev : sortedEvents)
	{
		lNone += std::log(probability(ev,0));
	}

	return (lBest-lNone);
}

//logLikelihood with output parameter outL and input parameters array par
void logLikelihood(int& npar, double* gin, double& outL, double* par, int iflag)
{
	double& nSignal = par[0]; //creating an alias for par[0]

	outL = 0;

	for(const Event& ev : sortedEvents)
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

	gFitter->SetParameter(0,"nSignal",nSignal,0.1,0,sortedEvents.size());

	double arglist[10] = {500,0.1}; //max iterations, step size
	gFitter->ExecuteCommand("MIGRAD",arglist,2); //last one num of prints
	nSignal      = gFitter->GetParameter(0);
	nSignalSigma = gFitter->GetParError(0);
}

int skyfit()
{
	TFile* pdf_input = TFile::Open("cos_theta.root","READ");
	TFile* prob_input = TFile::Open("prob.root","READ");

	TTree* filteredCascades;
	pdf_input->GetObject("f_spline",costheta);
	pdf_input->GetObject("filteredCascades",filteredCascades);

	vector<TF1*>* f_splines;
	prob_input->GetObject("f_splines",f_splines);

	vector<TF1*>* f_exps;
	prob_input->GetObject("f_exps",f_exps);

	costheta->SetNormalized(true);
	cout << "costheta(-0.5): " << (*costheta)(-0.5) << endl; 

	Event current_ev;
	TVector3* position = new TVector3();
	TVector3* mcPosition = new TVector3();
	TTimeStamp* eventTime = new TTimeStamp();

	filteredCascades->SetBranchAddress("seasonID", 			   &current_ev.m_seasonID);
	filteredCascades->SetBranchAddress("clusterID", 		   &current_ev.m_clusterID);
	filteredCascades->SetBranchAddress("runID", 			   &current_ev.m_runID);
	filteredCascades->SetBranchAddress("eventID", 			   &current_ev.m_eventID);
	filteredCascades->SetBranchAddress("nHits", 			   &current_ev.m_nHits);
	filteredCascades->SetBranchAddress("nHitsAfterCaus", 	   &current_ev.m_nHitsAfterCaus);
	filteredCascades->SetBranchAddress("nStringsAfterCaus",    &current_ev.m_nStringsAfterCaus);
	filteredCascades->SetBranchAddress("chi2AfterCaus", 	   &current_ev.m_chi2AfterCaus);
	filteredCascades->SetBranchAddress("nHitsAfterTFilter",    &current_ev.m_nHitsAfterTFilter);
	filteredCascades->SetBranchAddress("nStringsAfterTFilter", &current_ev.m_nStringsAfterTFilter);
	filteredCascades->SetBranchAddress("chi2AfterTFilter", 	   &current_ev.m_chi2AfterTFilter);
	filteredCascades->SetBranchAddress("energy", 			   &current_ev.m_energy);
	filteredCascades->SetBranchAddress("energySigma", 		   &current_ev.m_energySigma);
	filteredCascades->SetBranchAddress("theta",		 		   &current_ev.m_theta);
	filteredCascades->SetBranchAddress("thetaSigma",		   &current_ev.m_thetaSigma);
	filteredCascades->SetBranchAddress("phi", 				   &current_ev.m_phi);
	filteredCascades->SetBranchAddress("phiSigma", 			   &current_ev.m_phiSigma);
	filteredCascades->SetBranchAddress("directionSigma", 	   &current_ev.m_directionSigma);
	filteredCascades->SetBranchAddress("declination",		   &current_ev.m_declination);
	filteredCascades->SetBranchAddress("rightAscension",	   &current_ev.m_rightAscension);
	filteredCascades->SetBranchAddress("position", 			   &position);
	filteredCascades->SetBranchAddress("eventTime",			   &eventTime);
	filteredCascades->SetBranchAddress("time", 				   &current_ev.m_cascTime);
	filteredCascades->SetBranchAddress("mcEnergy", 			   &current_ev.m_mcEnergy);
	filteredCascades->SetBranchAddress("mcTheta", 			   &current_ev.m_mcTheta);
	filteredCascades->SetBranchAddress("mcPhi", 			   &current_ev.m_mcPhi);
	filteredCascades->SetBranchAddress("mcPosition", 		   &mcPosition);
	filteredCascades->SetBranchAddress("likelihood", 		   &current_ev.m_likelihood);
	filteredCascades->SetBranchAddress("likelihoodHitOnly",    &current_ev.m_likelihoodHitOnly);
	filteredCascades->SetBranchAddress("qTotal", 			   &current_ev.m_qTotal);
	filteredCascades->SetBranchAddress("nTrackHits", 		   &current_ev.m_nTrackHits);

	sortedEvents.reserve(filteredCascades->GetEntries());

	for (int i = 0; i < filteredCascades->GetEntries(); ++i)
	{
		filteredCascades->GetEntry(i);
		current_ev.m_position   = *position;
		current_ev.m_mcPosition = *mcPosition;
		current_ev.m_eventTime  = *eventTime;
		sortedEvents.push_back(current_ev);
	}

	sort(sortedEvents.begin(), sortedEvents.end(),Event::IsEarlier);

	double nSignal = 50;
	double nSignalSigma;

	SetFitter(1,false);

	int bins = 180;

	TH2F* h_nSignal  = new TH2F("h_nSignal","nSignal map",2*bins,-180,180,bins,-90,90);
	TH2F* h_testStat = new TH2F("h_testStat","Test statistic map",2*bins,-180,180,bins,-90,90);
	TH2F* h_prob     = new TH2F("h_prob","Probability map",2*bins,-180,180,bins,-90,90);

	for (int i = 0; i < (2*bins); ++i)
	{
		for (int j = 0; j < bins; ++j)
		{
			nSignal = 1;

			sigRa  = -180+360*i/(2*bins);
			sigDec = -90+180*j/bins;

			fit(nSignal,nSignalSigma);

			double ts = testStatistic(nSignal);
			double p;

			double dec   = (sigDec+90.0)/5.0;
			int low_dec  = floor(dec);
			int high_dec = ceil(dec);

			if(low_dec == high_dec)
			{
				if(ts <= fit_bounds[high_dec]) p = (*f_splines)[high_dec]->Eval(ts);
				else p = (*f_exps)[high_dec]->Eval(ts);
			}
			else
			{
				if(ts <= fit_bounds[high_dec])
					 p = (*f_splines)[low_dec]->Eval(ts)+(dec-low_dec)*((*f_splines)[high_dec]->Eval(ts)-(*f_splines)[low_dec]->Eval(ts));
				else p = (*f_exps)[low_dec]->Eval(ts)+(dec-low_dec)*((*f_exps)[high_dec]->Eval(ts)-(*f_exps)[low_dec]->Eval(ts));
			}

			// if(ts > 25) p = f_probfit2->Eval(ts);//cout << "p: " << p << endl;}
			// else p = f_probfit->Eval(ts);
			//cout << "p: " << p << ", ts: " << ts << ", dec: " << dec-low_dec << endl;
			h_nSignal->SetBinContent(i,j,nSignal);
			h_testStat->SetBinContent(i,j,ts);
			h_prob->SetBinContent(i,j,-log(p)/log(10));
		}
		cout << i << endl;
	}

	TCanvas* c1 = new TCanvas("c_nSig","nSignal map");
	gStyle->SetOptStat("0000");
	h_nSignal->Draw("z aitoff");
	h_nSignal->GetXaxis()->SetRangeUser(-200,200);
	h_nSignal->GetYaxis()->SetRangeUser(-105,105);
	// drawmap("");
	// drawLabels();

	TCanvas* c2 = new TCanvas("c_testStat","Test statistic map");
	h_testStat->Draw("z aitoff");

	TCanvas* c3 = new TCanvas("c_prob","Probability map");
	// gPad->SetLogz();
	h_prob->Draw("z aitoff");

	TFile* outputFile = new TFile("skyfit.root","RECREATE");
	c1->Write();
	c2->Write();
	c3->Write();

	return 0;
}

int main(int argc, char** argv)
{
	return skyfit();
}