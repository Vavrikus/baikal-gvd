#include "TCanvas.h"
#include "TString.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TVector3.h"
#include "TMath.h"
#include "TSystem.h"
#include "TTimeStamp.h"
#include "TGraph.h"

#include <iostream>
#include <fstream>
#include <map>
#include <cmath>

using namespace std;

// std::vector<double> stringXPositions = {5.22,52.13,57.54,25.17,-29.84,-53.6,-42.32,0};
// std::vector<double> stringYPositions = {62.32,37.15,-13.92,-52.01,-52.36,-7.49,42.74,0};

double xPos[40] = {-13.76,32.14,45.06,5.13,-45.03,-76.21,-59.85,-14.47,-195.19,-164.79,-180.08,-227.51,-276.24,-279.59,-248.17,-222.70,-270.25,-228.58,-220.89,-261.89,-309.86,-337.48,-319.74,-282.27,65.85,108.73,113.87,74.19,25.1,-2.48,16.08,58.37,-163.91,-119.26,-113.90,-152.28,-202.59,-230.83,-213.25,-170.30};
double yPos[40] = {-211.35,-235.88,-285.45,-325.83,-319.82,-281.63,-231.37,-270.17,-340.62,-384.09,-435.13,-450.13,-424.31,-372.59,-337.03,-391.09,-37.36,-65.26,-117.78,-153.57,-146.26,-101.43,-55.24,-96.82,-435.47,-462.39,-514.68,-549.90,-544.25,-500.53,-453,-491.97,-628.26,-656.49,-707.52,-744.24,-738.58,-694.13,-645.06,-685.35};

vector<vector<vector<int>>> ledMatrixRuns = {{{2,3,4,5,6,7,8,9,10,11,118,119,177,193,194,200,201,228,229,230,231,232,233,234,235,236,237,560,598}},{{},{}},{{7,117,412,429,443,459,474,490,505,520,548,564,579,595},{1,2,3,6,7,37,134,428,450,464,480,495,510,527,540,568,584,599,615,631,647,668},{35,36,117,120,131,151,412,429,443,459,474,489,504,519,520,547,575,591,607,623,644}},{{17,18,37,38,39,40,44,61,77,93,97,111,126,142,158,174,190,203,218,232,247,264,277,292,362,377,392,407,422,437,452,467,484,536,551,566,583,596,611,628,644,661,676,677,693},{8,41,54,56,60,61,77,92,107,123,138,154,169,184,201,215,231,245,260,276,306,375,391,406,421,436,451,466,481,498,553,571,586,603,616,631,648,663,679,694,709},{8,9,10,24,80,93,109,124,139,155,170,185,201,216,233,247,262,276,291,329,330,331,337,406,422,437,453,468,483,498,513,530,594,595,596,597,611,612,629,642,657,674,689,705,720,735},{13,23,36,51,67,82,100,116,131,146,162,179,193,208,222,237,251,268,283,350,367},{13,23,34,50,67,82,86,88,89,90,91,92,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,112,113,114,116,117,118,120,121,122,123,124,129,130,132,137,147,163,180,193,208,222,237,238,253,265,279,363,379}},{{3,19,32,42,51,52,62,71,82,92,102,122,145,156,165,180},{12,14,24,33,35,42,51,60,69,83,90,111,134,145,146,147,155,156,157,158,159,160,162,164},{9,13,14,17,132,143,153,164,165,167,169,172},{1,15,17,21,26,36,46,58,67,76,86,94,103,112,114},{2,12,17,19,23,24,26,36,44,55,63,73,82,89,98,106,117,131,143,151,160,166,168,175},{18,20,25,31,41,51,62,71,90,110,118,128,130,143,145,154,157,163,166,173,178,185,195,220,232,241,250,260,282,296,301,312,326,336,346,356,367,384,394},{7,10,12,16,17,22,30,40,49,58,67,76,84,93,102,105,113,115,129,131,143,144,149,152,159,165,169,174,177,207,219,228,237,245,254,264,277,281,290,301}}};

int unix1995 = 788918400; 

map<TString,TH1F*> flux_hist;
map<TString,THStack*> flux_stack;
map<TString,TCanvas*> flux_canv;

int seasonID, clusterID, runID, eventID, nHits, nHitsAfterCaus, nHitsAfterTFilter, nStringsAfterCaus, nStringsAfterTFilter, nTrackHits;
double energy,theta,phi,mcEnergy,mcTheta,mcPhi;
double energySigma,thetaSigma,phiSigma,directionSigma;
double chi2AfterCaus, chi2AfterTFilter, cascTime, likelihood, likelihoodHitOnly, qTotal;
double rightAscension, declination;
TVector3* position = new TVector3();
TVector3* mcPosition = new TVector3();
TTimeStamp* eventTime = new TTimeStamp();

void DrawResults()
{
	for(auto const& x : flux_hist)
	{
		flux_canv[x.first] = new TCanvas(x.first,"CascadeFlux",800,600);
		x.second->Draw();

		TString season = x.first(1,4);
		TString cluster = x.first(0,1);

		//if THStack with given key does not exist, create one, fill THStack
		if(flux_stack.find(season)!=flux_stack.end()) flux_stack[season]->Add(x.second);
		else
		{
			TString stack_name = "hs_cascFlux_y"+season(2,2);
			flux_stack[season] = new THStack(stack_name,"Cascade flux over time; Month; NoE [#]");

			flux_stack[season]->Add(x.second);
		}

		if(flux_stack.find(cluster)!=flux_stack.end()) flux_stack[cluster]->Add(x.second);
		else
		{
			TString stack_name = "hs_cascFlux_c"+cluster;
			flux_stack[cluster] = new THStack(stack_name,"Cascade flux over time; Month; NoE [#]");

			flux_stack[cluster]->Add(x.second);
		}		
	}

	for(auto const& x : flux_stack)
	{
		flux_canv[x.first] = new TCanvas(x.first,"CascadeFlux",800,600);
		x.second->Draw("nostack");

		x.second->GetXaxis()->SetTimeDisplay(1);
		x.second->GetXaxis()->SetTimeFormat("%m");
		x.second->Draw("nostack");	
	}
}

void SaveResults(int year, int cluster)
{
	TString outputFileName = Form("cascFlux_y%dc%d.root",year,cluster);
	TFile* outputFile = new TFile(outputFileName,"RECREATE");
	//for(auto const& x : flux_hist) x.second->Write();
	for(auto const& x : flux_stack) x.second->Write();
}

bool IsContained(TVector3* position, double distFromCluster = 0)
{
	if (TMath::Sqrt(TMath::Power(position->X(),2)+TMath::Power(position->Y(),2)) < 60+distFromCluster && TMath::Abs(position->Z()) < 265+distFromCluster)
		return true;
	else
		return false;
}

bool IsUncontained(TVector3* position, double near, double far)
{
	double horizontalDist = TMath::Sqrt(TMath::Power(position->X(),2)+TMath::Power(position->Y(),2));
	double verticalDist = TMath::Abs(position->Z());
	if ((horizontalDist < far && horizontalDist > near && verticalDist < 263) || (horizontalDist < far && verticalDist < 263+(far-60) && verticalDist > 263+(near-60)))
		return true;
	else
		return false;
}

bool IsLEDMatrixRun(int year, int cluster, int run)
{
	bool isLEDMatrixRun = false;
	for (int i = 0; i < ledMatrixRuns[year-16][cluster].size(); ++i)
	{
		if (run == ledMatrixRuns[year-16][cluster][i])
		{
			isLEDMatrixRun = true;
			break;
		}
	}
	return isLEDMatrixRun;
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
	return 1459468800+(season-2016)*(365+GetLeapYears(season))*86400-unix1995;
}

//returns unix (1995) time of 31/03/YYYY+1 23:59:59
int GetEndTime(int season)
{
	return 1491004799+(season-2016)*(365+GetLeapYears(season))*86400-unix1995;
}

int cascade_flux(bool val = false, int year = -1, int cluster = -1)
{
	TChain reconstructedCascades("Tree/t_RecCasc");

	TString filesDir;

	//choosing data based on cluster and year
	int startID = cluster!=-1?cluster:0;
	int endID = cluster!=-1?cluster+1:10;

	int startSeason = year!=-1?year:16;
	int endSeason = year!=-1?year+1:20+1;

	//path to data folder
	const char* env_p = val?"/home/vavrik/bajkal/recoCascades/v1.2":"/home/vavrik/work/data";

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

	TTree* filteredCascades = new TTree("filteredCascades","Filtered Cascades");

	reconstructedCascades.SetBranchAddress("seasonID", &seasonID);
	reconstructedCascades.SetBranchAddress("clusterID", &clusterID);
	reconstructedCascades.SetBranchAddress("runID", &runID);
	reconstructedCascades.SetBranchAddress("eventID", &eventID);
	reconstructedCascades.SetBranchAddress("nHits", &nHits);
	reconstructedCascades.SetBranchAddress("nHitsAfterCaus", &nHitsAfterCaus);
	reconstructedCascades.SetBranchAddress("nStringsAfterCaus", &nStringsAfterCaus);
	reconstructedCascades.SetBranchAddress("chi2AfterCaus", &chi2AfterCaus);
	reconstructedCascades.SetBranchAddress("nHitsAfterTFilter", &nHitsAfterTFilter);
	reconstructedCascades.SetBranchAddress("nStringsAfterTFilter", &nStringsAfterTFilter);
	reconstructedCascades.SetBranchAddress("chi2AfterTFilter", &chi2AfterTFilter);
	reconstructedCascades.SetBranchAddress("energy", &energy);
	reconstructedCascades.SetBranchAddress("energySigma", &energySigma);
	reconstructedCascades.SetBranchAddress("theta", &theta);
	reconstructedCascades.SetBranchAddress("thetaSigma", &thetaSigma);
	reconstructedCascades.SetBranchAddress("phi", &phi);
	reconstructedCascades.SetBranchAddress("phiSigma", &phiSigma);
	reconstructedCascades.SetBranchAddress("directionSigma", &directionSigma);
	reconstructedCascades.SetBranchAddress("declination",&declination);
	reconstructedCascades.SetBranchAddress("rightAscension",&rightAscension);
	reconstructedCascades.SetBranchAddress("position", &position);
	reconstructedCascades.SetBranchAddress("eventTime",&eventTime);
	reconstructedCascades.SetBranchAddress("time", &cascTime);
	reconstructedCascades.SetBranchAddress("mcEnergy", &mcEnergy);
	reconstructedCascades.SetBranchAddress("mcTheta", &mcTheta);
	reconstructedCascades.SetBranchAddress("mcPhi", &mcPhi);
	reconstructedCascades.SetBranchAddress("mcPosition", &mcPosition);
	reconstructedCascades.SetBranchAddress("likelihood", &likelihood);
	reconstructedCascades.SetBranchAddress("likelihoodHitOnly", &likelihoodHitOnly);
	reconstructedCascades.SetBranchAddress("qTotal", &qTotal);
	reconstructedCascades.SetBranchAddress("nTrackHits", &nTrackHits);

	filteredCascades->Branch("seasonID", &seasonID);
	filteredCascades->Branch("clusterID", &clusterID);
	filteredCascades->Branch("runID", &runID);
	filteredCascades->Branch("eventID", &eventID);
	filteredCascades->Branch("nHits", &nHits);
	filteredCascades->Branch("nHitsAfterCaus", &nHitsAfterCaus);
	filteredCascades->Branch("nStringsAfterCaus", &nStringsAfterCaus);
	filteredCascades->Branch("chi2AfterCaus", &chi2AfterCaus);
	filteredCascades->Branch("nHitsAfterTFilter", &nHitsAfterTFilter);
	filteredCascades->Branch("nStringsAfterTFilter", &nStringsAfterTFilter);
	filteredCascades->Branch("chi2AfterTFilter", &chi2AfterTFilter);
	filteredCascades->Branch("energy", &energy);
	filteredCascades->Branch("energySigma", &energySigma);
	filteredCascades->Branch("theta", &theta);
	filteredCascades->Branch("thetaSigma", &thetaSigma);
	filteredCascades->Branch("phi", &phi);
	filteredCascades->Branch("phiSigma", &phiSigma);
	filteredCascades->Branch("directionSigma", &directionSigma);
	filteredCascades->Branch("declination",&declination);
	filteredCascades->Branch("rightAscension",&rightAscension);
	filteredCascades->Branch("position", &position);
	filteredCascades->Branch("eventTime","TTimeStamp",&eventTime);
	filteredCascades->Branch("time", &cascTime);
	filteredCascades->Branch("mcEnergy", &mcEnergy);
	filteredCascades->Branch("mcTheta", &mcTheta);
	filteredCascades->Branch("mcPhi", &mcPhi);
	filteredCascades->Branch("mcPosition", &mcPosition);
	filteredCascades->Branch("likelihood", &likelihood);
	filteredCascades->Branch("likelihoodHitOnly", &likelihoodHitOnly);
	filteredCascades->Branch("qTotal", &qTotal);
	filteredCascades->Branch("nTrackHits", &nTrackHits);


	int nRecCasc = reconstructedCascades.GetEntries();

	cout << nRecCasc << endl;

	int nProcessedEvents = 0;

	for (int i = 0; i < reconstructedCascades.GetEntries(); ++i)
	{
		reconstructedCascades.GetEntry(i);

		//remove cascades from calibration
		if (IsLEDMatrixRun(seasonID-2000,clusterID,runID))
			continue;

		//selecting only contained cascades with likelihoodHitOnly <= 1.5
		if (!IsContained(position) || likelihoodHitOnly > 1.5) //|| theta/TMath::Pi()*180 > 80)
		// if (!IsContained(position) || likelihoodHitOnly > 3)
		// if (!IsContained(position,40) || likelihoodHitOnly > 3 || nHitsAfterTFilter < 50)
		// if (!IsContained(position,40) || likelihoodHitOnly > 1.5 || position->Z() > 200)
		// if (!IsUncontained(position,60,100) || likelihoodHitOnly > 3)
			continue;

		nProcessedEvents++;
		filteredCascades->Fill();

		//key for histogram identification
		TString hist_key = to_string(clusterID)+to_string(seasonID);

		//event before 01/01/2016 warning
		if(eventTime->GetSec() < 1451606400)
			cout << "Event " << eventID << " has low eventTime: " << eventTime << " runID: " << runID << "\n";

		//if histogram with given key does not exist, create one, fill histogram
		if(flux_hist.find(hist_key)!=flux_hist.end())
		{
			flux_hist[hist_key]->Fill(*eventTime-unix1995-GetStartTime(seasonID)+GetStartTime(2016)); //1970 unix to 1995 unix
		}

		else
		{
			TString hist_name = Form("h_cascFlux_y%dc%d",seasonID-2000,clusterID);
			flux_hist[hist_key] = new TH1F(hist_name,"Cascade flux over time; Month; NoE [#]",50,GetStartTime(2016),GetEndTime(2016));

			flux_hist[hist_key]->GetXaxis()->SetTimeDisplay(1);
			flux_hist[hist_key]->GetXaxis()->SetTimeFormat("%m");//("%m/%Y");
			flux_hist[hist_key]->SetLineColor(seasonID-2010+clusterID);	

			flux_hist[hist_key]->Fill(*eventTime-unix1995-GetStartTime(seasonID)+GetStartTime(2016)); //1970 unix to 1995 unix

			cout << "Year: " << seasonID << " Cluster: " << clusterID << "\n";
		}
	}

	gStyle->SetOptStat(111111);

	if(!val) DrawResults();
	SaveResults(year,cluster);

	cout << nProcessedEvents << endl;
	TString outputFileName = Form("filteredCascades_y%dc%d.root",year,cluster);
	TFile *newFile = new TFile(outputFileName,"recreate");
	filteredCascades->Write();
	newFile->Close();

	return 0;
}

int main(int argc, char** argv) 
{
	int val, year, cluster;

	if(argc < 3) cluster = -1;
	else cluster = stoi(argv[3]);

	if(argc < 2) year = -1;
	else year = stoi(argv[2]);
	
	if(argc < 1) val = 0;
	else val = stoi(argv[1]);

	cascade_flux(val,year,cluster);
	return 0;
}