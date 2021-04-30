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

// std::vector<double> stringXPositions = {5.22,52.13,57.54,25.17,-29.84,-53.6,-42.32,0};
// std::vector<double> stringYPositions = {62.32,37.15,-13.92,-52.01,-52.36,-7.49,42.74,0};

double xPos[40] = {-13.76,32.14,45.06,5.13,-45.03,-76.21,-59.85,-14.47,-195.19,-164.79,-180.08,-227.51,-276.24,-279.59,-248.17,-222.70,-270.25,-228.58,-220.89,-261.89,-309.86,-337.48,-319.74,-282.27,65.85,108.73,113.87,74.19,25.1,-2.48,16.08,58.37,-163.91,-119.26,-113.90,-152.28,-202.59,-230.83,-213.25,-170.30};
double yPos[40] = {-211.35,-235.88,-285.45,-325.83,-319.82,-281.63,-231.37,-270.17,-340.62,-384.09,-435.13,-450.13,-424.31,-372.59,-337.03,-391.09,-37.36,-65.26,-117.78,-153.57,-146.26,-101.43,-55.24,-96.82,-435.47,-462.39,-514.68,-549.90,-544.25,-500.53,-453,-491.97,-628.26,-656.49,-707.52,-744.24,-738.58,-694.13,-645.06,-685.35};

vector<vector<vector<int>>> ledMatrixRuns = {{{2,3,4,5,6,7,8,9,10,11,118,119,177,193,194,200,201,228,229,230,231,232,233,234,235,236,237,560,598}},{{},{}},{{7,117,412,429,443,459,474,490,505,520,548,564,579,595},{1,2,3,6,7,37,134,428,450,464,480,495,510,527,540,568,584,599,615,631,647,668},{35,36,117,120,131,151,412,429,443,459,474,489,504,519,520,547,575,591,607,623,644}},{{17,18,37,38,39,40,44,61,77,93,97,111,126,142,158,174,190,203,218,232,247,264,277,292,362,377,392,407,422,437,452,467,484,536,551,566,583,596,611,628,644,661,676,677,693},{8,41,54,56,60,61,77,92,107,123,138,154,169,184,201,215,231,245,260,276,306,375,391,406,421,436,451,466,481,498,553,571,586,603,616,631,648,663,679,694,709},{8,9,10,24,80,93,109,124,139,155,170,185,201,216,233,247,262,276,291,329,330,331,337,406,422,437,453,468,483,498,513,530,594,595,596,597,611,612,629,642,657,674,689,705,720,735},{13,23,36,51,67,82,100,116,131,146,162,179,193,208,222,237,251,268,283,350,367},{13,23,34,50,67,82,86,88,89,90,91,92,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,112,113,114,116,117,118,120,121,122,123,124,129,130,132,137,147,163,180,193,208,222,237,238,253,265,279,363,379}},{{3,19,32,42,51,52,62,71,82,92,102,122,145,156,165,180},{12,14,24,33,35,42,51,60,69,83,90,111,134,145,146,147,155,156,157,158,159,160,162,164},{9,13,14,17,132,143,153,164,165,167,169,172},{1,15,17,21,26,36,46,58,67,76,86,94,103,112,114},{2,12,17,19,23,24,26,36,44,55,63,73,82,89,98,106,117,131,143,151,160,166,168,175},{18,20,25,31,41,51,62,71,90,110,118,128,130,143,145,154,157,163,166,173,178,185,195,220,232,241,250,260,282,296,301,312,326,336,346,356,367,384,394},{7,10,12,16,17,22,30,40,49,58,67,76,84,93,102,105,113,115,129,131,143,144,149,152,159,165,169,174,177,207,219,228,237,245,254,264,277,281,290,301}}};

TH1F* h_nHits = new TH1F("h_nHits","Number of hits per Event; N_{hits} [#]; NoE [#]",250,0,250);
TH1F* h_nHitsFull = new TH1F("h_nHitsFull","Number of hits per Event; N_{hits} [#]; NoE [#]",250,0,250);
TH1F* h_nHitsAfterTFilter = new TH1F("h_nHitsAfterTFilter","Number of hits per Event used for reconstruction; N_{hits} [#]; NoE [#]",250,0,250);
TH1F* h_nHitsAfterTFilterFull = new TH1F("h_nHitsAfterTFilterFull","Number of hits per Event used for reconstruction; N_{hits} [#]; NoE [#]",250,0,250);
TH1F* h_nStringsAfterTFilter = new TH1F("h_nStringsAfterTFilter","Number of strings per Event used for reconstruction; N_{strings} [#]; NoE [#]",10,0,10);
TH1F* h_nStringsAfterTFilterFull = new TH1F("h_nStringsAfterTFilterFull","Number of strings per Event used for reconstruction; N_{strings} [#]; NoE [#]",10,0,10);
TH1F* h_energy = new TH1F("h_energy","Reconstructed energy; E [TeV]; NoE [#]",1000,0,1000);
TH1F* h_energyFull = new TH1F("h_energyFull","Reconstructed energy; E [TeV]; NoE [#]",1000,0,1000);
TH1F* h_theta = new TH1F("h_theta","Zenith angle (0 = up-going, 180 down-going); #theta [deg] ",180,0,180);
TH1F* h_thetaFull = new TH1F("h_thetaFull","Zenith angle (0 = up-going, 180 down-going); #theta [deg] ",180,0,180);
TH1F* h_cosTheta = new TH1F("h_cosTheta","Cosine of the zenith angle; cos(#theta) [1]; NoE [#]",180,-1,1);
TH1F* h_phi = new TH1F("h_phi","Azimuth angle (0 = East, 90 = North); #phi [deg] ",360,0,360);
TH1F* h_phiFull = new TH1F("h_phiFull","Azimuth angle (0 = East, 90 = North); #phi [deg] ",360,0,360);
TH1F* h_thetaContained = new TH1F("h_thetaContained","Zenith angle, contained cascades (0 = up-going, 180 down-going); #theta [deg] ",180,0,180);
TH2F* h_energyNHits = new TH2F("h_energyNHits","Energy vs NHits; log_{10}(E [TeV]); N_{hits} [#]",50,0,5,100,0,100);
TH1F* h_qTotal = new TH1F("h_qTotal","Overall deposited charge;Q [p.e.];NoE [#]",1000,0,10000);
TH1F* h_qTotalFull = new TH1F("h_qTotalFull","Overall deposited charge;Q [p.e.];NoE [#]",1000,0,10000);
TH1F* h_likelihood = new TH1F("h_likelihood","Likelihood;L [#];NoE [#]",1000,0,10);
TH1F* h_likelihoodHitOnly = new TH1F("h_likelihoodHitOnly","Likelihood;L [#];NoE [#]",100,0,10);
TH1F* h_likelihoodFull = new TH1F("h_likelihoodFull","Likelihood;L [#];NoE [#]",100,0,10);
TH1F* h_chi2 = new TH1F("h_chi2","Chi2;#chi_{2} [1];NoE [#]",100,0,100);
TH1F* h_z = new TH1F("h_z","Z position; z [m]; NoE [#]",60,-300,300);
TH1F* h_horDist = new TH1F("h_horDist","Horizontal distance; #rho [m]; NoE [#]",100,0,100);
TH1F* h_nTrackHits = new TH1F("h_nTrackHits","Number of Track Hits; #N_{trackHits} [#]; NoE [#]",30,0,30);

TGraph* g_cascadeXY = new TGraph();
// TGraph* g_stringPositions = new TGraph(8,&stringXPositions[0],&stringYPositions[0]);
TGraph* g_stringPositions = new TGraph(40,xPos,yPos);
TH2F* h_distQTotal = new TH2F("h_distQTotal","Horizontal distance vs. Qtotal; #rho [m]; QTotal [p.e.]",100,0,100,100,0,3000);
TH2F* h_distEnergy = new TH2F("h_distEnergy","Horizontal distance vs. Energy; #rho [m]; E [TeV]",100,0,100,100,0,1000);
TH2F* h_zenithEnergy = new TH2F("h_zenithEnergy","Zenith vs. Energy; #theta [deq]; E [TeV]",180,0,180,100,0,1000);

TH1F* h_nHitsUpGoing = new TH1F("h_nHitsUpGoing","N_{hits} for up-going events (#theta < 80 deg.); N_{hits} [#]; NoE [#]",100,0,100);
TH1F* h_energyUpGoing = new TH1F("h_energyUpGoing","Energy for up-going events (#theta < 80 deg.); E [TeV] ; NoE [#]",200,0,200);

TH1F* h_eventRuns = new TH1F("h_eventRuns","Event Run ID; Run ID [#]; NoE [#]",50,0,1000);
TH1F* h_cascFlux;

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
	TCanvas* c_nHits = new TCanvas("c_nHits","NHits",800,600);
	h_nHits->Draw();

	TCanvas* c_nHitsAfterTFilter = new TCanvas("c_nHitsAfterTFilter","NHitsAfterTFilter",800,600);
	h_nHitsAfterTFilter->Draw();

	TCanvas* c_energy = new TCanvas("c_energy","Energy",800,600);
	h_energy->Draw();

	TCanvas* c_theta = new TCanvas("c_theta","Theta",800,600);
	h_theta->Draw();
	h_thetaContained->SetLineColor(kRed);
	h_thetaContained->Draw("same");

	TCanvas* c_cosTheta = new TCanvas("c_cosTheta","CosTheta",800,600);
	h_cosTheta->Draw();

	TCanvas* c_energyNHits2 = new TCanvas("c_energyNHits2","EnergyVsNHits2",800,600);
	h_energyNHits->Draw("colz");

	TCanvas* c_z = new TCanvas("c_z","Z position",800,600);
	h_z->Draw();

	TCanvas* c_horDist = new TCanvas("c_horDist","HorizontalDistance",800,600);
	h_horDist->Draw();

	TCanvas* c_distQTotal = new TCanvas("c_distQTotal","DistQTotal",800,600);
	h_distQTotal->Draw("colz");

	TCanvas* c_distEnergy = new TCanvas("c_distEnergy","DistEnergy",800,600);
	h_distEnergy->Draw("colz");

	TCanvas* c_zenithEnergy = new TCanvas("c_zenithEnergy","ZenithEnergy",800,600);
	h_zenithEnergy->Draw("colz");

	TCanvas* c_cascadePositions = new TCanvas("c_cascadePositions","Results",800,600);

	g_stringPositions->SetMarkerStyle(20);
	g_stringPositions->SetMarkerColor(kRed);
	g_stringPositions->Draw("AP");
	g_stringPositions->SetTitle("Cascade XY positions; X [m]; Y [m]");
	g_cascadeXY->SetMarkerStyle(21);
	g_cascadeXY->SetTitle("Positions of reconstructed cascades;X [m];Y [m]");
	g_cascadeXY->Draw("PSame");

	TCanvas* c_nHitsUpGoing = new TCanvas("c_nHitsUpGoing","NHitsUpGoing",800,600);
	h_nHitsUpGoing->Draw();

	TCanvas* c_energyUpGoing = new TCanvas("c_energyUpGoing","EnergyUpGoing",800,600);
	h_energyUpGoing->Draw();

	TCanvas* c_EventRuns = new TCanvas("c_EventRuns","EventRuns",800,600);
	h_eventRuns->Draw();

	TCanvas* c_cascFlux = new TCanvas("c_cascFlux","CascadeFlux",800,600);
	h_cascFlux->Draw();
}

void SaveResults(int year, int cluster)
{
	TString outputFileName = Form("mcResults_data_y%dc%d.root",year,cluster);
	TFile* outputFile = new TFile(outputFileName,"RECREATE");

	h_nHits->Write();
	h_nHitsAfterTFilter->Write();
	h_nStringsAfterTFilter->Write();
	h_energy->Write();
	h_theta->Write();
	h_cosTheta->Write();
	h_phi->Write();
	h_qTotal->Write();
	h_likelihood->Write();
	h_likelihoodHitOnly->Write();
	h_chi2->Write();
	h_z->Write();
	h_horDist->Write();
	h_nTrackHits->Write();
	h_distQTotal->Write();
	h_distEnergy->Write();
	h_zenithEnergy->Write();

	h_nHitsFull->Write();
	h_nHitsAfterTFilterFull->Write();
	h_nStringsAfterTFilterFull->Write();
	h_energyFull->Write();
	h_thetaFull->Write();
	h_phiFull->Write();
	h_qTotalFull->Write();
	h_likelihoodFull->Write();

	h_nHitsUpGoing->Write();
	h_energyUpGoing->Write();
}

void SaveWikiResults(int year, int cluster, TTree* filteredCascades)
{
	TString outputFileName = Form("../../results/wiki_data_y%dc%d.wiki",year,cluster);
	ofstream wikiFile(outputFileName);

	wikiFile << "# Reconstructed Cascades - Season: " << year << " Cluster: " << cluster+1 << endl;
	wikiFile << "| CascadeID | RunID | EventID | No. Hits | No. Strings | Energy [TeV] | Zenith [deg] | Azimuth [deg] | x | y | z | Likelihood |" << endl;
	wikiFile << "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |" << endl;

	for (int i = 0; i < filteredCascades->GetEntries(); ++i)
	{
		filteredCascades->GetEntry(i);
		wikiFile << "| " << i << " | " << runID << " | " << eventID << " | " << nHitsAfterTFilter << " | " << nStringsAfterTFilter << " | " << energy << " | " << theta/TMath::Pi()*180 << " | " << phi/TMath::Pi()*180 << " | " << position->X() << " | " << position->Y() << " | " << position->Z() << " | " << likelihood << " |" << endl;
	}

	wikiFile.close();
}

void SaveInterestingCascades(int year, int cluster, TTree* filteredCascades)
{
	TString outputFileName = Form("../../results/recCasc_y%dc%d.txt",year,cluster);
	ofstream outputFile(outputFileName);

	for (int i = 0; i < filteredCascades->GetEntries(); ++i)
	{
		filteredCascades->GetEntry(i);
		outputFile << seasonID << "\t" << clusterID+1 << "\t" << runID << "\t" << eventID << "\t" << theta << "\t" << phi << "\t" << position->X() << "\t" << position->Y() << "\t" << position->Z() << "\t" << nHitsAfterTFilter << "\t" << energy << "\t" << likelihoodHitOnly << endl;
	}

	outputFile.close();
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

int DSRCprototype(int year, int cluster = -1, int folder = 0, bool upGoing = false, bool highEnergy = true)
{
	TChain reconstructedCascades("Tree/t_RecCasc");

	TString filesDir;

	//choosing data based on cluster and year
	int startID = cluster!=-1?cluster:0;
	int endID = cluster!=-1?cluster+1:10;

	int startSeason = year!=-1?year:16;
	int endSeason = year!=-1?year+1:19+1;

	//path to data folder
	const char* env_p = "/home/vavrik/work/data";

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

	//resizing buffer to ensure correct automatic range
	if(nRecCasc > 1000) TH1::SetDefaultBufferSize(nRecCasc);
	h_cascFlux = new TH1F("h_cascFlux","Cascade flux over time; Time; NoE [#]",50,0,0);

	int nProcessedEvents = 0;
	int nHighEnergyEvents = 0;

	for (int i = 0; i < reconstructedCascades.GetEntries(); ++i)
	{
		reconstructedCascades.GetEntry(i);
		// cout << runID << endl;

		//remove cascades from calibration
		if (IsLEDMatrixRun(seasonID-2000,clusterID,runID))
			continue;

			// cout << "Energy above 100 TeV - RunID: " << runID << " EventID: " << eventID << " E = " << energy << " L = " << likelihood << " S = " << directionSigma << " N = " << nHitsAfterTFilter << " T = " << theta/TMath::Pi()*180 << " P = " << phi/TMath::Pi()*180 << " (" << position->X() << "," << position->Y() << "," << position->Z() << ")" << endl;

		h_nHitsFull->Fill(nHits);
		h_nHitsAfterTFilterFull->Fill(nHitsAfterTFilter);
		h_nStringsAfterTFilterFull->Fill(nStringsAfterTFilter);
		h_energyFull->Fill(energy);
		h_thetaFull->Fill(theta/TMath::Pi()*180);
		h_phiFull->Fill(phi/TMath::Pi()*180);
		h_qTotalFull->Fill(qTotal);
		h_likelihoodFull->Fill(likelihood);

		//selecting only contained cascades with likelihoodHitOnly <= 1.5
		if (!IsContained(position) || likelihoodHitOnly > 1.5) //|| theta/TMath::Pi()*180 > 80)
		// if (!IsContained(position) || likelihoodHitOnly > 3)
		// if (!IsContained(position,40) || likelihoodHitOnly > 3 || nHitsAfterTFilter < 50)
		// if (!IsContained(position,40) || likelihoodHitOnly > 1.5 || position->Z() > 200)
		// if (!IsUncontained(position,60,100) || likelihoodHitOnly > 3)
			continue;

			//cout << "Up-going Event - SeasonID: " << seasonID << " ClusterID: " << clusterID << " RunID: " << runID << " EventID: " << eventID << " E = " << energy << " L = " << likelihood << " Lho = " << likelihoodHitOnly << " T = " << theta/TMath::Pi()*180 << " S = " << directionSigma << " N = " << nHitsAfterTFilter << " Nt = " << nTrackHits << endl;

		// if (nTrackHits > 10)
			// cout << seasonID << " ClusterID: " << clusterID << " RunID: " << runID << " EventID: " << eventID << endl;

		// if (directionSigma > 5 || energy < 10 || energySigma > 5)
			// continue;
		// if (mcEnergy < 20)
			// continue;

		// cout << runID << " " << eventID << endl;

		TVector3 cascDirRec(0,0,1);
		cascDirRec.SetTheta(theta);
		cascDirRec.SetPhi(phi);

		if (energy > 100 && highEnergy)
		{
			//cout << "Energy above 100 TeV - SeasonID: " << seasonID << " ClusterID: " << clusterID << " RunID: " << runID << " EventID: " << eventID << " E = " << energy << " L = " << likelihood << " Lho = " << likelihoodHitOnly << " S = " << directionSigma << " N = " << nHitsAfterTFilter << " T = " << theta/TMath::Pi()*180 << " P = " << phi/TMath::Pi()*180 << " Q = " << qTotal << endl;
			//cout << (*position).X() << " " << (*position).Y() << " " << (*position).Z() << endl;
		}
			g_cascadeXY->SetPoint(nHighEnergyEvents,position->X()+xPos[8*(clusterID+1)-1],position->Y()+yPos[8*(clusterID+1)-1]);
			nHighEnergyEvents++;

		// if (theta/TMath::Pi()*180 < 80 && upGoing && nHitsAfterTFilter > 20 && likelihoodHitOnly < 1.5)
		if (theta/TMath::Pi()*180 < 80 && upGoing && nHitsAfterTFilter > 20)
		{
			//cout << "Up-going Event - SeasonID: " << seasonID << " ClusterID: " << clusterID << " RunID: " << runID << " EventID: " << eventID << " E = " << energy << " L = " << likelihood << " Lho = " << likelihoodHitOnly << " T = " << theta/TMath::Pi()*180 << " S = " << directionSigma << " N = " << nHitsAfterTFilter << " Nt = " << nTrackHits << endl;
			//cout << (*position).X() << " " << (*position).Y() << " " << (*position).Z() << endl;
			h_nHitsUpGoing->Fill(nHitsAfterTFilter);
			h_energyUpGoing->Fill(energy);
		}

		double horizontalDist = TMath::Sqrt(TMath::Power(position->X(),2)+TMath::Power(position->Y(),2));

		// if (horizontalDist < 60)
		// {
		// 	cout << "Horizontal Distance - SeasonID: " << seasonID << " ClusterID: " << clusterID << " RunID: " << runID << " EventID: " << eventID << " E = " << energy << " L = " << likelihood << " Lho = " << likelihoodHitOnly << " S = " << directionSigma << " N = " << nHitsAfterTFilter << " T = " << theta/TMath::Pi()*180 << " P = " << phi/TMath::Pi()*180 << " Q = " << qTotal << endl;
		// 	cout << (*position).X() << " " << (*position).Y() << " " << (*position).Z() << endl;
		// 	g_cascadeXY->SetPoint(nHighEnergyEvents,position->X(),position->Y());
		// }

		// cout << i << " " << runID << " " << eventID << " " << theta  << " " << phi << " " << nHitsAfterTFilter << " " << nStringsAfterTFilter << " " << likelihood << " " << qTotal << endl;
		// cout << energy << " " << energySigma << " " << directionSigma << endl;

		h_nHits->Fill(nHits);
		h_nHitsAfterTFilter->Fill(nHitsAfterTFilter);
		h_nStringsAfterTFilter->Fill(nStringsAfterTFilter);
		h_energy->Fill(energy);
		h_theta->Fill(theta/TMath::Pi()*180);
		h_cosTheta->Fill(TMath::Cos(theta)*(-1));
		h_phi->Fill(phi/TMath::Pi()*180);
		h_qTotal->Fill(qTotal);
		h_likelihood->Fill(likelihood);
		h_likelihoodHitOnly->Fill(likelihoodHitOnly);
		h_chi2->Fill(chi2AfterTFilter);
		h_z->Fill(position->Z());

		nProcessedEvents++;
		filteredCascades->Fill();

		h_horDist->Fill(horizontalDist);
		h_nTrackHits->Fill(nTrackHits);
		h_distQTotal->Fill(horizontalDist,qTotal);
		h_distEnergy->Fill(horizontalDist,energy);
		h_zenithEnergy->Fill(theta/TMath::Pi()*180,energy);
		h_energyNHits->Fill(TMath::Log10(energy),nHitsAfterTFilter);
		h_eventRuns->Fill(runID);

		h_cascFlux->Fill(*eventTime-788918400); //1970 unix to 1995 unix
		h_cascFlux->GetXaxis()->SetTimeDisplay(1);
		h_cascFlux->GetXaxis()->SetTimeFormat("%m/%Y");

		if (IsContained(position))
		{
			h_thetaContained->Fill(theta/TMath::Pi()*180);
			// position->Print();
		}
	}

	DrawResults();
	SaveWikiResults(year,cluster,filteredCascades);
	SaveInterestingCascades(year,cluster,filteredCascades);
	SaveResults(year,cluster);

	cout << nProcessedEvents << endl;
	cout << nHighEnergyEvents << endl;
	TString outputFileName = Form("filteredCascades_y%dc%d.root",year,cluster);
	TFile *newFile = new TFile(outputFileName,"recreate");
	filteredCascades->Write();
	newFile->Close();

	return 0;
}