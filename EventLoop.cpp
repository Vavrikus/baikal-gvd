#include <deque>
#include <fstream>

#include "TFile.h"
#include "TSystem.h"

#include "Instrumentor.h"
#include "profilling.h"

std::ostream& operator<<(std::ostream& stream, const Event& ev)
{
	stream << "     =====================================  EVENT INFO  =====================================\n";
	stream << "     Time Stamp:\n          ";
	ev.m_eventTime.Print();

	stream << "\n     IDs:\n     ";
    stream << "     eventID: " << ev.m_eventID << ", coincidenceID: " << ev.m_coincidenceID;
    stream << ", seasonID: " << ev.m_seasonID << ", clusterID: " << ev.m_clusterID;
    stream << ", runID: " << ev.m_runID << "\n\n";

    stream << "     Reconstructed variables:\n     ";
    stream << "     energy = " << ev.m_energy << " TeV, sigma = " << ev.m_energySigma << " TeV\n     ";
    stream << "     theta = " << ev.m_theta/TMath::Pi()*180 << ", sigma = " << ev.m_thetaSigma/TMath::Pi()*180;
    stream << ", phi = " << ev.m_phi/TMath::Pi()*180 << ", sigma = " << ev.m_phiSigma/TMath::Pi()*180;
    stream << ", direction sigma = " << ev.m_directionSigma << "\n     ";
    stream << "     Position (XYZ): " << ev.m_position.X();
    stream << " " << ev.m_position.Y() << " " << ev.m_position.Z() << "\n     ";
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
    stream << "     Position (XYZ): " << ev.m_mcPosition.X();
    stream << " " << ev.m_mcPosition.Y() << " " << ev.m_mcPosition.Z() << endl;

    return stream;
}

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

std::ostream& operator<<(std::ostream& stream, const CoincidenceFinder::Coincidence& c)
{
	stream << "#######################################  COINCIDENCE  #######################################\n";
	stream << "ID: " << c.m_id << "\n";
	stream << "Number of events: " << c.m_indexes.size() << "\n";

	stream << "Time differences: ";

	for(int i = 1; i < c.m_indexes.size(); i++)
	{
		if(i!=1) cout << ", ";
		stream <<  c.cfinder->sortedEvents[c.m_indexes[i]].m_eventTime.GetSec()
				  -c.cfinder->sortedEvents[c.m_indexes[i-1]].m_eventTime.GetSec(); 
	}

	stream << "\nMinimal angular distance: ";
	double minAngDist = c.angDist(0,1);

	for (int i = 0; i < c.m_indexes.size(); i++)
	{
		for(int j = i+1; j < c.m_indexes.size(); j++)
		{
			if(c.angDist(i,j) < minAngDist)
			{
				minAngDist = c.angDist(i,j);
			}
		}
	}

	stream << minAngDist << "\n\n";

	for(int ev_index : c.m_indexes) stream << c.cfinder->sortedEvents[ev_index] << "\n" << endl;

    return stream;
}

//returns if events with indexes i,j are part coincidence with given time and angle differences
//and if both events have greater energy (TeV) than given 
bool CoincidenceFinder::IsTAEC(int i, int j, long int maxTimeDiff, double maxAngDist = 360, double minEnergy = 0)
{
	bool timeOK   = sortedEvents[j].m_eventTime.GetSec() - sortedEvents[i].m_eventTime.GetSec() <= maxTimeDiff;
	bool angleOK  = sortedEvents[i].angDist(sortedEvents[j]) <= maxAngDist;
	bool energyOK = (sortedEvents[i].m_energy >= minEnergy) and (sortedEvents[j].m_energy >= minEnergy);
	
	return timeOK and angleOK and energyOK;
}

//returns if events with indexes i,j satisfy time, run and position criterions
bool CoincidenceFinder::IsTRPC(int i, int j, long int maxTimeDiff = 3600, double maxDist = 5)
{
	bool timeOK = sortedEvents[j].m_eventTime.GetSec() - sortedEvents[i].m_eventTime.GetSec() <= maxTimeDiff;
	bool posOK  = sortedEvents[i].Dist(sortedEvents[j]) <= maxDist;
	bool yearOK	= sortedEvents[i].m_seasonID  == sortedEvents[j].m_seasonID;
	bool clusOK = sortedEvents[i].m_clusterID == sortedEvents[j].m_clusterID;
	bool runOK  = sortedEvents[i].m_runID	  == sortedEvents[j].m_runID;

	return timeOK and posOK and yearOK and clusOK and runOK;
}

//outputs coincidence counts into console
void CoincidenceFinder::WriteCStats()
{
	vector<int> counts; //how many coincidences with given number of events, 0th index = 2

	for(Coincidence* c : coincidences)
	{
		int cEvents = c->m_indexes.size();

		if(cEvents-1 > counts.size())
			for (int i = counts.size(); i < cEvents-1; ++i) counts.push_back(0);

		counts[cEvents-2]++;
	}

	cout << "#NoE   #NoC\n";
	cout << "=============\n";
	for (int i = 0; i < counts.size(); ++i)
	{
		random_coincidences->Fill(counts[i],i+2);
		if(i < 8)       cout << i+2 << "      " << counts[i] << "\n";
		else if(i < 98) cout << i+2 << "     " << counts[i] << "\n";
		else            cout << i+2 << "    " << counts[i] << "\n";
	}
	cout << "\nTotal coincidences: " << coincidences.size() << "\n\n";
}

//writes coincidences with stats into console
void CoincidenceFinder::WriteTAEC(long int maxTimeDiff, double maxAngDist, double minEnergy)
{	
	cout << "\n\nCoincidences with maximal time difference " << maxTimeDiff;
	cout << " seconds, maximal distance " << maxAngDist << " degrees\n";
	cout << "and minimal energy " << minEnergy << " TeV:\n" << endl;
	for(Coincidence* c : coincidences) cout << *c;

	WriteCStats();
}

//writes warning if two or more cascades are separated by smaller than selected amount of time
//and smaler than selected angle, saves coincidences into a vector
template<typename... args>
void CoincidenceFinder::FindCoincidences(bool(CoincidenceFinder::*IsCoin)(int, int, long int, args...), uint random_offset, long int maxdt, args... a)
{
	int numOfCoincidences = 0; //number of created coincidences, may be bigger than actual count

	coincidences.clear();

	//reset coincidence IDs
	for(int i = 0; i < sortedEvents.size(); ++i) sortedEvents[i].m_coincidenceID = -1;

	if(sortedEvents.size() > 0)
	{
		for(int i = 0; i < sortedEvents.size()-1; ++i)
		{
			Coincidence* c = new Coincidence(this);

			bool IsCoincidence = false;
			int currentID = sortedEvents[i].m_coincidenceID;
			if(currentID != -1) continue;

			int startEvent;

			//if random time shift is applied, find first later event and set it as start 
			if(random_offset != 0)
			{
				//add random shift
				//cout << "old " << i << " " << sortedEvents[i].m_eventTime.GetSec() << "\n";
				sortedEvents[i].m_eventTime.SetSec(sortedEvents[i].m_eventTime.GetSec()+random_offset);
				//cout << "new " << i << " " << sortedEvents[i].m_eventTime.GetSec() << "\n";

				int lower_bound  = i;
				int higher_bound = sortedEvents.size()-1;

				while(lower_bound != higher_bound-1)
				{
					int mid = floor((higher_bound+lower_bound)/2.0);
					//cout << lower_bound << " " << mid << " " << higher_bound << endl;
					if(sortedEvents[i].m_eventTime.GetSec() > sortedEvents[mid].m_eventTime.GetSec())
						lower_bound = mid;
					else higher_bound = mid;
				}

				//if should not happen for good shift
				if(higher_bound == i) 
				{
					startEvent = i+1;
					cout << "WARNING: Random shift insufficient.";
				}

				else startEvent = higher_bound;
				//cout << "startEvent " << startEvent << " " << sortedEvents[startEvent].m_eventTime.GetSec() << "\n";
			}

			else startEvent = i+1;

			//searching events coinciding with event i
			for (int j = startEvent; j < sortedEvents.size(); ++j)
			{
				if(sortedEvents[j].m_eventTime.GetSec()-maxdt > sortedEvents[i].m_eventTime.GetSec())
					break;

				if((this->*IsCoin)(i,j,maxdt,a...) and (sortedEvents[j].m_coincidenceID == -1))
				{	
					if(!IsCoincidence)
					{
						IsCoincidence = true;

						c->m_indexes.clear();
						c->m_id = numOfCoincidences;
						c->AddEvent(i);

						numOfCoincidences++;
					}

					c->AddEvent(j);					
				}
			}

			//adding new coincidence
			if(IsCoincidence and (currentID == -1)) coincidences.push_back(c);

			if(random_offset != 0) //remove random shift
			{
				//cout << "old2 " << i << " " << sortedEvents[i].m_eventTime.GetSec() << "\n";				
				sortedEvents[i].m_eventTime.SetSec(sortedEvents[i].m_eventTime.GetSec()-random_offset);
				//cout << "new2 " << i << " " << sortedEvents[i].m_eventTime.GetSec() << "\n";
			}
		}
	}


	for(Coincidence* c : coincidences) sort(c->m_indexes.begin(),c->m_indexes.end());

	//if((void*)IsCoin == (void*)IsTAEC) WriteTAEC(maxdt,a...);
}

void CoincidenceFinder::WarnLEDMatrixRun(int minCoinSize, long int maxTimeDiff, double maxDist)
{
	PROFILE_FUNCTION();

	cout << "\nPossible LED matrix runs detected:\n  coincidences with ";
	cout << minCoinSize << " or more events\n  position difference max ";
	cout << maxDist <<" meters\n  time difference max " << maxTimeDiff << " seconds\n";
	cout << "======================================\n";

	bool noLEDRunsDetected = true;

	FindCoincidences(&CoincidenceFinder::IsTRPC,0,maxTimeDiff,maxDist);

	for(int season = 2016; season < 2022; season++)
	{
		for(int cluster = 0; cluster < 10; cluster++)
		{
			vector<int> runs;

			for(Coincidence* c : coincidences)
			{
				if(std::find(runs.begin(), runs.end(), sortedEvents[c->m_indexes[0]].m_runID) == runs.end()) {
					if(sortedEvents[c->m_indexes[0]].m_seasonID == season && sortedEvents[c->m_indexes[0]].m_clusterID == cluster && c->m_indexes.size() >= minCoinSize)
					{
						noLEDRunsDetected = false;
						cout << "seasonID: "   << sortedEvents[c->m_indexes[0]].m_seasonID;
						cout << " clusterID: " << sortedEvents[c->m_indexes[0]].m_clusterID;
						cout << " runID: "     << sortedEvents[c->m_indexes[0]].m_runID << "\n";
						runs.push_back(sortedEvents[c->m_indexes[0]].m_runID);
					}
				}
			}
		}
	}

	if(noLEDRunsDetected) cout << "No runs detected." << endl;
	cout << endl;
}

void CoincidenceFinder::RandomCoincidences(long int maxTimeDiff,double maxAngDist)
{
	int iterations = 10000;

	for (int i = 0; i < iterations; ++i)
	{
		FindCoincidences(&CoincidenceFinder::IsTAEC,(i+1)*20,maxTimeDiff,maxAngDist,20.0);
	}

	random_coincidences->Fill(0.0,2.1,iterations-random_coincidences->GetEntries());
	random_coincidences->GetXaxis()->SetTitle("#NoC");
	random_coincidences->GetYaxis()->SetTitle("#NoE");

	TCanvas* c1 = new TCanvas("c1", "", 1000, 500);
	random_coincidences->Draw("Lego2");
}

#ifdef NEW_CASC_STRUCTURE
Event::Event(BRecoCascade* bcasc, BJointHeader* bhead)
{
	this->m_energy               = bcasc->GetEnergyRec();
	this->m_theta                = bcasc->GetThetaRec();
	this->m_phi                  = bcasc->GetPhiRec();
	this->m_rightAscension       = bcasc->GetSourceRARec();
	this->m_declination          = bcasc->GetSourceDecliRec();
	this->m_eventTime            = bhead->GetTimeCC();
	this->m_seasonID             = bhead->GetSeason();
	this->m_clusterID            = bhead->GetCluster();
	this->m_runID                = bhead->GetRun();
	this->m_eventID              = bhead->GetEventIDCC();
	this->m_nHits                = bcasc->GetNHits();
	this->m_nHitsAfterCaus       = bcasc->GetNHitsCaus();
	this->m_nHitsAfterTFilter    = bcasc->GetNHitsTFil();
	this->m_nStringsAfterCaus    = bcasc->GetNStringsCaus();
	this->m_nStringsAfterTFilter = bcasc->GetNStringsTFil();
	this->m_nTrackHits           = 0;//bcasc->Get();
	this->m_mcEnergy             = bcasc->GetEnergyMC();
	this->m_mcTheta              = bcasc->GetThetaMC();
	this->m_mcPhi                = bcasc->GetPhiMC();
	this->m_energySigma          = 0;//bcasc->Get();
	this->m_thetaSigma           = 0;//bcasc->Get();
	this->m_phiSigma             = 0;//bcasc->Get();
	this->m_directionSigma       = 0;//bcasc->Get();
	this->m_chi2AfterCaus        = 0;//bcasc->Get();
	this->m_chi2AfterTFilter     = 0;//bcasc->Get();
	this->m_cascTime             = bcasc->GetFitTime();
	this->m_likelihood           = bcasc->GetLikelihood();
	this->m_likelihoodHitOnly    = bcasc->GetLikelihoodHitOnly();
	this->m_qTotal               = bcasc->GetQTotal();
	this->m_position             = bcasc->GetFitPos();
	this->m_mcPosition           = bcasc->GetPosMC();
	this->m_distanceCS			 = bcasc->GetDistanceCS();
}
#endif

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

#ifndef NEW_CASC_STRUCTURE
	void EventLoop::LoadReco(const char* env_p, bool v379 = false)
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
		reconstructedCascades.SetBranchAddress("position", 			   &position);
		reconstructedCascades.SetBranchAddress("eventTime",			   &eventTime);
		reconstructedCascades.SetBranchAddress("time", 				   &current_ev.m_cascTime);
		reconstructedCascades.SetBranchAddress("mcEnergy", 			   &current_ev.m_mcEnergy);
		reconstructedCascades.SetBranchAddress("mcTheta", 			   &current_ev.m_mcTheta);
		reconstructedCascades.SetBranchAddress("mcPhi", 			   &current_ev.m_mcPhi);
		reconstructedCascades.SetBranchAddress("mcPosition", 		   &mcPosition);
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

#else
	void EventLoop::LoadReco(const char* env_p, bool v379 = false)
	{
		const char* versionFolder;
		if(v379) versionFolder = "v1.3-379";
		else     versionFolder = "v1.3-547";

		reconstructedCascades.SetName("Events");

		TString filesDir;

		for (int j = startSeason; j < endSeason; j++)
		{
			for (int i = startID; i < endID; ++i)
			{
				TString filesDir = Form("%s/exp20%d/%s/cluster%d/*.reco.cascade.root",env_p,j,versionFolder,i);
				cout << filesDir << "\n";
				reconstructedCascades.Add(filesDir);
			}
		}
	}

	void EventLoop::SetUpTTrees()
	{
		filteredCascades = new TTree("filteredCascades","Filtered Cascades");

		reconstructedCascades.SetBranchAddress("BRecoCascade.", &myCascade);
		reconstructedCascades.SetBranchAddress("BJointHeader.", &myHeader);

		filteredCascades->Branch("BRecoCascade.", &myCascade);
		filteredCascades->Branch("BJointHeader.", &myHeader);

		int nRecCasc = reconstructedCascades.GetEntries();

		cout << "\nnRecCasc: " << reconstructedCascades.GetEntries() << endl;

		sortedEvents.reserve(nRecCasc);
	}
#endif //NEW_CASC_STRUCTURE

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

void EventLoop::UseUpGoingFilter()
{
	FilterFn UpGoingFilter = [](const Event& e){if(e.m_theta > 90) return false; else return true;};
	filters.push_back(UpGoingFilter);
}

void EventLoop::RunLoop()
{
	for (int i = 0; i < reconstructedCascades.GetEntries(); ++i)
	{
		PrintProgress(i,reconstructedCascades.GetEntries());
		reconstructedCascades.GetEntry(i);
	#ifndef NEW_CASC_STRUCTURE
		current_ev.m_position   = *position;
		current_ev.m_mcPosition = *mcPosition;
		current_ev.m_eventTime  = *eventTime;
	#else
		current_ev = Event(myCascade,myHeader);
	#endif
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

	cfinder = new CoincidenceFinder(sortedEvents);
}

void EventLoop::SaveAll()
{
	TString outputFileName = Form("cascFlux_y%dc%d.root",year,cluster);
	TFile* outputFile = new TFile(outputFileName,"RECREATE");
	for(IDrawable* d : drawables) d->Save();
	for(auto x : flux_canv) x.second->Write();
	outputFile->Close();
	
	TString outputFileName2 = Form("filteredCascades_y%dc%d.root",year,cluster);
	TFile *newFile = new TFile(outputFileName2,"recreate");
	filteredCascades->Write();
	newFile->Close();		
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