#include "transformations.h"
#include "MCevents.h"
#include "threading.h"

#include "TRandom2.h"
#include "TF1.h"
#include "TTree.h"
#include "TFile.h"

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
			double energy   = energyDist->GetRandom();
			double altitude = altDist->GetRandom();
			double azimuth  = 360*rnd.Rndm();
	
			eqCoor pos 		= horToEq(horCoor(altitude,azimuth,uTime+time));
	
			//saving only events close to signal
			if(true)//angularDistance(ra,dec,pos.rAsc,pos.dec)<numOfSigma*posSigma)
			{
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
	std::vector<MCEvent>* mixed  = getBackround(MJDstart, timeWindow, numBackround, signalMean, signalSigma,
												ra, dec, posSigma, 7);
	std::vector<MCEvent>* signal = getSignalGaus(numSignal, ra, dec, posSigma, signalMean, signalSigma);

	*mixed += *signal;
	std::sort(mixed->begin(), mixed->end(), isSooner);

	return mixed;
}

void MCevents(const char* outpath = "experiments_data/MCevents.root", int back = 500, int sig = 0,
			  double timeWindow = 365.25, double timeSigma = 0.1, double posSigma = 5)
{
	TTree MC_events("MC_events","MC_events");
	
	MCEvent event{};
	auto branch1 = MC_events.Branch("type",   &event.type);
	auto branch2 = MC_events.Branch("mjd", 	  &event.mjd);
	auto branch3 = MC_events.Branch("ra", 	  &event.ra);
	auto branch4 = MC_events.Branch("dec",	  &event.dec);
	auto branch5 = MC_events.Branch("energy", &event.energy);

	std::vector<MCEvent>* b = getMixed(back,sig,58000-timeWindow/2,timeWindow,58000,timeSigma,posSigma,0,0);

	for(MCEvent ev : *b)
	{
		event = ev;
		MC_events.Fill();
	}

	TFile outfile(outpath,"RECREATE");

	MC_events.Write();
	outfile.Close();
}

// ./MCevents outpath back sig timeWindow timeSigma posSigma
int main(int argc, char** argv)
{
	if(argc > 1) MCevents(argv[1], std::stoi(argv[2]), std::stoi(argv[3]), std::stod(argv[4]),
	 std::stod(argv[5]), std::stod(argv[6]));
	else MCevents();

	return 0;
}

//plotting distributions
/*
	//std::vector<MCEvent>* b = getBackround(1,100000,58000);
	//for(auto x : *b) std::cout << x << '\n';

	double ra = 0;
	double dec = 89.999;

	std::vector<MCEvent>* b = getMixed(10000,1000,58000,1,58000.5,0.1,5,0,0);

	TCanvas* c1 = new TCanvas("","",1000,500);
	c1->cd();

	drawmap("Spatial distrubution");
	TGraph* gBack = new TGraph();

	for(auto x : *b)
	{
		double rAsc = x.ra;
        if(rAsc>180) rAsc -= 360;
        XY bCoor = toAitoff(rAsc, x.dec);

        gBack->SetPoint(gBack->GetN(),bCoor.x, bCoor.y);
	}

	gBack->SetMarkerColor(kBlue);
	gBack->Draw("P");
	drawLabels();

	TCanvas* c2 = new TCanvas();
	c2->cd();

	int bins = std::pow(b->size(),0.5);
	TH1F* h = new TH1F("Background","Time distrubution", bins, 0, 1);
	for(auto x : *b) h->Fill(x.mjd-58000);
	h->Scale(bins/(double)b->size());
	h->Draw("hist");

	TCanvas* c3 = new TCanvas();
	c3->cd();

	TH1F* h2 = new TH1F("Signal", "Distance from source", bins, 0, 6*10);
	for(auto x : *b)
	{
		h2->Fill(angularDistance(x.ra,x.dec,ra,dec));
	}
	h2->Scale(bins/(double)b->size());
	h2->Draw("hist");

*/


//energy distribution checking
/*
	//energy limits
	double minEnergy = 1;
	double maxEnergy = 10000;

	//energy distrubution for stmosferic neutrinos
	TF1* energyDist = new TF1("energyDist", "[0]*std::pow(x,[1])", minEnergy, maxEnergy);

	//normalization
	double power = -3.7; //-3.7 for atmospheric neutrinos
	double area  = powerLawArea(power, minEnergy, maxEnergy);

	energyDist->SetParameter(0,1/area);
	energyDist->SetParameter(1,power);

	int samples = 1000000;
	int bins = std::pow(samples,0.5);
	double binWidth = (maxEnergy-minEnergy)/bins;

	TH1F* h = new TH1F("Atmosferic neutrinos","Atmosferic neutrinos",bins, minEnergy, maxEnergy);
	for(int i = 0; i < samples; i++) h->Fill(energyDist->GetRandom());
	h->Scale(1./(samples*binWidth));
	
	TCanvas* c1 = new TCanvas();
	c1->cd();

	h->Draw("hist");
	energyDist->Draw("same");

	gPad->SetLogx();
	gPad->SetLogy();
*/