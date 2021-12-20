#include "THStack.h"

#include "BRecoCascade.h"
#include "BJointHeader.h"
#include "BRunInfo.h"

#define DEBUG() cout << "Current Line: " << __LINE__ << endl;

#define PROFILLING 0
#define NEW_CASC_STRUCTURE
#include "../EventLoop.h"
#include "../Instrumentor.h"


using namespace std;

int cascade_flux(int data = 0, int year = -1, int cluster = -1)
{
#if PROFILLING
	Instrumentor::Get().BeginSession("Session Name");
#endif

	int LCut;
	long int maxTimeDiff;
	double energyCut = 200;

	cout << "Contained + 40 cut is applied.\n";
	cout << "Apply likelihood cut? [1/0]\n";
	cin >> LCut;

	bool DoLikelihoodCut = LCut == 1;
	if(DoLikelihoodCut) cout << "Likelihood (Lho <= 1.5) cut is applied.\n";

	cout << "\nMaximal time difference for coincidences in seconds: \n";
	cin >> maxTimeDiff;
	cout << "Maximal time difference set to " << maxTimeDiff << " seconds.\n";

	typedef std::function<void(const Event&)> FillFn;
	typedef std::function<void()> DrawFn;

	//make aitoff map for energyCut TeV+
	DrawSingle<TGraph>* aitoff = new DrawSingle<TGraph>();

	FillFn fillAitoff = [&energyCut,&aitoff](const Event& e)
	{
		if(e.m_energy >= energyCut)
		{
			double ra = e.m_rightAscension;
			if(ra > 180) ra -= 360;
			XY pos = toAitoff(ra,e.m_declination);
			// cout << "ra: " << radToDeg(rightAscension) << " dec: " << radToDeg(declination) << endl;
			// cout << "aitoff x: " << pos.x << " aitoff y: " << pos.y << endl;
			aitoff->drawsingle->SetPoint(aitoff->drawsingle->GetN(),pos.x,pos.y);
		}
	};
	aitoff->SetFillFunc(fillAitoff);

	DrawFn drawAitoff = [&energyCut,&aitoff]()
	{
		aitoff->canvas = new TCanvas("c_aitoff", "", 1000, 500);
		drawmap(Form("Cascades with energy > %f TeV;Right ascension;Declination", energyCut));
		aitoff->drawsingle->SetMarkerColor(kBlue);
		aitoff->drawsingle->SetMarkerSize(1);
		aitoff->drawsingle->SetMarkerStyle(4);
		aitoff->drawsingle->Draw("*");
		drawLabels();
	};
	aitoff->SetDrawFunc(drawAitoff);

	//making flux histograms
	DrawMap<TH1F>* flux_hist = new DrawMap<TH1F>();

	FillFn fillFlux = [&flux_hist](const Event& e)
	{
		//key for histogram identification
		TString hist_key = to_string(e.m_clusterID)+to_string(e.m_seasonID);

		//if histogram with given key does not exist, create one, fill histogram
		if(flux_hist->drawmap.find(hist_key)==flux_hist->drawmap.end())
		{
			TString hist_name = Form("h_cascFlux_y%dc%d",e.m_seasonID-2000,e.m_clusterID);
			flux_hist->drawmap[hist_key] = new TH1F(hist_name,hist_name(11,5)+"; Month; NoE [#]",50,GetStartTime(2016),GetEndTime(2016));

			flux_hist->drawmap[hist_key]->GetXaxis()->SetTimeDisplay(1);
			flux_hist->drawmap[hist_key]->GetXaxis()->SetTimeFormat("%m");//("%m/%Y");

			// int color = e.m_seasonID-2014+e.m_clusterID;
			// if(color > 9) color += 30;

			// flux_hist->drawmap[hist_key]->SetLineColor(color);
			flux_hist->drawmap[hist_key]->SetLineWidth(3);	
		}

		flux_hist->drawmap[hist_key]->Fill(e.m_eventTime.GetSec()-unix1995-GetStartTime(e.m_seasonID)+GetStartTime(2016)); //1970 unix to 1995 unix
	};
	flux_hist->SetFillFunc(fillFlux);

	DrawFn drawFlux = [&flux_hist]()
	{
		map<TString,THStack*> flux_stack;
		for(auto const& x : flux_hist->drawmap)
		{
			// flux_canv[x.first] = new TCanvas(x.first,"CascadeFlux",800,600);
			// x.second->Draw();

			TString season = x.first(1,4);
			TString cluster = x.first(0,1);

			//if THStack with given key does not exist, create one, fill THStack
			if(flux_stack.find(season)!=flux_stack.end()) flux_stack[season]->Add(x.second);
			else
			{
				TString stack_name = "hs_cascFlux_y"+season(2,2);
				flux_stack[season] = new THStack(stack_name,"Cascade flux over time; Month; NoE [#] per 7.3 days");

				flux_stack[season]->Add(x.second);
			}

			if(flux_stack.find(cluster)!=flux_stack.end()) flux_stack[cluster]->Add(x.second);
			else
			{
				TString stack_name = "hs_cascFlux_c"+cluster;
				flux_stack[cluster] = new THStack(stack_name,"Cascade flux over time; Month; NoE [#] per 7.3 days");

				flux_stack[cluster]->Add(x.second);
			}		
		}

		for(auto const& x : flux_stack)
		{
			flux_hist->canvasmap[x.first] = new TCanvas("c_cascFlux_" + x.first,"CascadeFlux",800,600);
			x.second->Draw("nostack PLC");

			x.second->GetXaxis()->SetTimeDisplay(1);
			x.second->GetXaxis()->SetTimeFormat("%m");
			x.second->Draw("nostack PLC");
	  		gPad->BuildLegend(0.75,0.75,0.95,0.95,"");	
		}
	};
	flux_hist->SetDrawFunc(drawFlux);

	EventLoop* eloop = new EventLoop(year,cluster);

	const char* data_path;

	switch(data)
	{
		case 0:
			data_path = "/home/vavrik/work/data";
			break;
		case 1:
			data_path = "/home/vavrik/bajkal/recoCascades/v1.2";
			break;
		case 2:
			data_path = "/media/vavrik/Alpha/BaikalData/dataGidra_v1.3";
			break;
		case 3:
			data_path = "/media/vavrik/Alpha/BaikalData/dataAries_v1.3";
			break;
		case 4:
			data_path = "/media/vavrik/Alpha/BaikalData/dataAries_v1.5";
			break;
	}

	eloop->LoadReco(data_path);

	string logs_path;
	if(data == 0) logs_path = "/home/vavrik/work/Baikal-GVD/cascade_flux/logs/programOutput_";
	if(data == 1) logs_path = "/home/vavrik/storage/casc_flux/logs/programOutput_";
	eloop->LoadRunLogs(logs_path);
	eloop->SetUpTTrees();

	eloop->UseLEDfilter();
	eloop->UseContainedFilter(40);
	if(DoLikelihoodCut) eloop->UseLikelihoodFilter(1.5);
	// eloop->UseEnergyFilter(energyCut);
	eloop->UseUpGoingFilter();

	eloop->AddDrawable(aitoff);
	eloop->AddDrawable(flux_hist);
	eloop->RunLoop();
	eloop->DrawAll();
	eloop->DrawFluxGraphs();
	eloop->SaveAll();

	// eloop->cfinder->FindTAEC(maxTimeDiff,360.0,20.0);
	// eloop->cfinder->WriteTAEC(maxTimeDiff,360.0,20.0);
	// eloop->cfinder->FindTAEC(maxTimeDiff,20.0,20.0);
	// eloop->cfinder->WriteTAEC(maxTimeDiff,20.0,20.0);
	// eloop->cfinder->FindTAEC(maxTimeDiff,10.0,20.0);
	// eloop->cfinder->WriteTAEC(maxTimeDiff,10.0,20.0);

	// eloop->cfinder->RandomCoincidences(maxTimeDiff);

	eloop->cfinder->WarnLEDMatrixRun();

	//fixed window events (turn on energyCut)
	double window = 6; //angle window in degrees
	int step = 3;

	TString title = Form("Events E > %f TeV per %f degree window;Right ascension;Declination",energyCut,window);
	TH2F* h_density = new TH2F("h_density",title,360/step,-180,180,180/step+1,-90,91);
	TH1F* h_density_stats = new TH1F("h_density_stats","Bins with given number of events;#NoE;#NoB",16,0,16);

	for(int ra = -180; ra < 181; ra += step)
	{
		for(int dec = -90; dec < 91; dec += step)
		{
			int noe = 0;

			for(auto ev : eloop->sortedEvents)
			{
				if(ev.angDist(ra,dec) < window)
				{
					h_density->Fill(ra,dec);
					noe++;
				}
			}
			h_density_stats->Fill(noe);
		}
	}

	TCanvas* c2 = new TCanvas("c2", "", 1000, 500);
	h_density->Draw("z aitoff");

	TCanvas* c3 = new TCanvas("c3", "", 1000, 500);
	h_density_stats->Draw();

	//for(auto ev : sortedEvents) cout << ev << "\n";

#if PROFILLING
	Instrumentor::Get().EndSession();
#endif

	return 0;
}

int main(int argc, char** argv) 
{
	int data, year, cluster;

	if(argc < 4) cluster = -1;
	else cluster = stoi(argv[3]);

	if(argc < 3) year = -1;
	else year = stoi(argv[2]);
	
	if(argc < 2) data = 0;
	else data = stoi(argv[1]);

	return cascade_flux(data,year,cluster);
}