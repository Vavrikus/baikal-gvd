#include "../EventLoop.h"

#include "TF1.h"
#include "TSpline.h"

double GetIQR(TH1* hist)
{
	double quartils[2];
	double probs[2] = {0.25,0.75};
	hist->GetQuantiles(2,quartils,probs);
	return quartils[1] - quartils[0];
}

//uses Freedman-Diaconis rule to rebin histogram
void OptimalBins(TH1F* hist, EventLoop* eloop)
{
	int oldbins = hist->GetNbinsX();
	int newbins = round(pow(hist->GetEntries(),1./3.)/GetIQR(hist));

	cout << "Old bins: " << oldbins << " New bins:" << newbins << endl;

	if(oldbins != newbins)
	{
		delete hist;
		hist = new TH1F("h_costheta","Filtered cascades cosine of theta distribution",newbins,-1.0,1.0);
		eloop->sortedEvents.clear();
		eloop->SetUpTTrees();
		eloop->RunLoop();
	}
}

Double_t spline_5nodes(Double_t *x, Double_t *par)
{
	/*Fit parameters:
	par[0-4]=X of nodes (to be fixed in the fit!)
	par[5-9]=Y of nodes
	par[10-11]=first derivative at begin and end (to be fixed in the fit!)
	*/

	par[0] = -1;
	par[1] = -0.64;
	par[2] = -0.3;
	par[3] = 0;
	par[4] = 0.5;
	// par[10] = 460;
	// par[11] = -40./3.;

	Double_t xx = x[0];

	Double_t xn[5] = { par[0], par[1], par[2], par[3], par[4] };
	Double_t yn[5] = { par[5], par[6], par[7], par[8], par[9] };

	Double_t b1 = par[10];
	Double_t e1 = par[11];

	TSpline3 sp3("sp3", xn, yn, 5, "b1e1", b1, e1);

	return sp3.Eval(xx);
}

template<int N>
class NSpline 
{
	typedef std::function<double(double*,double*)> EvalFn;

private:
	double node_pos[N];
	bool FixDer1 = false;
	bool FixDer2 = false;
	double der1,der2;

public:
	NSpline(const double node_pos[N])
	{
		for (int i = 0; i < N; ++i) this->node_pos[i] = node_pos[i];
	}

	void SetDer1(double der1) {this->der1 = der1; FixDer1 = true;}
	void SetDer2(double der2) {this->der2 = der2; FixDer2 = true;}

	double Eval(double* x, double* par)
	{
		/*Fit parameters:
		par[0-N-1]=X of nodes (to be fixed in the fit!)
		par[N-2N-1]=Y of nodes
		par[2N-2N+1]=first derivative at begin and end (to be fixed in the fit!)
		*/

		for (int i = 0; i < N; ++i) par[i] = node_pos[i];

		Double_t xx = x[0];

		Double_t xn[N];
		Double_t yn[N];

		for (int i = 0; i < N; ++i) xn[i] = par[i];
		for (int i = 0; i < N; ++i) yn[i] = par[N+i];

		if(FixDer1) par[2*N]   = der1;
		if(FixDer2) par[2*N+1] = der2;

		Double_t b1 = par[2*N];
		Double_t e1 = par[2*N+1];

		TSpline3 sp3("sp3", xn, yn, N, "b1e1", b1, e1);

		return sp3.Eval(xx);		
	}

	EvalFn GetEval()
	{
		return [this](double* a, double* b)
		{
			return this->Eval(a,b);
		};
	}

};

int data_pdf(int data = 0, int year = -1, int cluster = -1)
{
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
	eloop->SetUpTTrees();
	eloop->UseLEDfilter();
	eloop->UseContainedFilter(40);
	eloop->UseLikelihoodFilter(1.5);
	// eloop->UseEnergyFilter(40);

	typedef std::function<bool(const Event&)> FilterFn;

	FilterFn HitFilter = [](const Event& e)
	{
		if(e.m_nHitsAfterTFilter < 30) return false;
		else return true;
	};
	eloop->AddFilter(HitFilter);

	typedef std::function<void(const Event&)> FillFn;
	typedef std::function<void()> DrawFn;

	DrawSingle<TH1F>* costheta = new DrawSingle<TH1F>("h_costheta","Filtered cascades cosine of theta distribution",46,-1.0,1.0);

	FillFn fill_costheta = [&costheta](const Event& e)
	{
		costheta->drawsingle->Fill(cos(e.m_theta));
	};
	costheta->SetFillFunc(fill_costheta);

	eloop->AddDrawable(costheta);
	eloop->RunLoop();

	// OptimalBins(costheta->drawsingle,eloop);

	const int nodes = 10;
	const double low = -1;
	const double high = 1;
	const double spline_nodes[nodes] = {low,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,high};

	NSpline<nodes>* s= new NSpline<nodes>(spline_nodes);
	s->SetDer2(0);

	costheta->drawsingle->Scale(1./costheta->drawsingle->Integral());

	TF1* f_spline = new TF1("f_spline", s->GetEval(), low, high, 2*nodes+2);
	costheta->drawsingle->Fit(f_spline,"M","",low, high);
	costheta->drawsingle->Fit(f_spline,"L","",low, high);

	TCanvas* c1 = new TCanvas("c_costheta","c_costheta");
	costheta->drawsingle->Draw();

	TFile* outputFile = new TFile("cos_theta.root","RECREATE");
	costheta->drawsingle->Write();
	f_spline->Write();
	eloop->filteredCascades->Write();

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

	return data_pdf(data,year,cluster);
}