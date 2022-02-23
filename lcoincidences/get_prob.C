#include <fstream>
#include <iostream>
#include <string>
#include <stdexcept>

using namespace std;

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

int get_prob()
{
	TH1F* prob_hist  = new TH1F("prob_hist","Test statistic distribution;Test statistic;Probability TS is bigger",10000,-1,22);
	TH1F* prob_hist2 = new TH1F("prob_hist2","nSignal distribution;nSignal;Probability nSign is bigger",10000,-1,5);

	ifstream  inf{"./data/merged/fdata_tStat_dec_-90.000000_all.txt"};
	ifstream inf2{"./data/merged/fdata_nSign_dec_-90.000000_all.txt"};

	if(!inf||!inf2) cout << "File not found.\n";

	while(inf)
	{
		string input;
		getline(inf,input);

		try 
		{
			if(input != "") prob_hist->Fill(stod(input));
    	} 
   
		catch (const std::invalid_argument&) 
		{
	        std::cerr << "Argument is invalid\n";
	        cerr << input << endl;
	        throw;
    	}
	}

	while(inf2)
	{
		string input;
		getline(inf2,input);

		try 
		{
			if(input != "") prob_hist2->Fill(stod(input));
    	} 
   
		catch (const std::invalid_argument&) 
		{
	        std::cerr << "Argument is invalid\n";
	        cerr << input << endl;
	        throw;
    	}
	}

	const int nodes = 50;
	const double low = -1;
	const double high = 33;
	double spline_nodes[nodes];

	for (int i = 0; i < nodes; ++i)
	{
		spline_nodes[i] = low+(high-low)*i/nodes;
	}

	NSpline<nodes>* s= new NSpline<nodes>(spline_nodes);
	s->SetDer1(0);
	s->SetDer2(0);
	TF1* f_spline = new TF1("f_probfit", s->GetEval(), low, high, 2*nodes+2);
	TF1* f_exp    = new TF1("f_probfit2", "[0]*exp(-[1]*x)",-1,50);
	f_exp->SetParameters(1,1);

	TCanvas* c1 = new TCanvas("c1","TS distribution");

	gStyle->SetOptStat(111111);

	gPad->SetLogy();
	TH1* prob_hist_cumul = prob_hist->GetCumulative(kFALSE);
	prob_hist_cumul->Scale(1./prob_hist->GetEntries());
	prob_hist_cumul->Fit(f_spline,"M","",low,high);
	prob_hist_cumul->Fit(f_exp,"M","",15,35);
	prob_hist_cumul->Draw();

	TCanvas* c2 = new TCanvas("c2","nS distribution");
	gPad->SetLogy();
	TH1* prob_hist2_cumul = prob_hist2->GetCumulative(kFALSE);
	prob_hist2_cumul->Scale(1./prob_hist2->GetEntries());
	prob_hist2_cumul->Draw();


	TFile* outputFile = new TFile("prob.root","RECREATE");
	prob_hist_cumul->Write();
	prob_hist2_cumul->Write();
	f_spline->Write();
	f_exp->Write();

	return 0;
}