#include <fstream>
#include <iostream>
#include <string>
#include <stdexcept>
#include <vector>

#include <TF1.h>
#include <TH1.h>
#include <TFile.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TPad.h>

using namespace std;

string data_folder = "./data/merged/";//"/media/vavrik/Elements/Baikal-GVD/lcoincidences/fit_results/test_03_dec_5degstep_full/"; //
double prob_bound = 1-1e-5;

//1e-5 probability bounds
const vector<double> fit_bounds = {10.9983,10.888,10.6235,10.4007,10.3206,10.3551,10.3797,10.5331,10.6906,10.8716,11.1583,11.6488,12.2753,13.1082,14.1282,15.2548,16.4592,17.9392,19.7688,22.1002,24.8501,28.3045,31.8804,35.5551,39.5829,44.0158,49.259,54.5009,58.8381,61.8463,63.9316,65.1841,66.6971,68.6378,70.9946,72.7542,73.3589};
vector<TF1*> f_splines;
vector<TF1*> f_exps;

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

void ReadAndFill1D(double dec, double ra, TH1F* prob_hist, TH1F* prob_hist2)
{
	string inpath, inpath2;

	if(ra == -1000)
	{
		inpath  = data_folder + "fdata_tStat_dec_" + to_string(dec) + "_all.txt";
		inpath2 = data_folder + "fdata_nSign_dec_" + to_string(dec) + "_all.txt";
	}
	else
	{
		inpath  = data_folder + "fdata_tStat_dec_" + to_string(dec) + "_" + to_string(ra) + "_all.txt";
		inpath2 = data_folder + "fdata_nSign_dec_" + to_string(dec) + "_" + to_string(ra) + "_all.txt";		
	}

	ifstream  inf{inpath};
	ifstream inf2{inpath2};

	if(!inf||!inf2) cout << "File not found.\n";
	int line_err = 0;

	while(inf)
	{
		string input;
		getline(inf,input);
		line_err++;

		try 
		{
			if(input != "")
			{
				prob_hist->Fill(stod(input));

				// if(stod(input) > 100)
				// {
				// 	cout << "VERY LARGE TSTAT INPUT: " << input << endl;
				// 	cout << "File: " << inpath << " Line: " << line_err << endl;
				// }
			}
    	} 
   
		catch (const std::invalid_argument&) 
		{
	        std::cerr << "Argument is invalid\n";
	        cerr << input << endl;
	        cerr << "File: " << inpath << " Line: " << line_err << endl;
	        throw;
    	}
	}

	line_err = 0;
	while(inf2)
	{
		string input;
		getline(inf2,input);
		line_err++;

		try 
		{
			if(input != "")
			{
				prob_hist2->Fill(stod(input));
				// if(stod(input) > 100)
				// {
				// 	cout << "VERY LARGE NSIGN INPUT: " << input << endl;
				// 	cout << "File: " << inpath2 << " Line: " << line_err << endl;
				// }
			}
    	} 
   
		catch (const std::invalid_argument&) 
		{
	        std::cerr << "Argument is invalid\n";
	        cerr << input << endl;
	        cerr << "File: " << inpath2 << " Line: " << line_err << endl;
	        throw;
    	}
	}
}

int get_prob()
{
	TFile* outputFile = new TFile("prob.root","RECREATE");
	THStack* hs  = new THStack("hs", "Test statistic distribution");
	THStack* hs2 = new THStack("hs2","nSignal distribution");

	for(double sigDec = -90; sigDec <= 90; sigDec += 5)
	{
		cout << "Current declination: " << sigDec << endl;

		string name1 = "prob_hist_"  + to_string(sigDec);
		string name2 = "prob_hist2_" + to_string(sigDec);

		TH1F* prob_hist  = new TH1F(name1.c_str(),"Test statistic distribution;Test statistic;Probability TS is bigger",10000,-1,100);
		TH1F* prob_hist2 = new TH1F(name2.c_str(),"nSignal distribution;nSignal;Probability nSign is bigger",10000,-1,100);
		ReadAndFill1D(sigDec,-1000,prob_hist,prob_hist2);

		TH1* prob_hist_cumul = prob_hist->GetCumulative(kFALSE);
		prob_hist_cumul->Scale(1./prob_hist->GetEntries(),"nosw2");
		TH1* prob_hist2_cumul = prob_hist2->GetCumulative(kFALSE);
		prob_hist2_cumul->Scale(1./prob_hist2->GetEntries());

		//FIT PART
		int i = (sigDec+90)/5;
		gStyle->SetOptStat(111111);

		int min_bin  = prob_hist_cumul->GetMinimumBin();
		double min_x = -1+0.0101*min_bin; //1e-8 probability when using larger set

		const int nodes = 50;
		const double low = -1;
		const double high = min_x;
		double spline_nodes[nodes];

		for (int i = 0; i < nodes; ++i)
		{
			spline_nodes[i] = low+(high-low)*i/nodes;
		}

		NSpline<nodes>* s = new NSpline<nodes>(spline_nodes);
		s->SetDer1(0);
		s->SetDer2(0);
		string fname1 = "f_probfit_"  + to_string(i);
		string fname2 = "f_probfit2_" + to_string(i);
		f_splines.push_back(new TF1(fname1.c_str(), s->GetEval(), low, high, 2*nodes+2));
		f_exps.push_back(new TF1(fname2.c_str(), "[0]*exp(-[1]*x)",-1,100));

		double tangent = (log(1e-6)-log(1e-5))/(min_x-fit_bounds[i]);
		double par0 = exp(log(1e-5)-tangent*fit_bounds[i]);
		double par1 = -tangent;

		//cout << "Parameter 0: " << par0 << ", parameter 1: " << par1 << endl;
		f_exps[i]->SetParameter(0,par0);
		f_exps[i]->SetParameter(1,par1);

		prob_hist_cumul->Fit(f_exps[i],"M","",fit_bounds[i],100);

		string cname1 = "c_spline_"  + to_string(i);
		string cname2 = "c_exp_" + to_string(i);
		TCanvas* c = new TCanvas(cname2.c_str(),cname2.c_str());
		gPad->SetLogy();
		prob_hist_cumul->Draw();
		c->Write();

		prob_hist_cumul->Fit(f_splines[i],"M","",low,high);

		TCanvas* c2 = new TCanvas(cname1.c_str(),cname1.c_str());
		gPad->SetLogy();
		prob_hist_cumul->Draw();
		c2->Write();

		hs->Add(prob_hist_cumul);
		hs2->Add(prob_hist2_cumul);
	}

	TCanvas* c1 = new TCanvas("a1","TS distribution");
	gStyle->SetOptStat(111111);
	gPad->SetLogy();
	hs->Draw("nostack plc");

	TCanvas* c2 = new TCanvas("a2","nS distribution");
	gPad->SetLogy();
	hs2->Draw("nostack plc");

	hs->Write();
	hs2->Write();
	outputFile->WriteObject(&f_exps,"f_exps");
	outputFile->WriteObject(&f_splines,"f_splines");
	outputFile->Close();

	// TFile *fin = TFile::Open("prob.root", "READ");
 //   	std::vector<TF1*> *yy;
 //   	fin->GetObject("f_exps", yy);
 //   	(*yy)[0]->Draw();

	return 0;
}