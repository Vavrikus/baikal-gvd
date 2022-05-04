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

string bounds_folder = "/media/vavrik/Elements/Baikal-GVD/lcoincidences/fit_results/test_03_dec_5degstep_full/"; //"./data/merged/";//
string data_folder = "/media/vavrik/Elements/Baikal-GVD/lcoincidences/fit_results/test_06_dec_5degstep_signal_0_to_20/"; //"./data/merged/";//

//1e-5 probability bounds
//const vector<double> fit_bounds = {10.9983,10.888,10.6235,10.4007,10.3206,10.3551,10.3797,10.5331,10.6906,10.8716,11.1583,11.6488,12.2753,13.1082,14.1282,15.2548,16.4592,17.9392,19.7688,22.1002,24.8501,28.3045,31.8804,35.5551,39.5829,44.0158,49.259,54.5009,58.8381,61.8463,63.9316,65.1841,66.6971,68.6378,70.9946,72.7542,73.3589};

//0.05 probability bounds
//const vector<double> ts_p005    = {0.287219,0.34949,0.512512,0.675874,0.80989,0.964257,1.10151,1.15451,1.22432,1.31015,1.41855,1.56639,1.76518,2.03345,2.38835,2.84556,3.41799,4.12712,5.04739,6.27374,7.86129,9.77794,11.9026,14.1414,16.5465,19.3086,22.4561,25.6454,28.4524,30.6246,32.1805,33.3783,34.629,36.2366,38.0266,39.4606,40.0134};

//TS bounds for given p-value
vector<double> ts_bounds;

//data cosine of theta distribution
static TF1* costheta;

//signal distribution degrees
static double sigRa 	   = 0;
static double sigDec	   = 0;
static double sigPosSigma  = 5;

void GetTSBound(double dec, double pval)
{
	string inpath, inpath2;
	double prob_bound = 1-pval;
	std::vector<double> tStat;

	inpath  = bounds_folder + "fdata_tStat_dec_" + to_string(dec) + "_all.txt";
	inpath2 = bounds_folder + "fdata_nSign_dec_" + to_string(dec) + "_all.txt";

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
				tStat.push_back(stod(input));
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

	sort(tStat.begin(),tStat.end());
	int low_index = floor(prob_bound*tStat.size());
	int high_index = ceil(prob_bound*tStat.size());
	double shift = (prob_bound*tStat.size()-low_index)/(high_index-low_index);
	if (low_index == high_index) shift = 0;
	double TS_bound = tStat[low_index]+shift*(tStat[high_index]-tStat[low_index]);
	ts_bounds.push_back(TS_bound);
}

void ReadAndFill1D(int sig_evs, int dec_index, double dec_step, vector<vector<double>>& dir, vector<vector<double>>& exc)
{
	string inpath = data_folder + "fdata_tStat_dec_";
	inpath += to_string(-90+dec_index*dec_step) + "_SE" + to_string(sig_evs) + "_all.txt";

	ifstream  inf{inpath};

	if(!inf) cout << "File not found.\n";
	int line_err = 1;

	while(inf)
	{
		bool p_bound_reached = false;
		
		for (int i = 0; i <= sig_evs; ++i)
		{
			string input;
			getline(inf,input);

			try 
			{
				if(input != "")
				{
					double d_input = stod(input);
					if(d_input > ts_bounds[dec_index])
					{
						dir[dec_index][i]++;
						if(!p_bound_reached)
						{
							exc[dec_index][i]++;
							p_bound_reached = true;
						}
					}
				}
				else break;
	    	} 
	   
			catch (const std::invalid_argument&) 
			{
		        std::cerr << "Argument is invalid\n";
		        cerr << input << endl;
		        cerr << "File: " << inpath << " Line: " << line_err << endl;
		        throw;
	    	}

	    	if((!p_bound_reached) && (i == sig_evs))
				exc[dec_index][sig_evs+1]++;
		}

		line_err += sig_evs+1;
	}
}

int disc_pot()
{
	int sig_evs = 20;
	double dec_step = 5;
	int nSimul = 1000000;
	double pval = 0.05;

	for(double sigDec = -90; sigDec <= 90; sigDec += 5)
	{
		cout << "Current declination: " << sigDec << endl;
		GetTSBound(sigDec,pval);
	}

	//n[declination_index][signal events]
	vector<vector<double>> n_direct;
	vector<vector<double>> n_exclusive;

	//vector initialization, from now on values will be incremented
	for (int j = 0; j <= floor(180.0/dec_step); ++j)
	{
		vector<double> direct, exclusive;

		for(int i = 0; i <= sig_evs+1; ++i) {direct.push_back(0); exclusive.push_back(0);}
		n_direct.push_back(direct); n_exclusive.push_back(exclusive);
	}

	for (int j = 0; j <= floor(180.0/dec_step); ++j)
	{
		ReadAndFill1D(sig_evs,j,dec_step,n_direct,n_exclusive);
	}

	vector<TH1F*> vh_percent;
	for (int i = 0; i <= sig_evs; ++i)
	{
		string h_name  = "h_percent"+to_string(i);
		string h_title = "\% of simulations with p-value "+to_string(pval)+" or lower ("+to_string(i)+" signal events);Declination [°];\% simulations"; 

		vh_percent.push_back(new TH1F(h_name.c_str(),h_title.c_str(),floor(180/dec_step)+1,-90,90));
	}

	string h_title = "Mean value of signal events needed to get bellow p-value "+to_string(pval)+";Declination;Mean value";
	TH1F* h_mean = new TH1F("h_mean",h_title.c_str(),floor(180/dec_step)+1,-90,90);

	for(int j = 0; j <= floor(180.0/dec_step); ++j)
	{
		double mean = 0;

		for (int i = 0; i <= sig_evs; ++i)
		{
			vh_percent[i]->SetBinContent(j+1,100*n_direct[j][i]/nSimul);
			mean += i*n_exclusive[j][i]/nSimul;
		}

		cout << "Overflow for declination " << -90+j*dec_step << ": " << n_exclusive[j][sig_evs+1] << endl;

		h_mean->SetBinContent(j+1,mean);
	}

	TCanvas* c1 = new TCanvas("c1");
	h_mean->Draw();

	for(auto h : vh_percent)
	{
		TCanvas* c = new TCanvas();
		h->Draw(); 
	}

	cout << endl << "N_DIRECT:" << endl;
	for(auto v : n_direct) {for(auto d : v) {cout << d << " ";} cout<<endl;}
	cout << endl << "N_EXCLUSIVE:" << endl;
	for(auto v : n_exclusive) {for(auto d : v) {cout << d << " ";} cout<<endl;}

	TFile* outputFile = new TFile("potential.root","RECREATE");
	for(auto h : vh_percent) h->Write();
	h_mean->Write();

	return 0;
}

int main(int argc, char** argv)
{
	return disc_pot();
}