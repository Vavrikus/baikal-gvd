#include "../transformations.h"

#include <vector>
#include <fstream>
#include <string>

//returns minimal (x) and maximal (y) element of array in XY struct
XY minmax(std::vector<double> vec)
{
    double min = vec[0];
    double max = vec[0];

    for (int i = 0; i < vec.size(); ++i)
    {
        if(vec[i] > max) max = vec[i];
        if(vec[i] < min) min = vec[i];
    }

    return XY{min,max};
}

//read test statistics from txt file
std::vector<double> readData(const char* path)
{
	std::vector<double> dataOut;
	std::ifstream inf{ path };

	if(!inf) std::cerr << "File could not be opened.\n";
	
	std::string line;

	while(inf)
	{
		std::getline(inf, line);
		
		try{if(line != "") dataOut.push_back(std::stod(line));}
		catch (...){std::cout << "ERROR: " << line << "\n";}
	}

	return dataOut;
}

//input file paths separated by commas
void plotHist(std::string paths)
{
	std::vector<std::string> files = split(paths, ',');

	bool cumulative = std::stoi(files[0]);
	bool drawSigmas = std::stoi(files[1]);

	THStack* hs = new THStack("hs","");
	hs->SetTitle("Cumulative distribution of test statistics;Test statistics;Probability");


	for(int i = 2; i < files.size(); i++)
	{
		std::vector<double> data = readData(files[i].c_str());
		std::cout << data.size() << " numbers in " << files[i] << '\n';
		
		double bins;

		if(cumulative) bins = 1e+5;
		else bins = std::floor(std::pow(data.size(),0.5));

		XY extrems = minmax(data);
	
		TH1D* hist = new TH1D(std::to_string(i).c_str(),files[i].c_str(),bins, extrems.x, extrems.y);
	
		for(double d : data) hist->Fill(d);
	
		hist->Scale(1/hist->Integral());
		
		TH1D* hist2;

		if(cumulative)
		{
			hist2 = (TH1D*)hist->GetCumulative(kFALSE);
			delete hist;
		}
		else hist2 = hist;
	
		hist2->SetLineColor(i);
		hs->Add(hist2, "l");		
	}

	hs->Draw("nostack l");
	gPad->BuildLegend(0.75,0.75,0.95,0.95,"");

	double xMax = hs->GetXaxis()->GetXmax();

	if(drawSigmas)
	{
		TLine* l1 = new TLine(0,0.31731050786291415,xMax,0.31731050786291415);
		l1->Draw();
		TLine* l2 = new TLine(0,0.04550026389635842,xMax,0.04550026389635842);
		l2->Draw();
		TLine* l3 = new TLine(0,0.002699796063260207,xMax,0.002699796063260207);
		l3->Draw();
		TLine* l4 = new TLine(0,0.00006334248366624,xMax,0.00006334248366624);
		l4->Draw();
		TLine* l5 = new TLine(0,5.733031436805369e-7,xMax,5.733031436805369e-7);
		l5->Draw();
	}

	if(cumulative) gPad->SetLogy();
}