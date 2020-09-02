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
		
		if(line != "") dataOut.push_back(std::stod(line));
	}

	return dataOut;
}

void plotHist(const char* path)
{
	std::vector<double> data = readData(path);
	std::cout << data.size() << '\n';

	double bins = 400;
	XY extrems = minmax(data);

	double binWidth = (extrems.y-extrems.x)/bins;

	TH1F* hist = new TH1F("Histogram","Cumulative distribution of test statistics for background",
						  bins, extrems.x, extrems.y);

	for(double d : data) hist->Fill(d);

	hist->Scale(1/hist->Integral());
	//hist = (TH1F*)hist->GetCumulative(kFALSE);

	hist->Draw("hist");
	hist->SetTitle("Cumulative distribution of test statistics for background;Test statistics;Probability");
	hist->SetLineColor(kRed);

	/*
	TLine* l1 = new TLine(0,0.31731050786291415,5,0.31731050786291415);
	l1->Draw();
	TLine* l2 = new TLine(0,0.04550026389635842,7,0.04550026389635842);
	l2->Draw();
	TLine* l3 = new TLine(0,0.002699796063260207,7,0.002699796063260207);
	l3->Draw();
	TLine* l4 = new TLine(0,0.00006334248366624,7,0.00006334248366624);
	l4->Draw();
	TLine* l5 = new TLine(0,5.733031436805369e-7,7,5.733031436805369e-7);
	l5->Draw();
	*/

	//gPad->SetLogy();
}