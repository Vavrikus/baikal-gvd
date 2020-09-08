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

	//creating title
	std::string title;
	std::string TimeSigma;
	if(cumulative) title = "Cumulative distribution of test statistics for background;Test statistics;Probability";
	else 		   
	{
		std::vector<std::string> info = split(files[2], '_');
		double timeSigma = stod(info[4]);
		std::string timeUnit;

		if(timeSigma < 1)
		{
			timeSigma *= 24;

			if(timeSigma < 1)
			{
				timeSigma *= 60;

				if(timeSigma < 1)
				{
					timeSigma *= 60;
					timeUnit = " s";
				}
				else timeUnit = " min";
			}
			else timeUnit = " h";
		}
		else timeUnit = " d";

		TimeSigma = std::to_string((int)std::round(timeSigma)) + timeUnit;
		title = "Distribution of test statistics with time sigma " + TimeSigma + ";Test statistics;Probability density";
	}

	TCanvas* c = new TCanvas();

	THStack* hs = new THStack("hs","");
	hs->SetTitle(title.c_str());


	for(int i = 2; i < files.size(); i++)
	{
		std::vector<double> data = readData(files[i].c_str());
		std::cout << data.size() << " numbers in " << files[i] << '\n';

		//create legend
		std::string legend;
		std::vector<std::string> info = split(files[i], '_');
		
		if(cumulative)
		{
			double timeSigma = stod(info[4]);
			std::string timeUnit;

			if(timeSigma < 1)
			{
				timeSigma *= 24;

				if(timeSigma < 1)
				{
					timeSigma *= 60;

					if(timeSigma < 1)
					{
						timeSigma *= 60;
						timeUnit = " s";
					}
					else timeUnit = " min";
				}
				else timeUnit = " h";
			}
			else timeUnit = " d";

			legend = "time sigma " + std::to_string((int)std::round(timeSigma)) + timeUnit;
		}

		else
		{
			int numSignal = std::stoi(info[2]);

			if(numSignal == 1) legend = std::to_string(numSignal) + " signal event";
			else			   legend = std::to_string(numSignal) + " signal events";
		}
		
		double bins;

		if(cumulative) bins = 1e+5;
		else bins = std::floor(std::pow(data.size(),0.5));

		XY extrems = minmax(data);
	
		TH1D* hist = new TH1D(std::to_string(i).c_str(),legend.c_str(),bins, extrems.x, extrems.y);
	
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
		hs->Add(hist2, "hist");		
	}

	hs->Draw("nostack");
	double yLower = 0.89 - (files.size()-2)*0.04;
	gPad->BuildLegend(0.7,yLower,0.89,0.89,"");


	if(drawSigmas)
	{
		if(cumulative)
		{
			hs->GetXaxis()->SetRangeUser(0,25);
			double xMax = hs->GetXaxis()->GetXmax();
			
			/*/TLine* l1 = new TLine(0,0.31731050786291415,xMax,0.31731050786291415);
			l1->Draw();
			TLine* l2 = new TLine(0,0.04550026389635842,xMax,0.04550026389635842);
			l2->Draw();
			TLine* l4 = new TLine(0,0.00006334248366624,xMax,0.00006334248366624);
			l4->Draw();*/

			double sigma3 = 0.002699796063260207;
			double sigma5 = 5.733031436805369e-7;

			double xOffset = 0.92;
			double yOffset = 0.3;

			TLine* l3 = new TLine(0,sigma3,xMax,sigma3);
			l3->Draw();
			
			TLatex* label3 = new TLatex(xMax*xOffset, sigma3*yOffset, "3\\sigma");
			label3->Draw();
			
			TLine* l5 = new TLine(0,sigma5,xMax,sigma5);
			l5->Draw();

			TLatex* label5 = new TLatex(xMax*xOffset, sigma5*yOffset,"5\\sigma");
			label5->Draw();
		}

		else
		{
			hs->SetMaximum(0.01);
			hs->GetXaxis()->SetRangeUser(0,80);
			c->Update();
			c->GetFrame()->GetY2();
			double yMax   = gPad->GetUymax();

			double sigma3 = 0;
			double sigma5 = 0;

			switch(timeSigma)
			{
				case "1 d":
					sigma3 = 3.37;
					sigma5 = 15.12;
					break;
				case "15 min":
					sigma3 = 0;
					sigma5 = 0;
					break;
				case "10 s":
					sigma3 = 0;
					sigma5 = 0;
					break;
			}

			double xOffset = 0.5;
			double yOffset = 0.9;

			TLine* l3 = new TLine(sigma3,0,sigma3,yMax);
			l3->Draw();

			l3->SetLineWidth(2);
			
			TLatex* label3 = new TLatex(sigma3+xOffset, yMax*yOffset, "3\\sigma");
			label3->Draw();
			
			TLine* l5 = new TLine(sigma5,0,sigma5,yMax);
			l5->Draw();
			
			l5->SetLineWidth(2);

			TLatex* label5 = new TLatex(sigma5+xOffset, yMax*yOffset, "5\\sigma");
			label5->Draw();			
		}
	}

	if(cumulative) gPad->SetLogy();
}