#include "transformations.h"

constexpr int samples = 200;
constexpr double eventTime = 1502973664; // 17.8.2017 12:41:04 UTC

//templated struct for two arrays of same size used in sortSplitted
template<int N>
struct XYarr
{
	double x[N],y[N];
	int size() {return N;}

	//manually copying arrays in constructor
    XYarr(const double (&x_)[N], const double (&y_)[N]) 
    {
        memcpy(&x,x_,sizeof (double)*N);
        memcpy(&y,y_,sizeof (double)*N);
    }
};

//checks if splitted array is ascending using first 4 elements
bool isAscending(double arr[])
{
	int num = (int)(arr[0] < arr[1]) + (int)(arr[1] < arr[2]) + (int)(arr[2] < arr[3]);
	return num > 1; 
}

//sort x and y array with splitted data
template<int size>
XYarr<size> sortSplitted(double x[], double y[])
{
	int splitIndex = 0;

	if(isAscending(x))
	{
		for (int i = 0; i < size-1; ++i) 
			if(x[i] > 0 and x[i+1] < 0) splitIndex = i+1;
	}

	else
	{
		for (int i = 0; i < size-1; ++i) 
			if(x[i] < 0 and x[i+1] > 0) splitIndex = i+1;		
	}

	double x2[size], y2[size];

	for (int i = splitIndex; i < size; ++i)
	{
		x2[i - splitIndex] = x[i];
		y2[i - splitIndex] = y[i];
	}

	for (int i = 0; i < splitIndex; ++i)
	{
		x2[i + size - splitIndex] = x[i];
		y2[i + size - splitIndex] = y[i];
	}

	return XYarr<size>{x2, y2};
}

//draws horizon of certain location in certain time
TGraph* drawHorizon(const double time, const double latDet = latB, 
				 const double lonDet = lonB, EColor color = kRed)
{
	
	double ra[samples], dec[samples], azm[samples];

	for (int i = 0; i < samples; ++i) azm[i] = ((double)i/(double)samples)*360;
	for (int i = 0; i < samples; ++i)
	{
		eqCoor horizon = horToEq(horCoor(0,azm[samples-i-1],time),latDet,lonDet);
		ra[i]  = horizon.rAsc;
		dec[i] = horizon.dec;
	}


	for (int i = 0; i < samples; ++i) 
	{
		//ra[i] = 180-ra[i];
		if(ra[i]>180) 
		{
			ra[i] -= 360;
		}
	}

	auto h = sortSplitted<samples>(ra, dec);

	double x[samples], y[samples];

	for (int i = 0; i < samples; ++i)
	{
		XY data = toAitoff(h.x[i], h.y[i]);
		x[i] = data.x;
		y[i] = data.y;
	}

	TGraph* horizon = new TGraph(samples, x, y);

	horizon->SetLineColor(color);
	horizon->SetLineWidth(2);
	horizon->Draw("L");

	return horizon;

}

void horizon()
{
	TCanvas* c1 = new TCanvas("c1", "", 1000, 500);
	drawmap("Neutrino detectors horizons GW170817;Right ascension;Declination");

	TGraph* bai = drawHorizon(eventTime);
	TGraph* ice = drawHorizon(eventTime, latI, lonI, kBlue);
	TGraph* ant = drawHorizon(eventTime, latA, lonA, kGreen);

	TLegend* leg = new TLegend(.76,.725,.9,.9,"");

	leg->AddEntry(bai, "Baikal-GVD");
	leg->AddEntry(ice, "IceCube");
	leg->AddEntry(ant, "ANTARES");
	leg->Draw("Same");

	drawLabels();
}


//basic projection

	/*TGraph* horizon = new TGraph(samples, ra2, dec2);
	TAxis* axis = horizon->GetXaxis();

	axis->SetLimits(0.,360.);
	horizon->GetHistogram()->SetMinimum(-90.);
	horizon->GetHistogram()->SetMaximum(90.);

	horizon->SetTitle("Baikal-GVD horizon;Right ascension [degrees];Declination [degrees]");
	horizon->SetLineColor(kRed);
	horizon->Draw("AL");*/