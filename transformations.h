#pragma once

#include "TVector3.h"
#include "TLatex.h"
#include "TH2F.h"
#include "TGraph.h"

extern constexpr double PI = 3.14159265358979323;

//coordinates of Baikal-GVD converted to radians (assuming east longitude)
extern constexpr double latB = PI * 51.768 / 180.0;  //about 3.6 km south from 107 km station
extern constexpr double lonB = PI * 104.399 / 180.0; //about 3.6 km south from 107 km station

//IceCube coordinates
extern constexpr double latI = PI * 89.99 / 180.0;
extern constexpr double lonI = PI * -63.453056 / 180.0;

//ANTARES coordinates
extern constexpr double latA = PI * 42.8 / 180.0;
extern constexpr double lonA = PI * 6.167 / 180.0;

inline double radToDeg(const double& r) {return 180*r/PI;} //convert radians to degrees
inline double degToRad(const double& d) {return PI*d/180;} //convert degrees to radians

//spherical cap angle in steradians
inline double capAngle(const double& sideAngle) {return 2*PI*(1-std::cos(sideAngle));}

//distance between two points with given equatorial coordinates (degrees)
double angularDistance(const double& ra, const double& dec, const double& ra2, const double& dec2)
{
    TVector3 v1(0,0,1);
    v1.SetTheta(degToRad(90+dec));
    v1.SetPhi(degToRad(ra));

    TVector3 v2(0,0,1);
    v2.SetTheta(degToRad(90+dec2));
    v2.SetPhi(degToRad(ra2));

    return radToDeg(v1.Angle(v2));
}

//function for splitting string by whitespace into substrings 
std::vector<std::string> split(const std::string& str, char separator = ' ')
{
    std::vector<std::string> output;
    int currentPos = 0;
    int lastSpacePos = -1;

    for (char c : str)
    {
        if (c == separator) 
        {
            std::string substr = str.substr(lastSpacePos + 1, currentPos - lastSpacePos - 1);
            lastSpacePos = currentPos;
            if (substr != "") output.push_back(substr);
        }

        currentPos++;
    }

    std::string substr = str.substr(lastSpacePos + 1, currentPos - lastSpacePos - 1);
    if (substr != "") output.push_back(substr);

    return output;
}

//========================== TIME FORMAT CONVERSION ==========================//

//convert UTC time to Unix
double UTCtoUnix(const int& year, const int& month, const int& day, const int& hours,
                  const int& minutes, const double& seconds)
{
    int leapYears = std::floor((year - 1969) / 4);
    leapYears += std::floor((year - 2000) / 100);
    leapYears -= std::floor((year - 2000) / 400);

    int monthDays = std::ceil((month - 1) * 30.5);

    if (month > 2)
    {
        monthDays -= 2;
        if (year % 4 == 0 and year % 400 != 0) monthDays++;
        if (month == 9 or month == 11) monthDays++;
    }

    return ((year - 1970) * 365 + leapYears + monthDays + day - 1) * 86400 + hours * 3600 
            + minutes * 60 + seconds;
}

//convert unix time to modified julian date
inline double UnixToMJD(const double& unix) {return unix/86400 + 40587;}

//convert modified julian date to unix time
inline double MJDtoUnix(const double& mjd) {return (mjd - 40587) * 86400;}

//convert UTC time to modified julian date 
inline double UTCtoMJD(const int& year, const int& month, const int& day, const int& hours,
                 		const int& minutes, const double& seconds)
{
    return UnixToMJD(UTCtoUnix(year, month, day, hours, minutes, seconds));
}

//==================== ASTRONOMICAL COORDINATES TRANSFORM ====================//

//struct for storing horizontal coordinates of event
//altitude and azimuth in degrees, time in Unix format
struct horCoor
{
	double alt, azm, uTime;

	horCoor() {}
	horCoor(double alt, double azm, double uTime)
		: alt(alt), azm(azm), uTime(uTime) {}
};

//struct for storing equatorial coordinates of event
//right ascension in hours and declination in degrees
struct eqCoor
{
	double rAsc, dec;

	eqCoor() {}
	eqCoor(double rAsc, double dec)
		: rAsc(rAsc), dec(dec) {}
};

//calculating local sidereal time in degrees acording to https://www.aa.quae.nl/en/reken/sterrentijd.html
double LST(double unixTime, double longitude)
{
    //some constants from the website
    constexpr double l0 = 99.967794687;
    constexpr double l1 = 360.98564736628603;
    constexpr double l2 = 2.907879e-13;
    constexpr double l3 = -5.302e-22;

    double daysSince2000 = (unixTime - 946684800.0) / 86400.0;

    //faster but can be less precise than fmod on each
    double localSiderealTime = l0 + l1 * daysSince2000 + l2 * std::pow(daysSince2000, 2)
        + l3 * std::pow(daysSince2000, 3) + radToDeg(longitude);

    //converting to hours
    localSiderealTime = std::fmod(localSiderealTime, 360.0);
    if (localSiderealTime < 0) localSiderealTime += 360.0;

    return localSiderealTime;
}

//trasform from horizontal to equatorial coordinates
eqCoor horToEq(const horCoor& hor, const double& latDet = latB, const double& lonDet = lonB)
{
    using namespace std;

    //converting azimuth and altitude to radians
    double alt = degToRad(hor.alt);
    double azm = degToRad(hor.azm);

    //transform using equations from http://star-www.st-and.ac.uk/~fv/webnotes/chapter7.htm
    double sinDec = sin(alt) * sin(latDet) + cos(alt) * cos(azm) * cos(latDet);

    //making sure asin does not return NaN
    if(sinDec > 1) sinDec = 1;
    else if (sinDec < -1) sinDec = -1;

    double declination = asin(sinDec);

    double cosLHA = (sin(alt) - sin(declination) * sin(latDet)) / (cos(declination) * cos(latDet));

    //making sure acos does not return NaN
    if(cosLHA > 1) cosLHA = 1;
    else if (cosLHA < -1) cosLHA = -1;

    double localHourAngle = acos(cosLHA);

    //expressing localHourAngle in degrees
    localHourAngle = radToDeg(localHourAngle);

    if (sin(azm) > 0) localHourAngle = 360.0 - localHourAngle;

    double localSiderealTime = LST(hor.uTime, lonDet);

    //calculating right ascension
    double rightAscension = localSiderealTime - localHourAngle;
    if (rightAscension < 0) rightAscension += 360.0;

    return eqCoor{ rightAscension, radToDeg(declination) };
}

//trasform from equatorial to horizontal coordinates
horCoor eqToHor(const eqCoor& eq, const double& uTime, const double& latDet = latB, const double& lonDet = lonB)
{
    using namespace std;

    //converting right ascension and declination to radians
    double rAsc = degToRad(eq.rAsc);
    double dec = degToRad(eq.dec);

    double localSiderealTime = LST(uTime, lonDet);
    double LHA = degToRad(localSiderealTime) - rAsc; //local hour angle

    if (LHA < 0) LHA += 2 * PI;

    //transform using equations from http://star-www.st-and.ac.uk/~fv/webnotes/chapter7.htm
    double sinAlt = sin(dec) * sin(latDet) + cos(dec) * cos(latDet) * cos(LHA);

    //making sure asin does not return NaN
    if(sinAlt > 1) sinAlt = 1;
    else if (sinAlt < -1) sinAlt = -1;

    double altitude = asin(sinAlt);

    double cosAz = (sin(dec) - sin(latDet) * sin(altitude)) / (cos(latDet) * cos(altitude));

    //making sure acos does not return NaN
    if(cosAz > 1) cosAz = 1;
    else if (cosAz < -1) cosAz = -1;

    double azimuth  = acos(cosAz);

    if (sin(LHA) > 0) azimuth = 2 * PI - azimuth;

    return horCoor{ radToDeg(altitude), radToDeg(azimuth), uTime };
}

//transform from shifted spherical coordinates given coordinates of north pole (npoleDec = 90 fails)
eqCoor shiftSpherTrans(const double& zenithAngle, const double& Azm, const double& npoleRa,
					   const double& npoleDec)
{
    using namespace std;

    //converting azimuth and altitude to radians
    double alt = degToRad(90-zenithAngle);
    double azm = degToRad(Azm);

    //transform using equations from http://star-www.st-and.ac.uk/~fv/webnotes/chapter7.htm
    double sinDec = sin(alt) * sin(npoleDec) + cos(alt) * cos(azm) * cos(npoleDec);

    //making sure asin does not return NaN
    if(sinDec > 1) sinDec = 1;
    else if (sinDec < -1) sinDec = -1;

    double declination = asin(sinDec);

    double cosLHA = (sin(alt) - sin(declination) * sin(npoleDec)) / (cos(declination) * cos(npoleDec));

    //making sure acos does not return NaN
    if(cosLHA > 1) cosLHA = 1;
    else if (cosLHA < -1) cosLHA = -1;

    double localHourAngle = acos(cosLHA);

    //expressing localHourAngle in degrees
    localHourAngle = radToDeg(localHourAngle);

    if (sin(azm) > 0) localHourAngle = 360.0 - localHourAngle;

    //calculating right ascension
    double rightAscension = radToDeg(npoleRa) - localHourAngle;
    if (rightAscension < 0) rightAscension += 360.0;

    return eqCoor{ rightAscension, radToDeg(declination) };
}

//============================== DRAWING AITOFF ==============================//

struct XY { double x,y; };

//aitoff transform from original drawmap.C file converted to function
XY toAitoff(double x, double y)
{
	double z = sqrt(1+cos(degToRad(y))*cos(degToRad(x/2)));
	double x2 = 180*cos(degToRad(y))*sin(degToRad(x/2))/z;
	double y2 = 90*sin(degToRad(y))/z;

	return XY{x2,y2};
}

//draws labels on aitoff map
void drawLabels()
{
	//declination labels
	int stepDec = 15;

	for(int i = -90; i <= 90; i += stepDec)
	{
		double xOffset = -14.5;
		double yOffset = -5+i/12;

		if(i < 0) xOffset -= 3;

		XY coor = toAitoff(-180,i);
		TLatex* label = new TLatex(coor.x + xOffset, coor.y + yOffset, 
							("#scale[0.75]{"+std::to_string(i)+"\260}").c_str());
		label->Draw();
	}

	//right ascension labels
	int stepRa = 45;

	for(int i = -180; i <= 180; i += stepRa)
	{
		double xOffset = 2;
		double yOffset = -9;

		XY coor = toAitoff(i,0);
		int hours = i/15;
		if(hours < 0) hours += 24;

		TLatex* label = new TLatex(coor.x + xOffset, coor.y + yOffset, 
							("#scale[0.7]{"+std::to_string(hours)+" h}").c_str());
		label->Draw();
	}
}

//https://root-forum.cern.ch/t/how-to-plot-skymaps-with-root/4663/3 with minor changes
TH2F* drawmap(const char* title)
{
	TH2F* skymap = new TH2F("Blank skymap","Blank skymap",40,-210,210,20,-105,105);

	//gStyle->SetOptStat(0000);
		
	const int Nl=13; // Number of drawn latitudes
	const int NL=25; // Number of drawn longitudes
	int M=30;

	TGraph*	latitudes[Nl];
	TGraph*	longitudes[NL];

	for(int j=0;j<Nl;++j)
	{
		latitudes[j]=new TGraph();
		double la=-90+180/(Nl-1)*j;
		
		for(int i=0;i<M+1;++i)
		{
			double lo= -180+360/M*i;
			XY aitoff = toAitoff(lo,la);
			latitudes[j]->SetPoint(i,aitoff.x,aitoff.y);
		}
	}  

	for(int j=0;j<NL;++j)
	{
		longitudes[j]=new TGraph();
		double lo=-180+360/(NL-1)*j;
		
		for(int i=0;i<M+1;++i)
		{
			double la= -90+180/M*i;
			XY aitoff = toAitoff(lo,la);
			longitudes[j]->SetPoint(i,aitoff.x,aitoff.y);
		}
	}

	skymap->Draw("AH");
	skymap->SetTitle(title);

	for(int j=0;j<Nl;++j) latitudes[j]->Draw("c");
	for(int j=0;j<NL;++j) longitudes[j]->Draw("c");

	return skymap;
}