#pragma once
#include <iostream>
#include <cmath>

struct MCEvent
{
	double energy;
	double ra, dec;
	double mjd;
	char type; //B for background, S for signal
};

//operator overloading for printing
std::ostream& operator<<(std::ostream& stream, const MCEvent& b)
{
	stream << "MJD: " << b.mjd << '\n';
	stream << "Energy: " << b.energy << '\n';
	stream << "Ra: " << b.ra << " Dec: " << b.dec << '\n';

	return stream;
}

//return definite integral of power law
double powerLawArea(double power, double xmin, double xmax)
{
	return std::pow(xmax, power+1)/(power+1) - std::pow(xmin, power+1)/(power+1);
}

//comparison for sorting
bool isSooner(const MCEvent& a, const MCEvent& b)
{
	return a.mjd < b.mjd;
}