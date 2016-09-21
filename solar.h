#ifndef SOLAR_H
#define SOLAR_H

#include <string>

using namespace std;

class SolarData {
	public:
		SolarData(string file);
		double get_area(double theta, double phi);

	private:
		double* data;
		int n;
};

#endif
