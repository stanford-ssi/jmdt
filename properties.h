#ifndef PROPERTIES_H
#define PROPERTIES_H

#include <string>

using namespace std;

class SatelliteProperties {
	public:
		SatelliteProperties(string solar, string drag);
		double get_solar_area(double theta, double phi);
		double get_drag_area(double theta, double phi);

	private:
		double* solar_data;
		double* drag_data;
		int n;
		int drag_n;
		bool solar;
		bool drag;
};

#endif
