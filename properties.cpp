#include <iostream>
#include <fstream>
#include <cmath>

#include <Eigen/Dense>

#include "properties.h"

using namespace Eigen;

SatelliteProperties::SatelliteProperties(string solar_file, string drag_file) {
	if (solar_file != "none") {
		ifstream file (solar_file.c_str(),
			ios::in | ios::binary | ios::ate);
		int bytes = file.tellg();

		char* mem = new char[bytes];

		file.seekg(0, ios::beg);
		file.read(mem, bytes);
		file.close();

		solar_data = (double*) mem;
		n = sqrt(bytes/sizeof(double));
		solar = true;
	} else {
		solar = false;
	}

	if (drag_file != "none") {
		cout << "wat m8" << endl;
		ifstream file (drag_file.c_str(),
			ios::in | ios::binary | ios::ate);
		int bytes = file.tellg();

		char* mem = new char[bytes];

		file.seekg(0, ios::beg);
		file.read(mem, bytes);
		file.close();

		drag_data = (double*) mem;
		drag_n = sqrt(bytes/sizeof(double));
		drag = true;
	} else {
		drag = false;
	}
}

double SatelliteProperties::get_solar_area(double theta, double phi) {
	if (!solar) return 42.0;

	// There's more resolution in theta than in phi for no good reason.

	if (phi >= 2*M_PI) {
		phi = fmod(phi, 2*M_PI);
	}
	if (theta >= M_PI) {
		theta = fmod(theta, M_PI); // Doesn't make a lot of sense.
	}
	
	int row1 = phi/(2*M_PI)*(n-1); double phi1 = 2*M_PI/(n-1)*row1;
	int col1 = theta/M_PI*(n-1); double theta1 = M_PI/(n-1)*col1;
	int row2 = row1+1; double phi2 = 2*M_PI/(n-1)*row2;
	int col2 = col1+1; double theta2 = M_PI/(n-1)*col2;

	double f11 = solar_data[row1*n+col1];
	double f12 = solar_data[row2*n+col1];
	double f21 = solar_data[row1*n+col2];
	double f22 = solar_data[row2*n+col2];
	
	double factor = 1.0/((theta2-theta1)*(phi2-phi1));
	Matrix<double, 1, 2> v1;
	v1 << (theta2-theta), (theta-theta1);
	Matrix<double, 2, 2> m1;
	m1 << f11, f12, f21, f22;
	Matrix<double, 2, 1> v2;
	v2 << (phi2-phi), (phi-phi1);

	return factor*v1*m1*v2;
}

double SatelliteProperties::get_drag_area(double theta, double phi) {
	if (!drag) return 42.0;
	// There's more resolution in theta than in phi for no good reason.

	if (phi >= 2*M_PI) {
		phi = fmod(phi, 2*M_PI);
	}
	if (theta >= M_PI) {
		theta = fmod(theta, M_PI); // Doesn't make a lot of sense.
	}
	
	int row1 = phi/(2*M_PI)*(drag_n-1);
	double phi1 = 2*M_PI/(drag_n-1)*row1;
	int col1 = theta/M_PI*(drag_n-1);
	double theta1 = M_PI/(drag_n-1)*col1;
	int row2 = row1+1;
	double phi2 = 2*M_PI/(drag_n-1)*row2;
	int col2 = col1+1;
	double theta2 = M_PI/(drag_n-1)*col2;

	double f11 = drag_data[row1*drag_n+col1];
	double f12 = drag_data[row2*drag_n+col1];
	double f21 = drag_data[row1*drag_n+col2];
	double f22 = drag_data[row2*drag_n+col2];
	
	double factor = 1.0/((theta2-theta1)*(phi2-phi1));
	Matrix<double, 1, 2> v1;
	v1 << (theta2-theta), (theta-theta1);
	Matrix<double, 2, 2> m1;
	m1 << f11, f12, f21, f22;
	Matrix<double, 2, 1> v2;
	v2 << (phi2-phi), (phi-phi1);

	return factor*v1*m1*v2;
}
