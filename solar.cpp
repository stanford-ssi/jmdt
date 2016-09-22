#include <iostream>
#include <fstream>
#include <cmath>

#include <Eigen/Dense>

#include "solar.h"

using namespace Eigen;

SolarData::SolarData(string filename) {
	ifstream file (filename.c_str(), ios::in | ios::binary | ios::ate);
	int bytes = file.tellg();

	char* mem = new char[bytes];

	file.seekg(0, ios::beg);
	file.read(mem, bytes);
	file.close();

	data = (double*) mem;
	n = sqrt(bytes/sizeof(double));
}

double SolarData::get_area(double theta, double phi) {
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

	double f11 = data[row1*n+col1];
	double f12 = data[row2*n+col1];
	double f21 = data[row1*n+col2];
	double f22 = data[row2*n+col2];
	
	double factor = 1.0/((theta2-theta1)*(phi2-phi1));
	Matrix<double, 1, 2> v1;
	v1 << (theta2-theta), (theta-theta1);
	Matrix<double, 2, 2> m1;
	m1 << f11, f12, f21, f22;
	Matrix<double, 2, 1> v2;
	v2 << (phi2-phi), (phi-phi1);

	return factor*v1*m1*v2;
}
