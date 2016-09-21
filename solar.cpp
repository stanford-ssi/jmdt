#include <iostream>
#include <fstream>
#include <cmath>

#include "solar.h"

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

	int row = phi/(2*M_PI)*(n-1);
	int col = theta/M_PI*(n-1);

	return data[row*n+col];
}
