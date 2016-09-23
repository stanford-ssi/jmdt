#ifndef UTIL_H
#define UTIL_H

#include <Eigen/Dense>
#include <GeographicLib/GravityModel.hpp>

#include "properties.h"

const double EARTH_MU = 3.986004418e14;
const double EARTH_RADIUS = 6371009.0;
const Eigen::Vector3d EARTH_OMEGA(0, 0, 7.2921150e-5);
const double AU = 149597870700;

const double DEG2RAD = M_PI/180.0;

struct IntegratorParams {
	int earth;
	int atmosphere;
	GeographicLib::GravityModel* gravity_object;
	double Cd;
	double A;
	double mass;
	int power;
	SatelliteProperties* properties;
	double t0;
	Eigen::Vector3d orientation;
	int drag;
	double output_power;
	double output_BC;
};

Eigen::Vector3d earth_sun_vector(double jd);

#endif
