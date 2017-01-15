#ifndef UTIL_H
#define UTIL_H

#include <Eigen/Dense>
#include <GeographicLib/GravityModel.hpp>
#include <GeographicLib/Geocentric.hpp>

#include "properties.h"

const double EARTH_MU = 3.986004418e14;
const double EARTH_RADIUS = 6371009.0;
const Eigen::Vector3d EARTH_OMEGA(0, 0, 7.2921150e-5);
const double AU = 149597870700;
const double SOLAR_CONSTANT = 1360.8;

const double DEG2RAD = M_PI/180.0;

typedef Eigen::Matrix<double, 13, 1> StateVector;

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
	double solar_efficiency;
	string orientation_str;
	int two_satellites;
	StateVector* lover;

	int i_am_leader; // Unused
	int separation_target;
	int controller_behavior;

	double f107A;
	double f107;
	double ap;

	double filter; // Unused
	double doy;
	double sec;

	double output_power;
	double output_BC;
};

Eigen::Vector3d earth_sun_vector(double jd);

bool satellite_in_shade(Eigen::Vector3d& r_sat, Eigen::Vector3d r_sun);

void jd_to_date(double jd, double& year, double& month, double& day, double& hour, double& minute, double& second);

const GeographicLib::Geocentric EARTH = GeographicLib::Geocentric::WGS84();

#endif
