#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <Eigen/Dense>
#include <GeographicLib/GravityModel.hpp>

const double EARTH_MU = 3.986004418e14;
const double EARTH_RADIUS = 6371009.0;
const Eigen::Vector3d EARTH_OMEGA(0, 0, 7.2921150e-5);

struct IntegratorParams {
	int earth;
	int atmosphere;
	GeographicLib::GravityModel* gravity_object;
	double Cd;
	double A;
	double mass;
};

#endif
