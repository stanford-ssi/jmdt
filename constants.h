#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <GeographicLib/GravityModel.hpp>

const double EARTH_MU = 3.986004418e14;

struct IntegratorParams {
	int earth;
	int atmosphere;
	GeographicLib::GravityModel* gravity_object;
};

#endif
