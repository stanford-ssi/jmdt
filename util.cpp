#include <cmath>

#include <Eigen/Dense>
#include <Eigen/Geometry>

#include "util.h"

using namespace Eigen;

/* See: David A. Vallado. Fundamentals of Astrodynamics and Applications. */
Vector3d earth_sun_vector(double jd) {
	double ut1 = (jd-2451545.0)/36525.0;
	double lambdam0 = 280.4606184+36000.77005361*ut1;
	double m0 = 357.5277233+35999.05034*ut1;
	double lambdaecl = lambdam0 + 1.914666471*sin(m0*DEG2RAD) +
					0.019994643*sin(2*m0*DEG2RAD);
	double r0 = 1.000140612 - 0.016708617*cos(m0*DEG2RAD) -
					0.000139589*cos(2*m0*DEG2RAD);
	double epsilon = 23.439291-0.0130042*ut1;

	lambdaecl = lambdaecl*DEG2RAD;
	epsilon = epsilon*DEG2RAD;
	r0 = r0*AU;
		
	Vector3d out(r0*cos(lambdaecl), r0*cos(epsilon)*sin(lambdaecl),
		r0*sin(epsilon)*sin(lambdaecl));

	return out;
}

/* r_sat is the geocentric vector of the satellite, r_sun is a Sun-Earth
 * vector. This function just assumes that the shadow is a cylinder, and
 * tests whether the satellite is inside said cylinder. */
bool satellite_in_shade(Vector3d& r_sat, Vector3d r_sun) {

	r_sun = r_sun.normalized();

	double dot = r_sat.dot(r_sun);
	if (dot < 0.0) return false;
	else {
		double dst2 = r_sat.squaredNorm() - dot*dot;
		if (dst2 < EARTH_RADIUS*EARTH_RADIUS) return true;
		else return false;
	}
}

void jd_to_date(double jd, double& year, double& month, double& day, double& hour, double& minute, double& second) {
	jd += 0.5;
	double z = floor(jd);
	double f = jd-z;

	double a;
	if (z < 2299161) {
		a = z;
	} else {
		double alpha = floor((z-1867216.25)/36524.25);
		a = z + 1 + alpha - floor(alpha/4);
	}
	double b = a+1524;
	double c = floor( (b-122.1)/365.25);
	double d = floor(365.25*c);
	double e = floor((b-d)/30.6001);
	double _dd = b - d - floor(30.6001*e) + f;
	day = floor(_dd);
	double _hh = (_dd-day)*24;
	hour = floor(_hh);
	double _minu = (_hh-hour)*60;
	minute = floor(_minu);
	month = (e < 13.5) ? (e-1) : (e-13);
	year = (month > 2.5) ? (c-4716) : (c-4715);
}
