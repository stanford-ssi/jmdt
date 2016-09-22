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
