#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <Eigen/Dense>

#include "util.h"

typedef Eigen::Matrix<double, 13, 1> StateVector;

class ABMIntegrator {
	public:
		ABMIntegrator(StateVector (*func)(StateVector, double, IntegratorParams*), StateVector x0, double step, IntegratorParams* param);
		StateVector step();

		StateVector x;
		double t;

	private:
		static const int SIZE = 5;
		StateVector history[5];
		int index;
		double dt;
		StateVector (*function)(StateVector, double, IntegratorParams*);
		IntegratorParams* params;
};

#endif
