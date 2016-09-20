#include "integrator.h"
#include <iostream>

using namespace std;

ABMIntegrator::ABMIntegrator(StateVector (*func)(StateVector, double, IntegratorParams*), StateVector x0, double step, IntegratorParams* param) {
	dt = step;
	t = 0;
	function = func;
	params = param;

	/* Fill in the first four steps with Runge Kutta. */

	x = x0;
	StateVector k1, k2, k3, k4;

	for (index=0; index<4; index++) {
		k1 = func(x, t, params);
		k2 = func(x+(dt/2.)*k1, t+(dt/2.), params);
		k3 = func(x+(dt/2.)*k2, t+(dt/2.), params);
		k4 = func(x+dt*k3, t+dt, params);
		x = x + (dt/6.) * (k1 + 2*k2 + 2*k3 + k4);
		t = t + dt;
		history[index] = func(x, t, params);
	}
	index--;

}

StateVector ABMIntegrator::step() {
	StateVector fkm3 = history[(index + SIZE - 3) % SIZE];
	StateVector fkm2 = history[(index + SIZE - 2) % SIZE];
	StateVector fkm1 = history[(index + SIZE - 1) % SIZE];
	StateVector fk = history[index];

	StateVector pkp1 = x + (dt/24.) * (-9*fkm3 + 37*fkm2 - 59*fkm1 + 55*fk);

	StateVector fkp1 = function(pkp1, t+dt, params);

	x = x + (dt/24.) * (fkm2 - 5*fkm1 + 19*fk + 9*fkp1);

	index = (index + 1) % SIZE;
	history[index] = fkp1;


	//cout << "err:" << endl;
	//cout << (-19./270. * (x - pkp1)).array()/x.array() << endl;;
	//cout << endl << endl;
	// Error: -19./270. * (x - pkp1);

	t += dt;

	return x;
}
