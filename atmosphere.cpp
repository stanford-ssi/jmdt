#include <cmath>

#include "atmosphere.h"

double density_us1976(double h) {
	if (h > 1e6) return 0.0;
	else if (h < 0) return us1976[0];

	h = h/1000.0; // Meters to kilometers
	int a = h;
	int b = a+1;

	int da = exp(us1976[a]);
	int db = exp(us1976[b]);

	return da + (db-da)*(h-a);
}
