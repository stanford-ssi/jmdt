#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <string>

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <GeographicLib/GravityModel.hpp>
#include "nrlmsise-00/nrlmsise-00.h"

#include "util.h"
#include "integrator.h"
#include "atmosphere.h"
#include "solar.h"

using namespace std;
using namespace Eigen;
using namespace GeographicLib;

// x y z vx vy vz

typedef Matrix<double, 6, 1> StateVector;

StateVector func(StateVector x, double t, IntegratorParams* params) {
	StateVector out;

	/* Velocity. */
	out[0] = x[3];
	out[1] = x[4];
	out[2] = x[5];

	Vector3d rvec(x[0], x[1], x[2]);
	double rmag = rvec.norm();
	
	/* Gravity. */
	switch (params->earth) {
	case 1:
	case 2:
		params->gravity_object->V(x[0], x[1], x[2],
					  out[3], out[4], out[5]);
	case 0:
	default:
		double factor = -EARTH_MU/pow(rmag, 3);
		out[3] = factor * x[0];
		out[4] = factor * x[1];
		out[5] = factor * x[2];
	}

	/* Drag. */
	if (params->atmosphere != 0) {
		double rho;
		
		/* Get geopotential height. Not implemented properly yet. */
		double h;
		switch (params->earth) {
		case 1: // GeographicLib would go here
		case 2: // GeographicLib would go here
		case 0:
		default:
			h = rmag - EARTH_RADIUS;
		}
		
		Vector3d vec(x[3], x[4], x[5]);
		Vector3d vrel = vec - EARTH_OMEGA.cross(rvec);

		switch (params->atmosphere) {
		case 1:
		default:
			rho = density_us1976(h);
		}
		
		Vector3d drag = -0.5*params->Cd*params->A*rho*
					vrel.squaredNorm()*vrel.normalized();
		drag = drag/params->mass;
		
		out[3] += drag[0];
		out[4] += drag[1];
		out[5] += drag[2];
	}

	/* Simulate solar panels and power. */
	if (params->power == 1) {
		Vector3d r_es = earth_sun_vector(params->t0 + t/86400.);
		
		/* Vector love triangle. */
		Vector3d rsun = r_es - rvec;

		/* Simulation works with the Sun vector that "hits" the
		 * satellite, not with the satellite-Sun vector. */
		rsun = -rsun;
		rsun.normalize();

		/* Rotate orientation so that the satellite faces (1, 0, 0),
		 * which is the nominal orientation that was used in the
		 * computation of sun incidence. Remember the rotation matrix
		 * and use it to transform the Sun-satellite vector; */
		const Vector3d b(1, 0, 0);
		params->orientation.normalize();

		Vector3d v = params->orientation.cross(b);
		double s = v.norm();
		double c = params->orientation.dot(b);

		Matrix<double, 3, 3> vx;
		vx << 0, -v[2], v[1], v[2], 0, -v[0], -v[1], v[0], 0;

		Matrix<double, 3, 3> R = MatrixXd::Identity(3,3) +
						vx + vx*vx*(1./(1+c));

		Vector3d rsunp = R*rsun;
		double theta = acos(rsunp[2]);
		double phi = atan2(rsunp[1], rsunp[0]);
		if (phi < 0.0) {
			phi = 2*M_PI + phi;
		}

		double area = params->solar_data->get_area(theta, phi);
		params->output_power = area;
	} else {
		params->output_power = 0.0;
	}

	return out;
}

double get_energy(StateVector& x) {
       return -EARTH_MU/pow(x[0]*x[0]+x[1]*x[1]+x[2]*x[2],0.5) +
		0.5*(x[3]*x[3]+x[4]*x[4]+x[5]*x[5]);
}

int main () {
	ios::sync_with_stdio(false);

	double dt;
	int report_steps, atmosphere, earth;
	double x, y, z, vx, vy, vz, t0, tf;
	int output_size;
	double Cd, A, mass;
	int power;
	string solar_file;

	cin >> dt >> report_steps >> atmosphere >> earth >>
		x >> y >> z >> vx >> vy >> vz >> t0 >> tf >> output_size >>
		Cd >> A >> mass >> power >> solar_file;

	StateVector x0;
	x0 << x, y, z, vx, vy, vz;

	/* Well, that was quite stupid. So at first I was like, oh, I want to
	 * communicate with Python! I have the power of Unix, I'll use pipes
	 * and a subprocess! Right? Well as it turns out cout is really slow -
	 * it was making everything a few times slower than without printing.
	 * Which was sad. So I tried to instead compile this as a library, and
	 * run it from Python with ctypes, and it was working fine, but for
	 * some reason it was still really slow. My suspicion is that there
	 * couldn't be as many optimizations because of position independent 
	 * code or something like that. So instead I shared memory between a 
	 * Python array and C++ with mmap, which I had never used before, and
	 * after some pain, achieve really fast IPC. Why, you might ask?
	 * 
	 * To get a few milliseconds of performance. That's why. Argh! */
	FILE* fmap = fopen("output.mmap", "r+");
	cout << 8*output_size*ceil(tf/dt/report_steps) << endl;
	double* output = (double*) mmap(0,
					8*output_size*ceil(tf/dt/report_steps),
					PROT_WRITE, MAP_SHARED,
					fileno(fmap), 0);

	IntegratorParams params;
	params.atmosphere = atmosphere;
	params.earth = earth;
	if (earth == 1) {
		params.gravity_object = new GravityModel("wgs84");
	} else if (earth == 2) {
		params.gravity_object = new GravityModel("egm96");
	}
	params.Cd = Cd;
	params.A = A;
	params.mass = mass;
	params.power = power;
	if (power == 1) {
		params.solar_data = new SolarData(solar_file);
	}
	params.orientation << 1, 0, 0;

	params.t0 = t0;

	clock_t start = clock();
	ABMIntegrator integrator(func, x0, dt, &params);

	int steps = 0;
	while (integrator.t < tf) {
		integrator.step();

		if (steps % report_steps == 0) {
			int arg0 = output_size*steps/report_steps;
			output[arg0+0] = integrator.t;
			output[arg0+1] = integrator.x[0];
			output[arg0+2] = integrator.x[1];
			output[arg0+3] = integrator.x[2];
			output[arg0+4] = integrator.x[3];
			output[arg0+5] = integrator.x[4];
			output[arg0+6] = integrator.x[5];
			output[arg0+7] = params.output_power;
		}
		steps++;
	}

	cerr << "Time elapsed: ";
	double elapsed = (clock()-start)/((double) CLOCKS_PER_SEC)*1000;
	cerr << elapsed << " ms (" << integrator.t/(elapsed/1000.)/86400.
			<< " days per second)" << endl;
	cerr << "Number of steps: " << steps << endl;


}
