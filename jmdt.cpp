#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <string>
#include <queue>

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <GeographicLib/GravityModel.hpp>
#include "nrlmsise-00/nrlmsise-00.h"

#include "util.h"
#include "integrator.h"
#include "atmosphere.h"
#include "properties.h"

using namespace std;
using namespace Eigen;
using namespace GeographicLib;

typedef Matrix<double, 13, 1> StateVector;

/* 	[0] = x position in ECEF frame [m]
		[1] = y position in ECEF frame [m]
		[2] = z position in ECEF frame [m]
		[3] = x velocity in ECEF frame [m/s]
		[4] = y velocity in ECEF frame [m/s]
		[5] = z velocity in ECEF frame [m/s]
		[6] = Quaternion rotation angle [rad]
		[7] = Quaternion vector x component [unitless]
		[8] = Quaternion vector y component [unitless]
		[9] = Quaternion vector z component [unitless]
		[10] = Satellite frame x axis angular velocity [rad/s]
		[11] = Satellite frame y axis angular velocity [rad/s]
		[12] = Satellite frame z axis angular velocity [rad/s]

		See:
			S. D'Amico's AA279A notes in SSI Google Drive - ECEF discussed in Lecture 5
			https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
			https://en.wikipedia.org/wiki/Quaternion
			http://run.usc.edu/cs520-s14/quaternions/quaternions-cs520.pdf
*/

StateVector func(StateVector x, double t, IntegratorParams* params) {
	StateVector out;

	/* TEMPORARY */
	out[6] = 0.0;
	out[7] = 0.0;
	out[8] = 0.0;
	out[9] = 0.0;
	out[10] = 0.0;
	out[11] = 0.0;
	out[12] = 0.0;

	/* Velocity. */
	out[0] = x[3];
	out[1] = x[4];
	out[2] = x[5];

	Vector3d rvec(x[0], x[1], x[2]);
	double rmag = rvec.norm();

	Vector3d r_es = earth_sun_vector(params->t0 + t/86400.);
	Vector3d rsun = r_es - rvec;

	/* Gravity. */
	switch (params->earth) {
	case 1:
	case 2:
		params->gravity_object->V(x[0], x[1], x[2],
					  out[3], out[4], out[5]);
		break;
	case 0:
	default:
		double factor = -EARTH_MU/pow(rmag, 3);
		out[3] = factor * x[0];
		out[4] = factor * x[1];
		out[5] = factor * x[2];
	}

	Vector3d vec(x[3], x[4], x[5]);

	if (params->drag != 0 || params->power != 0) {
		if ((params->orientation_str[0] == 'c') && (params->lover != NULL)) {
			StateVector rtgtv = *(params->lover);
			Vector3d diff(rtgtv[0]-rvec[0], rtgtv[1]-rvec[1], rtgtv[2]-rvec[2]);
			if (diff.squaredNorm() > 100000.0*100000.0) {
				params->orientation_str = "p";
			}
		}
		if ((params->orientation_str == "t")
			&& (params->lover != NULL)) {
			StateVector rtgt = *(params->lover);
			Vector3d tgt(rtgt[0], rtgt[1], rtgt[2]);
			params->orientation = tgt-rvec;
		} else if ((params->orientation_str == "z") || ((params->orientation_str == "c1") && (params->lover != NULL))) {
			params->orientation = rvec;
		} else if (params->orientation_str == "n") {
			params->orientation = -rvec;
		} else if ((params->orientation_str == "p") || ((params->orientation_str == "c2") && (params->lover != NULL))) {
			params->orientation = vec;
		} else if (params->orientation_str == "r") {
			params->orientation = -vec;
		} else { // look at the Sun
			params->orientation = rsun;
		}
	}

	Matrix<double, 3, 3> R;
	if ((params->atmosphere != 0 && params->drag != 0) ||
		params->power == 1) {
		/* For solar panel and realistic drag, we need the satellite
		 * to face head on, using the nominal orientation according to
		 * which the GPU script generates the model. */
		const Vector3d b(1, 0, 0);
		params->orientation.normalize();

		Vector3d v = params->orientation.cross(b);
		double s = v.norm();
		double c = params->orientation.dot(b);

		R = MatrixXd::Identity(3, 3);

		if (fabs(c + 1) > 0.00001) {
			Matrix<double, 3, 3> vx;
			vx << 0, -v[2], v[1], v[2], 0, -v[0], -v[1], v[0], 0;

			R = R + vx + vx*vx*(1./(1+c));
		}
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

		Vector3d vrel = vec - EARTH_OMEGA.cross(rvec);

		switch (params->atmosphere) {
		case 2:
			nrlmsise_input inp;
			inp.year = 2016; // ignored by the library
			inp.doy = fmod(params->doy + (t/86400.0), 365.25);
			inp.sec = fmod(params->sec + t, 86400);

			double lat, lon;
			EARTH.Reverse(x[0], x[1], x[2], lat, lon, h);
			inp.alt = h/1000.;
			inp.g_lat = lat;
			inp.g_long = lon;
			inp.lst = inp.sec/3600 + inp.g_long/15;

			/* TODO: Use space weather data. */
			// Consider http://omniweb.gsfc.nasa.gov/form/dx1.html
			inp.f107A = params->f107A;
			inp.f107 = params->f107;
			inp.ap = params->ap;

			nrlmsise_flags flags;
			for (int i=0;i<24;i++) {
  				flags.switches[i]=1;
  			}
  			nrlmsise_output output;
  			gtd7d(&inp, &flags, &output);
  			rho = output.d[5];
  			break;
		case 1:
		default:
			rho = density_us1976(h);
		}

			Vector3d vrelp = R*(-vrel.normalized());
			double theta = acos(vrel.normalized()[2]);
			double phi = atan2(vrelp[1], vrelp[0]);
			if (phi < 0.0) {
				phi = 2*M_PI + phi;
			}
			/*if (phi > M_PI) {
				phi = 2*M_PI - phi;
			}*/
		double area;
		if (params->drag == 0) {
			area = params->A;
		} else {

			area = params->properties->get_drag_area(theta, phi);
		}

		Vector3d drag = -0.5*params->Cd*area*rho*
					vrel.squaredNorm()*vrel.normalized();
		drag = drag/params->mass;

		out[3] += drag[0];
		out[4] += drag[1];
		out[5] += drag[2];

		params->output_BC = drag.norm();//params->mass/(params->Cd*area);
	}

	/* Simulate solar panels and power. */
	if (params->power == 1) {
		if (!satellite_in_shade(rvec, -r_es)) {
			/* Vector love triangle, shout-out to SSP :'). */
			double dist2 = rsun.squaredNorm();

			/* Simulation works with the Sun vector that "hits" the
			 * satellite, not with the satellite-Sun vector. */
			rsun = -rsun;
			rsun.normalize();

			Vector3d rsunp = R*rsun;
			double theta = acos(rsunp[2]);
			double phi = atan2(rsunp[1], rsunp[0]);
			if (phi < 0.0) {
				phi = 2*M_PI + phi;
			}

			double area = params->properties->get_solar_area(theta,
								phi);
			params->output_power = params->solar_efficiency *
							area *
							SOLAR_CONSTANT *
							dist2/AU/AU;
		} else {
			params->output_power = 0.0;
		}

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
	string drag_file;
	double solar_efficiency;
	int two_satellites, separation_target;
	double x2, y2, z2, vx2, vy2, vz2;
	string first_orientation;
	string second_orientation;
	double f107A, f107, ap;
	int differential_drag;
	double k_i_drag, k_p_drag;
	double threshold_drag;
	int propulsion;
	double k_i_prop, k_p_prop;
	double threshold_prop;

	cin >> dt >> report_steps >> atmosphere >> earth >>
		x >> y >> z >> vx >> vy >> vz >> t0 >> tf >> output_size >>
		Cd >> A >> mass >> power >> solar_file >> drag_file >>
		solar_efficiency >> two_satellites >> separation_target >>
		x2 >> y2 >> z2 >> vx2 >> vy2 >> vz2 >> first_orientation >>
		second_orientation >> f107A >> f107 >> ap >>
		differential_drag >> k_i_drag >> k_p_drag >> threshold_drag >>
		propulsion >> k_i_prop >> k_p_prop >> threshold_prop;

	StateVector x0;
	x0 << x, y, z, vx, vy, vz, 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0;

	StateVector second_x0;
	second_x0 << x2, y2, z2, vx2, vy2, vz2 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0;

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
	params.drag = (drag_file == "none") ? 0 : 1;
	params.properties = new SatelliteProperties(solar_file, drag_file);
	params.orientation << 1, 0, 0;
	params.orientation_str = first_orientation;
	params.solar_efficiency = solar_efficiency;
	params.t0 = t0;
	params.f107A = f107A;
	params.f107 = f107;
	params.ap = ap;
	params.two_satellites = two_satellites;
	params.i_am_leader = true;
	params.separation_target = separation_target;
	params.controller_behavior = 0;
	params.lover = NULL;
	params.filter = 0.0;

	double year, month, day, hour, minute, second;
	//jd_to_date(t0, year, month, day, hour, minute, second);
	//params.doy = month*30.44 + day; // more or less...
	//params.sec = hour*3600 + minute*60 + second;

	IntegratorParams second_params;
	if (two_satellites == 1) {
		second_params = params;

		// No need to simulate power for the second one
		second_params.power = 0;

		// Second satellite is follower
		second_params.i_am_leader = false;

		second_params.orientation_str = second_orientation;
	}

	clock_t start = clock();
	ABMIntegrator integrator(func, x0, dt, &params);

	ABMIntegrator* second_integrator;
	if (two_satellites == 1) {
		second_integrator = new ABMIntegrator(func, second_x0,
							dt, &second_params);
		params.lover = &second_integrator->x;
		second_params.lover = &integrator.x;
	}

	int steps = 0;

	// Yeah, I took CS106X
	queue<double> dist_deriv_lpf;
	int lpf_max_samples = 5000;
	double lpf_sum = 0;

	while (integrator.t < tf) {
		integrator.step();
		if (two_satellites == 1) {
			second_integrator->step();
		}

		// Begin world's jankiest janktroller
		if(differential_drag){
			if(params.two_satellites == 1){
				Vector3d rvec_lead(integrator.x[0], integrator.x[1], integrator.x[2]);
				Vector3d rvec_follow(second_integrator->x[0], second_integrator->x[1], second_integrator->x[2]);
				Vector3d diff(rvec_lead[0]-rvec_follow[0], rvec_lead[1]-rvec_follow[1], rvec_lead[2]-rvec_follow[2]);

				Vector3d vec_lead(integrator.x[3], integrator.x[4], integrator.x[5]);
				Vector3d vec_follow(second_integrator->x[3], second_integrator->x[4], second_integrator->x[5]);
				Vector3d diff_vec(vec_lead[0]-vec_follow[0], vec_lead[1]-vec_follow[1], vec_lead[2]-vec_follow[2]);

				// See http://math.stackexchange.com/questions/1481701/time-derivative-of-the-distance-between-2-points-moving-over-time

				Vector3d bearing(rvec_lead[0]-rvec_follow[0], rvec_lead[1]-rvec_follow[1], rvec_lead[2]-rvec_follow[2]);
				bearing.normalize();
				double dist_deriv = (bearing[0]*diff_vec[0] + bearing[1]*diff_vec[1] + bearing[2]*diff_vec[2]);

				// Low Pass Filtering

				dist_deriv_lpf.push(dist_deriv);
				lpf_sum += dist_deriv;

				if(dist_deriv_lpf.size() > lpf_max_samples){
					lpf_sum -= dist_deriv_lpf.front();
					dist_deriv_lpf.pop();
				}

				double dist_deriv_filtered = lpf_sum / dist_deriv_lpf.size();

				double err_i = diff.norm() - params.separation_target; // Integral error term
				double err_p = dist_deriv_filtered; // Proportional error term, which is just velocity

				double controller_response = (k_i_drag * err_i) + (k_p_drag * err_p);

				if(controller_response < -threshold_drag){
					params.orientation_str = 'r';
					second_params.orientation_str = 'n';
					params.controller_behavior = 1;
					second_params.controller_behavior = 0;
				}else if(controller_response > threshold_drag){
					params.orientation_str = 'n';
					second_params.orientation_str = 'r';
					params.controller_behavior = 0;
					second_params.controller_behavior = 1;
				}else{
					params.orientation_str = 'n';
					second_params.orientation_str = 'n';
					params.controller_behavior = 0;
					second_params.controller_behavior = 0;
				}
			}
		}
		// End janktroller



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
			output[arg0+8] = params.output_BC;
			if (two_satellites == 1) {
				output[arg0+9] = second_integrator->x[0];
				output[arg0+10] = second_integrator->x[1];
				output[arg0+11] = second_integrator->x[2];
				output[arg0+12] = second_integrator->x[3];
				output[arg0+13] = second_integrator->x[4];
				output[arg0+14] = second_integrator->x[5];
				output[arg0+15] = second_params.output_power;
				output[arg0+16] = second_params.output_BC;
			}
			output[arg0+17] = params.controller_behavior;
			output[arg0+18] = second_params.controller_behavior;
			output[arg0+19] = params.filter;
			output[arg0+20] = second_params.filter;
		}
		steps++;
	}

	cerr << "Time elapsed: ";
	double elapsed = (clock()-start)/((double) CLOCKS_PER_SEC)*1000;
	cerr << elapsed << " ms (" << integrator.t/(elapsed/1000.)/86400.
			<< " days per second)" << endl;
	cerr << "Number of steps: " << steps << endl;


}
