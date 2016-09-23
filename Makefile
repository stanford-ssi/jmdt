all:
	g++ -march=native -Ofast `pkg-config --cflags eigen3` nrlmsise-00/nrlmsise-00.cpp nrlmsise-00/nrlmsise-00_data.c integrator.cpp atmosphere.cpp properties.cpp util.cpp jmdt.cpp -lgsl -lgslcblas -lm -lGeographic -o jmdt
