#include <cstdlib>
#include <sstream>
#include <ctime>
#include "simGRO.h"
#include "arguments.h"

using namespace std;

int main(int argc, char *argv[])
{
	// Read program arguments
	string prefFN;
	int eq_time=600;
	string outFN;
	int bin=10;
	int deltaT = 50;
	float eq_speed = 5.0;

	arg::push("-p","Preference file name",prefFN);
	arg::push("-o","Output file name",outFN);
	arg::push("-e","Equillibrium time (default=600s)",eq_time,true);
	arg::push("-s","Equillibrium speed (default=5x)",eq_speed,true);
	arg::push("-b","Bin size (default=10)",bin,true);
	arg::push("-r","Simulation resolution (default=50 ms)",deltaT,true);
	if(!arg::get(argc, argv)) return 0;
	
	// Initialize simulator
	InitSimGRO((char *)prefFN.c_str(), (char *)outFN.c_str());

	// Set timer
	clock_t start_time, cur_time;
	int elapsedTime, prevTime = 0, prevRepTime = 0;
	start_time = clock();
	
	// Equalibrate
	int totalEqTime = 0;
	while(totalEqTime < eq_time*1000) {
		simTime = 0; 					// Force simulated time to be zero
		UpdateSimGRO(deltaT*eq_speed);			// Update simulator every clock cycle
		totalEqTime += deltaT*eq_speed;			// Equillibration timer
		cur_time = clock();				// Real time tracing
		elapsedTime	= (int)((float)(cur_time - start_time)/CLOCKS_PER_SEC*1000);
		// Report status every 2 seconds
		if(elapsedTime>prevRepTime+2000) {
			cout<<"Equillibrating "<<(float)(totalEqTime/100)/10<<" seconds"<<endl;
			PrintStatus(elapsedTime, elapsedTime - prevTime);
			prevRepTime=elapsedTime;
		}
		prevTime=elapsedTime;
	}
	
	// Simulate
	simTime = 0;
	while(true) {
		UpdateSimGRO(deltaT);			// Update simulator
		cur_time = clock();				// Trace time
		elapsedTime	= (int)((float)(cur_time - start_time)/CLOCKS_PER_SEC*1000);
		Print(elapsedTime, elapsedTime - prevTime, bin);	// Report status and record
		prevTime = elapsedTime;
		if(currentRecord==countRecord) break;	// Simulate until the last recording
	}
	return 1;
}

