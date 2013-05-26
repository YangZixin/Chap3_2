#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <fstream>

using namespace std;
const double r = 160.7;
const double sigma = 10;
const double b = 8/3;
const double delta_t = 0.0001;
const int STEP = 500000;

struct coordinates
{
	double t;
	double x;
	double y;
	double z;
} coordinate={0,1,0,0};

struct velocities
{
	double vx;
	double vy;
	double vz;
} velocity;

void EulerMethod()
{
	coordinate.x += velocity.vx * delta_t; 
	coordinate.y += velocity.vy * delta_t; 
	coordinate.z += velocity.vz * delta_t;

	velocity.vx = sigma * (coordinate.y - coordinate.x);
	velocity.vy =(- coordinate.x * coordinate.z + r * coordinate.x - coordinate.y);
	velocity.vz = (coordinate.x * coordinate.y - b * coordinate.z);
	
	coordinate.t += delta_t;	 
}

void SORungeKutta()
{
	coordinates tempC;
	velocities tempV;
	
	tempC.x = coordinate.x + velocity.vx * delta_t / 2;
	tempC.y = coordinate.y + velocity.vy * delta_t / 2;
	tempC.z = coordinate.z + velocity.vz * delta_t / 2;
	tempC.t = coordinate.t + delta_t / 2;
	
	tempV.vx = sigma * (tempC.y - tempC.x);
	tempV.vy =(- tempC.x * tempC.z + r * tempC.x - tempC.y);
	tempV.vz = (tempC.x * tempC.y - b * tempC.z);

	coordinate.x += tempV.vx * delta_t;
	coordinate.y += tempV.vy * delta_t;
	coordinate.z += tempV.vz * delta_t;
	coordinate.t += delta_t;	 

	velocity.vx = sigma * (coordinate.y - coordinate.x);
	velocity.vy =(- coordinate.x * coordinate.z + r * coordinate.x - coordinate.y);
	velocity.vz = (coordinate.x * coordinate.y - b * coordinate.z);
}

void FORungeKutta()
{
	velocity.vx = sigma * (coordinate.y - coordinate.x);
	velocity.vy =(- coordinate.x * coordinate.z + r * coordinate.x - coordinate.y);
	velocity.vz = (coordinate.x * coordinate.y - b * coordinate.z);

	coordinates tempC1, tempC2, tempC3;
	velocities tempV1, tempV2, tempV3;

	tempC1.x = coordinate.x + velocity.vx * delta_t / 2;
	tempC1.y = coordinate.y + velocity.vy * delta_t / 2;
	tempC1.z = coordinate.z + velocity.vz * delta_t / 2;
	tempC1.t = coordinate.t + delta_t / 2;

	tempV1.vx = sigma * (tempC1.y - tempC1.x);
	tempV1.vy =(- tempC1.x * tempC1.z + r * tempC1.x - tempC1.y);
	tempV1.vz = (tempC1.x * tempC1.y - b * tempC1.z);

	tempC2.x = coordinate.x + tempV1.vx * delta_t / 2;
	tempC2.y = coordinate.y + tempV1.vy * delta_t / 2;
	tempC2.z = coordinate.z + tempV1.vz * delta_t / 2;
	tempC2.t = coordinate.t + delta_t / 2;

	tempV2.vx = sigma * (tempC2.y - tempC2.x);
	tempV2.vy =(- tempC2.x * tempC2.z + r * tempC2.x - tempC2.y);
	tempV2.vz = (tempC2.x * tempC2.y - b * tempC2.z);

	tempC3.x = coordinate.x + tempV2.vx * delta_t;
	tempC3.y = coordinate.y + tempV2.vy * delta_t;
	tempC3.z = coordinate.z + tempV2.vz * delta_t;
	tempC3.t = coordinate.t + delta_t;

	tempV3.vx = sigma * (tempC3.y - tempC3.x);
	tempV3.vy =(- tempC2.x * tempC3.z + r * tempC3.x - tempC3.y);
	tempV3.vz = (tempC2.x * tempC3.y - b * tempC3.z);

	coordinate.x += (velocity.vx + 2 * tempV1.vx + 2 * tempV2.vx + tempV3.vx) * delta_t / 6;
	coordinate.y += (velocity.vy + 2 * tempV1.vy + 2 * tempV2.vy + tempV3.vy) * delta_t / 6;
	coordinate.z += (velocity.vz + 2 * tempV1.vz + 2 * tempV2.vz + tempV3.vz) * delta_t / 6;
	coordinate.t += delta_t;	 

}



int main()
{
	ofstream myfile;

	/*myfile.open ("Euler_Method.dat");
	for(int i = 0; i < STEP; i++)
	{
		EulerMethod();
		myfile << coordinate.t <<"\t"<< coordinate.z <<endl;
	}
	myfile.close();

	myfile.open ("2-order_RK_Method.dat");
	for(int i = 0; i < STEP; i++)
	{
		SORungeKutta();
		myfile << coordinate.t <<"\t"<< coordinate.z <<endl;
	}
	myfile.close();
*/	
	myfile.open ("4-order_RK_Method_more.dat");
	for(int i = 0; i < STEP; i++)
	{
		FORungeKutta();
		myfile << coordinate.t << "    "<< coordinate.z <<endl;
	}

	myfile.close();
	return 0;
}
