#include<cstdio>
#include<cstdlib>
#include<cmath>
#include<vector>
#include<iomanip>
#include<fstream>
#include<iostream>


#include"md_fec.cpp"		// plik z definicjami klas i procedurami numerycznymi

using namespace std;




int main(){
		
	double xmax=25.4;
	double ymax=26.1;
	double zmax=26.1;

	double fgrav=0.000;

	double dt=0.01;
	int nit_evol=4E+6;

	double vbar_x=0.;
	double vbar_y=0.;
	double vbar_z=0.0;
	double delta_bar=1.0;

	MD_FEC ob;
	ob.xmax=xmax;
	ob.ymax=ymax;
	ob.zmax=zmax;
	ob.fgrav=fgrav;
	ob.dt=dt;
	ob.nit_evol=nit_evol;
	ob.vbar_x=vbar_x;
	ob.vbar_y=vbar_y;
	ob.vbar_z=vbar_z;
	ob.delta_bar=delta_bar;

	ob.use_thermostat=false;
	ob.init_from_json("./input/particles_initial.json");

	ob.evolution();
	
	
	

	
	
}
