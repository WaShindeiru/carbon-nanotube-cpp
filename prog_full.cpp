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
		
	int N=500;
	double xmax=22.;
	double ymax=22.;
	double zmax=22.;
	
	double fgrav=0.005;	//przyspieszenie grawitacyjne w kierunku "z"
	
	int N_sim_anneal=1000;	//liczba iteracji w symulowanym wyrzazaniu - przygotowanie ukladu
	
	double dt=0.05;		//krok czasowy
	int nit_evol=1E+6;	//liczba iteracji
	
	
	//UZYWAMY PERIODYCZNE WB wiec w x/y nie ma barier
	double vbar_x=0.;	//potencjal sklaujacy na barierach w kierunku x
	double vbar_y=0.;	//potencjal sklaujacy na barierach w kierunku y
	double vbar_z=0.0;	//potencjal sklaujacy na barierach w kierunku z
	double delta_bar=1.0;	//zasieg potencjalu barierowego
	
	
	MD_FEC ob;
	ob.N=N;
	ob.xmax=xmax;
	ob.ymax=ymax;
	ob.zmax=zmax;
	ob.fgrav=fgrav;
	ob.N_sim_anneal=N_sim_anneal;
	ob.dt=dt;
	ob.nit_evol=nit_evol;
	
	
	ob.vbar_x=vbar_x; 	
	ob.vbar_y=vbar_y; 	
	ob.vbar_z=vbar_z; 	
	ob.delta_bar=delta_bar;	
	
	ob.init();
	
	ob.evolution();
	
	
	

	
	
}
