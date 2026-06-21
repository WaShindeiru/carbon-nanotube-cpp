/************************************************
*
*
*
*************************************************/
#include<cmath>
#include<vector>
#include<cstdio>
#include<cstdlib>
#include<iomanip>
#include<fstream>
#include<iostream>
using namespace std;

class MD_FEC{
	

public: class DATA{
		public:
		vector<double> D0={1.5, 6.0, 4.82645134};
		vector<double> r0={2.29, 1.39, 1.47736510};	
		vector<double> beta={1.4, 2.1, 1.63208170};
		vector<double> S={2.0693109, 1.22, 1.43134755};
		vector<double> gamma={0.0115751, 2.0813E-4, 0.00205862};
		vector<double> c={1.2898716, 330.0, 8.95583221};
		vector<double> d={0.3413219, 3.5, 0.72062047};
		vector<double> h={-0.26, 1.0, 0.87099874};
		vector<double> R={3.15, 1.85, 2.5};
		vector<double> D={0.2, 0.15, 0.2};
		vector<double> rf={0.95, 0.6, 1.0};
		vector<double> bf={2.9, 8.0, 10.};
		vector<double> mass={55.845, 12.0};
		
		double kb=1.38E-23;
		double mass_unit=1.66E-27;
		double length_unit=1.0E-10;
		double energy_unit=1.602E-19;
		double temperature_unit=energy_unit/kb;
		double velocity_unit=sqrt(energy_unit/mass_unit);
		double force_unit=energy_unit/length_unit;
		double acceleration_unit=energy_unit/mass_unit/length_unit;
		double time_unit=sqrt(mass_unit/energy_unit)*length_unit;
		
	};

public:  class PARTICLE{
		public:
		vector<double> r_vec={0.,0.,0.};
		vector<double> v_vec={0.,0.,0.};
		vector<double> vxp_vec={0.,0.,0.};
		vector<double> f_vec={0.,0.,0.};
		vector<double> tmp_vec={0.,0.,0.};
		vector<int> cell={-1,-1,-1};		//polozenie w komorce obliczeniowej
		double vpot;
		int type=-1;  //type of particle: 0-Fe, 1-C (must be set)
	};
	

public:
	DATA data;	
	
	int N_sim_anneal=100; //liczba iteracji w symulowanym wyzarzaniu (szukamy optymalnego rozkladu poczatkowego czastek)
	double temperature_anneal=0.1; //liczba iteracji w symulowanym wyzarzaniu (szukamy optymalnego rozkladu poczatkowego czastek)
	double temperature_reference=0.13;
	double fgrav=0.0; //sila grawitacyjna 
	double vbar_x=0.;		//skalowanie potencjalu oddzialywania z barierami
	double vbar_y=0.;		//skalowanie potencjalu oddzialywania z barierami
	double vbar_z=0.;		//skalowanie potencjalu oddzialywania z barierami
	double delta_bar=0.5; 	//zasieg oddzialywania z barierami
	
	int N=-1;	//liczba czastek
	double xmax=-1,ymax=-1,zmax=-1;
	double delta_cell_min;
	double dx_cell,dy_cell,dz_cell;
	int mx_cell,my_cell,mz_cell;

	double dt=0.01;
	int nit_evol=-1;
	double vpot_tot=0;
	double ekin_tot=0;
	double force_max;
	
	vector<PARTICLE> particle;
	vector<vector<vector<int>>> header;  //tablica komorek do metody sortowania linked-cell
	vector<int> link;				//tablica indeksow czastek do metody sortowania linked-cell
	
	
public:	
	
	void init(){
		DATA data;
		particle.resize(N);
		
		// okreslamy minimalna szerokosc komorki na podstawie zakresu oddzialywania: R+D
		// szerokosc komorki moze byc wieksza
		delta_cell_min=0;
		for(int i=0;i<2;i++){
			delta_cell_min=max(delta_cell_min,data.D[i]+data.R[i]);
		}
		
		mx_cell=int(xmax/delta_cell_min);
		my_cell=int(ymax/delta_cell_min);
		mz_cell=int(zmax/delta_cell_min);
		dx_cell=xmax/mx_cell;
		dy_cell=ymax/my_cell;
		dz_cell=zmax/mz_cell;


		cout<<"================================================="<<endl;
		cout<<"N = "<<setw(20)<<N<<endl;
		cout<<"xmax = "<<setw(20)<<xmax<<endl;
		cout<<"ymax = "<<setw(20)<<ymax<<endl;
		cout<<"zmax = "<<setw(20)<<zmax<<endl;
		
		cout<<"mx_cell = "<<setw(20)<<mx_cell<<endl;
		cout<<"my_cell = "<<setw(20)<<my_cell<<endl;
		cout<<"mz_cell = "<<setw(20)<<mz_cell<<endl;
		cout<<"dx_cell = "<<setw(20)<<dx_cell<<endl;
		cout<<"dy_cell = "<<setw(20)<<dy_cell<<endl;
		cout<<"dz_cell = "<<setw(20)<<dz_cell<<endl;

		cout<<"================================================="<<endl;
		
		
		//tworzymy tablice do sortowania
		header.resize(mx_cell,vector<vector<int>>(my_cell,vector<int>(mz_cell)));
		link.resize(N);
		
		//okreslamy typ czastek
		for(int i=0;i<N;i++){
			if(i<=N/4)particle[i].type=1;
			else	particle[i].type=0;
		}
		// rozmieszczamy czastki w pudle losowo - pozniej  optymalizacja "simulated annealing"
		init_particle_distribution();
				
		
		return;
	}
	
	
	
	
	
/************************************************************************************
 * ewolucja czasowa
 *
 ************************************************************************************/	
	void evolution(){
		
		double vtot;
		double ekin;
		double etot;
		double temp;
		double ekin_1;
		double temp_avg;
		double temp_d=1500./data.temperature_unit;
		int ksr;
		double gamma_temp;
		int ile_print;
		
		//przed pierwszym wywolaniem verleta trzeba posortowac i policzyc sily
		sort();
		compute_forces();
		compute_forces_barrier_grav();
		
		double t=0.;
		
				
		ksr=0;
		ekin_1=0.;
		ile_print=0;
		
		for(int it=0;it<=nit_evol;it++){
		
			
			
			ksr++;
			ekin_1+=ekin_tot;  //compute_total_kinetic_energy();
			temp_avg=2./3.*ekin_1/N/ksr;
			gamma_temp=sqrt(temp_d/temp_avg);
			
			
			
			if(it%500==0){
				ile_print++;
				
				double vmax=0.;
				for(int i=0;i<N;i++){
					double vx=particle[i].v_vec[0];	
					double vy=particle[i].v_vec[1];
					double vz=particle[i].v_vec[2];
					double v=sqrt(vx*vx+vy*vy+vz*vz);
					
					vmax=max(vmax,v);
					
				}
				
				
				ekin=ekin_tot;  //ekin_tot=compute_total_kinetic_energy();
				vtot=vpot_tot;  //vpot_tot=compute_vtot();
				etot=ekin+vtot;
				temp=2./3.*ekin/N;
				
				cout<<setw(10)<< t <<"\t";
				cout<<setw(12)<< temp*data.temperature_unit <<"\t";
				cout<<setw(12)<< temp_avg*data.temperature_unit <<"\t";
				cout<<setw(12)<< vtot <<"\t";
				cout<<setw(12)<< ekin <<"\t";
				cout<<setw(12)<< etot <<"\t";
				cout<<setw(12)<< vmax <<"\t";
				cout<<setw(12)<< vmax*dt <<"\t";
				cout<<setw(12)<< force_max <<"\t";
				cout<<endl;
				write_positions_to_file("pos.dat");
				
				
				//zmiana temperatury
				if(ile_print>1 && ile_print<100){
					ekin_1=0;
					ksr=0;
					for(int i=0;i<N;i++){
						for(int j=0;j<3;j++)	particle[i].v_vec[j]*=gamma_temp;
					}
					cout<<ile_print<<"\t"<< temp_avg<<"\t" << temp_d<<"\t"<<gamma_temp<<endl;
				}
				
				
				
				 
			}
			
			// tu wykonujemy 1-krok VERLETEM
			verlet_one_step();
			
			t+=dt;	
						
		}
		
		
		return;
	}
	
	
	
	
	
	
	
	
	
/************************************************************************************
 *	velocity-Verlet: one-step
 * 	calculated: 
 * 			single-particle forces (storage: array particles)
 * 			single particle potential energy (storage: array particles)
 * 			total potential energy (vpot_tot)
 ************************************************************************************/
void verlet_one_step(){
	
	//calculations of: v(t+dt/2), r(t+dt)
	for(int i=0;i<N;i++){
		
		int itype=particle[i].type;
		double mass=data.mass[itype];
		
		//liczymy predkosc w dla t+dt/2
		add_vectors(3, 1.0,particle[i].v_vec, dt/2./mass,particle[i].f_vec, 0.0,particle[i].vxp_vec);
		
		//liczymy polozenie w t+dt
		add_vectors(3, 0.00000,particle[i].v_vec, dt,particle[i].vxp_vec, 1.0,particle[i].r_vec);

		//periodyczne WB - x/y (predkosc zachowujemy)
		if(particle[i].r_vec[0]<0)	particle[i].r_vec[0]+=xmax;
		if(particle[i].r_vec[0]>xmax) particle[i].r_vec[0]-=xmax;
		
		if(particle[i].r_vec[1]<0)	particle[i].r_vec[1]+=ymax;
		if(particle[i].r_vec[1]>ymax) particle[i].r_vec[1]-=ymax;

	//w kierunku z: odbicie wzgledem dolnej/gornej scianki niezaleznie od pozostalych WB
	//jesli wprowadzony jest potencjal barierorwy to ponizsze WB nie beda aktywowane

		if(particle[i].r_vec[2]<0){
			particle[i].r_vec[2]*=-1;
			particle[i].vxp_vec[2]*=-1;
		}
		if(particle[i].r_vec[2]>zmax){
			particle[i].r_vec[2]=zmax-(particle[i].r_vec[2]-zmax);
			particle[i].vxp_vec[2]*=-1;
		}

	
		
	}
	
	//sortowanie
	sort();
		
	//calculations of forces for: t+dt
	compute_forces();
	compute_forces_barrier_grav(); //do sztywnych WB
	
	//calculations of: v(t+dt)
	for(int i=0;i<N;i++){
		int itype=particle[i].type;
		double mass=data.mass[itype];
		add_vectors(3, 1.0,particle[i].vxp_vec, dt/2/mass,particle[i].f_vec, 0.0,particle[i].v_vec);
		
	}
	
	//energie kinetyczna liczymy w tej samej chwili co potencjal (t+dt)
	ekin_tot=compute_total_kinetic_energy();
	
	
	return;
}	
	
	
	
/***********************************************************************************
 *	wklad do sil i energii od barier i pola grawitacyjnego
 * 
 ***********************************************************************************/
void	compute_forces_barrier_grav()
{
	double epot_bar;
	double epot_grav;
	for(int i=0;i<N;i++){
		compute_single_energy_forces_barrier_grav(i,epot_bar,epot_grav);
		vpot_tot+=epot_bar+epot_grav;
	}
	
	
	
}
	
	
/***********************************************************************************
 *	1) wklad do sily dzialajacej na pojedyncza czastke od barier
 * 	2) wklad do energii od barier
 * 	3) wklad do sily grawitacyjnej
 * 	4) wklad do potencjalu grawitacyjnego
 *  
 * 	UWAGA: WKLADY DO SIL DODAWANE OD RAZU DO TABLICY PARTICLE  !!!!!!!
 * 
 * 		 ENERGIE ZWRACANE  (UZYWANE TEZ W SYMULOWANYM WYGRZEWANIU)
 * 
 * 	vbar_x,vbar_y, vbar_z - zmienne globalna
 * 	delta_bar -zmienna globalna
 * 
 ***********************************************************************************/
void	compute_single_energy_forces_barrier_grav(const int & i, double & epot_bar, double & epot_grav)
{
	
	//polozenie czastki
	double x=particle[i].r_vec[0];
	double y=particle[i].r_vec[1];
	double z=particle[i].r_vec[2];
	
	//oddzialywanie z barierarmi
	double fx=0.;	//skladowe sil pochodzacych od barier
	double fy=0.;
	double fz=0.;
	double de=0.;	//energia oddzialywania z barierami
	
/*
	if(x<delta_bar && fabs(vbar_x)>0.01){
		fx=-2*vbar_x/x*log(x/delta_bar);
		de+=vbar_x*pow(log(x/delta_bar),2);
	}
	
	if(x>(xmax-delta_bar) && fabs(vbar_x)>0.01){
		fx=2*vbar_x/(xmax-x)*log((xmax-x)/delta_bar);
		de+=vbar_x*pow(log((xmax-x)/delta_bar),2);
	}
	
	if(y<delta_bar && fabs(vbar_y)>0.01){
		fy=-2*vbar_y/y*log(y/delta_bar);
		de+=vbar_y*pow(log(y/delta_bar),2);
	}
			
	if(y>(ymax-delta_bar && fabs(vbar_y)>0.01)){
		fy=2*vbar_y/(ymax-y)*log((ymax-y)/delta_bar);
		de+=vbar_y*pow(log((ymax-y)/delta_bar),2);
	}
*/
	
	
	if(z<delta_bar && fabs(vbar_z)>0.01){
		fz=-2*vbar_z/z*log(z/delta_bar);
		de+=vbar_z*pow(log(z/delta_bar),2);
	}
	
	if(z>(zmax-delta_bar && fabs(vbar_z)>0.01)){
		fz=2*vbar_z/(zmax-z)*log((zmax-z)/delta_bar);
		de+=vbar_z*pow(log((zmax-z)/delta_bar),2);
	}
			
	//energia oddzialywania z barierami		
	epot_bar=de;
	//dodajemy sily od barier
	particle[i].f_vec[0]+=fx;
	particle[i].f_vec[1]+=fy;
	particle[i].f_vec[2]+=fz;
				
			
	//sila i energia grawitacyjna
	int itype=particle[i].type;
	double mass=data.mass[itype];
	double fgrav_z=-fgrav*mass;
	
	particle[i].f_vec[2]+=fgrav_z;
	epot_grav=fgrav*mass*z;  
	
	
	return;	
}
	
	
	
	
	
	
	
	
/************************************************************************************
 *	calculate: forces, total potential 
 ************************************************************************************/
void	compute_forces(){
	
	int n_neigh;
	vector<int> neigh_index(5000,0);  //tablica z indeksami sasiadow do liczenia potencjalu
	vector<vector<double>> neigh_r_vec(5000,vector<double>(3,0.));  //tablica z polozeniami sasiadow z uwzglednieniem periodycznych WB
	
	int type;
	double R,D,S,D0,beta,r0,gamma,c,d,h; //data
	
	vector<double> rij_vec(3.,0);
	vector<double> rik_vec(3.,0);
	vector<double> rj_vec(3.,0);
	vector<double> rk_vec(3.,0);
	vector<double> grad_fc_ij_ri(3.,0);
	vector<double> grad_va_ij_ri(3.,0);
	vector<double> grad_vr_ij_ri(3.,0);
	
	vector<double> grad_fc_ik_ri(3.,0);
	vector<double> grad_gik_ri(3.,0);
	vector<double> grad_gik_rj(3.,0);
	vector<double> grad_gik_rk(3.,0);
	
	
	vector<double> grad_bij_ri(3.,0);
	vector<double> grad_bij_rj(3.,0);
	vector<double> vec_empty(3.,0); 	//empty vector
	
	
	/************************************
	* reset force array for each atom
	*************************************/ 
	for(int i=0;i<N;i++){
		vector_zero(3,particle[i].f_vec);	
	}	
			
	vpot_tot=0.;

	for(int i=0;i<N;i++){
		int itype=particle[i].type;
		
		//szukamy sasiadow
		find_nearest_particles_and_positions(i,n_neigh,neigh_index,neigh_r_vec);
		
		/************************************
		* loop over FIRST neighbour: j-th atom
		*************************************/ 
		for(int l=0;l<n_neigh;l++){
			int j=neigh_index[l]; 			//global index of neighbour
			//define type of j-th atom and read data
			int jtype=particle[j].type;
			type=(2*(itype+jtype))%3; //0-FeFe,1-CC,2-FeC
			R=data.R[type];
			D=data.D[type];
			S=data.S[type];
			D0=data.D0[type];
			beta=data.beta[type];
			r0=data.r0[type];
			//calculatate quantities depending on: rij
			add_vectors(3, 1.0, neigh_r_vec[l], -1.0,particle[i].r_vec, 0., rij_vec);
			double rij=vector_norm(3,rij_vec);
			double fc_ij=compute_fc(rij,R,D);
			if(fc_ij<1.0E-15)continue;
			double va_ij=compute_va(rij,D0,S,beta,r0);
			double vr_ij=compute_vr(rij,D0,S,beta,r0);
			compute_grad_fc_ri(rij,R,D,rij_vec,grad_fc_ij_ri);
			compute_grad_vr_ri(rij,D0,S,beta,r0,rij_vec,grad_vr_ij_ri);
			compute_grad_va_ri(rij,D0,S,beta,r0,rij_vec,grad_va_ij_ri);
		
			/*****************************************
			* loop over SECOND neighbour: k-th atom
			******************************************/	
			vector_zero(3,grad_bij_ri);
			vector_zero(3,grad_bij_rj);
			
			double chi_ij=0.;
			
			for(int m=0;m<n_neigh;m++){
				int k=neigh_index[m]; 			//global index of neighbour
				if(k==i || k==j)continue;
			//define type of k-th atom and read data
				int ktype=particle[k].type;
				type=(2*(itype+ktype))%3; //0-FeFe,1-CC,2-FeC
				R=data.R[type];
				D=data.D[type];
				S=data.S[type];
				D0=data.D0[type];
				beta=data.beta[type];
				r0=data.r0[type];
				gamma=data.gamma[type];
				c=data.c[type];
				d=data.d[type];
				h=data.h[type];
			//calculate quantities depending on: rik, ijk
				add_vectors(3, 1.0, neigh_r_vec[m], -1.0,particle[i].r_vec, 0., rik_vec);
				double rik=vector_norm(3,rik_vec);
				double fc_ik=compute_fc(rik,R,D);
				//if(fc_ik<1.0E-15)continue;
				double cos_teta_ijk=compute_cos_teta_ijk(rij_vec,rij,rik_vec,rik);
				double gik=compute_gik(gamma,c,d,h,cos_teta_ijk);
				compute_grad_fc_ri(rik,R,D,rik_vec,grad_fc_ik_ri);
				compute_grad_gik(gamma,c,d,h,cos_teta_ijk,rij_vec,rij,rik_vec,rik,grad_gik_ri,grad_gik_rj,grad_gik_rk);
				
			//add contributions of grad_bij_ri 	
				add_vectors(3, gik, grad_fc_ik_ri, fc_ik, grad_gik_ri, 1.0, grad_bij_ri);
			//add contributions of grad_bij_rj 		
				add_vectors(3, 0.0, vec_empty, fc_ik, grad_gik_rj, 1.0, grad_bij_rj);
			
			//add contribution of grad_bij_rk  (-gik because:  grad_fc_ik_ri= -1*grad_fc_ik_rk)
				add_vectors(3, -gik, grad_fc_ik_ri, fc_ik, grad_gik_rk,0.0,particle[k].tmp_vec);
				
				
			//calculate chi_ij	
				chi_ij+=fc_ik*gik;
			}		
			
			double bij=1./sqrt(1.+chi_ij);
			double dbij_chi_ij=-pow(1+chi_ij,-1.5)/2.;
			
			
			//scale vectors with derivative of bij w.r.t. chi_ij 
			scale_vector(3,dbij_chi_ij, grad_bij_ri);
			scale_vector(3,dbij_chi_ij, grad_bij_rj);
			
			for(int m=0;m<n_neigh;m++){
				int k=neigh_index[m]; 			//global index of neighbour
				if(k==i || k==j)continue;
				scale_vector(3,dbij_chi_ij, particle[k].tmp_vec);		
			}
						
			
			//add contribution to vector F_i:
			for(int m=0;m<3;m++){
				particle[i].f_vec[m]-=0.5*(grad_fc_ij_ri[m]*(vr_ij-bij*va_ij)
								  +fc_ij*(grad_vr_ij_ri[m]-grad_bij_ri[m]*va_ij-bij*grad_va_ij_ri[m]) );
			}
			
			//add contribution to vector F_j (grad_rj <= -1*grad_i):
			for(int m=0;m<3;m++){
				particle[j].f_vec[m]-=0.5*(-grad_fc_ij_ri[m]*(vr_ij-bij*va_ij)
								  +fc_ij*(-grad_vr_ij_ri[m]-grad_bij_rj[m]*va_ij-bij*(-grad_va_ij_ri[m])) );
			}
			
			//add contribution to vector F_k
			for(int m=0;m<n_neigh;m++){
				int k=neigh_index[m]; 			//global index of neighbour
				if(k==i || k==j)continue;
				for(int m=0;m<3;m++){
					particle[k].f_vec[m]-=0.5*fc_ij*(-1.*particle[k].tmp_vec[m])*va_ij;
				}
			}
			
			vpot_tot+=0.5*fc_ij*(vr_ij-bij*va_ij);

		}
	}
	
	
	
	
	//sprawdzamy sume wszystkich sil w ukladzie
	double fx=0,fy=0,fz=0.;
	force_max=0.; //zmienna globalna
	
	for(int i=0;i<N;i++){
		fx+=particle[i].f_vec[0];
		fy+=particle[i].f_vec[1];
		fz+=particle[i].f_vec[2];
		double force=sqrt(pow(particle[i].f_vec[0],2)+pow(particle[i].f_vec[1],2)+pow(particle[i].f_vec[2],2));
		force_max=max(force_max,force);
	}
	if(fabs(fx)>1.0E-10  || fabs(fy)>1.0E-10  || fabs(fz)>1.0E-10 ){
		cout<<"nieskompensowane sily wewnetrzne"<<endl;
		cout<<fx<<endl;
		cout<<fy<<endl;
		cout<<fz<<endl;
	}
	
	

	
	
	
	return;
}	
	
	
	
	
	
	
	
/**************************************************************
 * zapisujemy polozenia czastek do pliku
 **************************************************************/	
	void write_positions_to_file(string filename){
		ofstream file2;
		file2.open(filename);
		file2.setf(ios::scientific|ios::right);
		file2.precision(10);
		for(int i=0;i<N;i++){
			file2<<setw(20)<< particle[i].r_vec[0]<<"\t";
			file2<<setw(20)<< particle[i].r_vec[1]<<"\t";
			file2<<setw(20)<< particle[i].r_vec[2]<<"\t";
			file2<<setw(4)<< particle[i].type<<"\t";
			file2<<setw(4)<< particle[i].cell[0]<<"\t";
			file2<<setw(4)<< particle[i].cell[1]<<"\t";
			file2<<setw(4)<< particle[i].cell[2]<<"\t";
			file2<<endl;
		}
		file2.close();
	
		return;
	}
	
	
/**************************************************************
 * generujemy poczatkowe polozenia czastek
 * 1: losowo rozkladamy w pudle
 * 2: metoda symulowanego wyzarzania oddalamy je od siebie
 * 3: na dole i na gorze czastki odbijane od scianek do srodka
***************************************************************/	
	void init_particle_distribution(){	
			
	
		vector<int> neigh_index(5000,0);  //tablica z indeksami sasiadow do liczenia potencjalu
		vector<vector<double>> neigh_r_vec(5000,vector<double>(3,0.));  //tablica z polozeniami sasiadow z uwzglednieniem periodycznych WB
		
		int n_neigh;
		
		//losowe rozlozenie czastek w pudle	
		for(int i=0;i<N;i++){
			particle[i].r_vec[0]=xmax*uniform();
			particle[i].r_vec[1]=ymax*uniform();
			particle[i].r_vec[2]=zmax*uniform();
		}
		
						
		//symulowane wyzarzanie - petla glowna		
		for(int it=0;it<=N_sim_anneal;it++){
			
			//sortowanie czastek
			sort();
	
			//energia oddzialywania
			double vtot=0.;
			for(int i=0;i<N;i++){
				//double vi_bar=0;
				//double vi_grav=0;
				//compute_single_energy_forces_barrier_grav(i,vi_bar,vi_grav);
				//find_nearest_particles_and_positions(i,n_neigh,neigh_index,neigh_r_vec);
				//vtot+=compute_potential_vi(i,n_neigh,neigh_index,neigh_r_vec);
				//vtot+=vi_grav;
			}

			compute_forces();
			compute_forces_barrier_grav();
			vtot+=vpot_tot;


				
			if(it%100==0){
				cout<<"simulated annealing, it, Vtot = "<<it<<"\t"<<vtot<<endl;	
				write_positions_to_file("xyz_init0.dat");				
			}
			
			//petla po czastkach
			for(int i=0;i<N;i++){
								
				//szukamy sasiadow czastki
				find_nearest_particles_and_positions(i,n_neigh,neigh_index,neigh_r_vec);
				
				//stare polozenia
				double xold=particle[i].r_vec[0];
				double yold=particle[i].r_vec[1];
				double zold=particle[i].r_vec[2];
				
				//liczymy potencjal oddzialywania + grawitacyjny
				//oddzialywanie z innymi czastkami
				double vi_old=compute_potential_vi(i,n_neigh,neigh_index,neigh_r_vec);				
				//oddzialywanie z barierami i polem grawitacyjnym
				double vi_grav_old=0.;
				double vi_bar_old=0.;
				compute_single_energy_forces_barrier_grav(i,vi_bar_old,vi_grav_old);
				
				//nowe polozenia: losowo przesuwamy ze starej pozycji
				double dshift_max=0.1;
				
				double xnew=xold+(2*uniform()-1.0)*dshift_max;
				double ynew=yold+(2*uniform()-1.0)*dshift_max;
				double znew=zold+(2*uniform()-1.0)*dshift_max;
				
				//uwzgledniamy periodyczne WB
				if(xnew<0)xnew+=xmax;
				if(xnew>xmax)xnew-=xmax;
				if(ynew<0)ynew+=ymax;
				if(ynew>ymax)ynew-=ymax;
				if(znew<0)znew=zold;
				if(znew>zmax)znew=zold;
				
				particle[i].r_vec[0]=xnew;
				particle[i].r_vec[1]=ynew;				
				particle[i].r_vec[2]=znew;				
				
				//oddzialywanie z innymi czastkami
				double vi_new=compute_potential_vi(i,n_neigh,neigh_index,neigh_r_vec);	
				//oddzialywanie z barierami i polem grawitacyjnym
				double vi_grav_new=0.;
				double vi_bar_new=0.;
				compute_single_energy_forces_barrier_grav(i,vi_bar_new,vi_grav_new);
				//algorytm Metropolisa: jesli nie akceptujemy nowego polozenia to przywracamy stare
				
				double eold=vi_old+vi_grav_old+vi_bar_old;
				double enew=vi_new+vi_grav_new+vi_bar_new;
				
				double pacc=exp(-(enew-eold)/temperature_anneal);
				if(uniform()>pacc){
					particle[i].r_vec[0]=xold;
					particle[i].r_vec[1]=yold;				
					particle[i].r_vec[2]=zold;				
				}
			}
		}
		
		
		return;
	}

	
	
	
	
	
/************************************************************************************
 *	kinetic energy:
 ************************************************************************************/
inline double compute_total_kinetic_energy(){
	double ekin_tot=0.;
	for(int i=0;i<N;i++){
		int type=particle[i].type;
		double mass=data.mass[type];
		double vx=particle[i].v_vec[0];
		double vy=particle[i].v_vec[1];
		double vz=particle[i].v_vec[2];
		ekin_tot+=mass*(vx*vx+vy*vy+vz*vz)/2;
	}
	return ekin_tot;	
}	
	
	
/**************************************************************
 * liczymy calkowity potencjal: vtot 
***************************************************************/	
	double compute_vtot(){	
	
		double vtot=0;
		vector<int> neigh_index(5000,0);  //tablica z indeksami sasiadow do liczenia potencjalu
		vector<vector<double>> neigh_r_vec(5000,vector<double>(3,0.));  //tablica z polozeniami sasiadow z uwzglednieniem periodycznych WB
		int n_neigh;
		for(int i=0;i<N;i++){
			double vi_grav=compute_potential_grav(particle[i].r_vec[2]);
			find_nearest_particles_and_positions(i,n_neigh,neigh_index,neigh_r_vec);		
			vtot+=0.5*compute_potential_vi(i,n_neigh,neigh_index,neigh_r_vec);
			vtot+=vi_grav;
		}	
		return vtot;
	}
	

/************************************************************************************
 *	potencjal grawitacyjny:  U = Fgrav*z
 ************************************************************************************/
	double compute_potential_grav(const double & z){
		return fgrav*z;
	}	
	
	



/************************************************************************************
 *	potencjal oddzialywania i-tej czastki
 ************************************************************************************/
	double compute_potential_vi(const int & i, const int & n_neigh, const vector<int> & neigh_index, const vector<vector<double>> & neigh_r_vec ){
	
	int type;
	double R,D,S,D0,beta,r0,gamma,c,d,h; //data
	
	vector<double> rij_vec(3.,0);
	vector<double> rik_vec(3.,0);
			
	double vi=0.;
	int itype=particle[i].type;
	
	/***************************************
	* loop over FIRST neighbour: j-th atom
	****************************************/
	for(int l=0;l<n_neigh;l++){
		int j=neigh_index[l]; 			//global index of neighbour
		//define type of j-th atom and read data
		int jtype=particle[j].type;
		type=(2*(itype+jtype))%3; //0-FeFe,1-CC,2-FeC
		R=data.R[type];
		D=data.D[type];
		S=data.S[type];
		D0=data.D0[type];
		beta=data.beta[type];
		r0=data.r0[type];
		//calculate quantities depending on: rij
		add_vectors(3, 1.0, neigh_r_vec[l], -1.0,particle[i].r_vec, 0., rij_vec);
		double rij=vector_norm(3,rij_vec);
		double fc_ij=compute_fc(rij,R,D);
		
		if(fc_ij<1.0E-15) continue;
		
		double va_ij=compute_va(rij,D0,S,beta,r0);
		double vr_ij=compute_vr(rij,D0,S,beta,r0);
		
		
		
		/*****************************************
		* loop over SECOND neighbour: k-th atom
		*****************************************/	
		double chi_ij=0.;
		
		for(int m=0;m<n_neigh;m++){
			int k=neigh_index[m]; 			//global index of neighbour
			if(k==i || k==j)continue;
			//define type of k-th atom and read data
			int ktype=particle[k].type;
			type=(2*(itype+ktype))%3; //0-FeFe,1-CC,2-FeC
			R=data.R[type];
			D=data.D[type];
			S=data.S[type];
			D0=data.D0[type];
			beta=data.beta[type];
			r0=data.r0[type];
			gamma=data.gamma[type];
			c=data.c[type];
			d=data.d[type];
			h=data.h[type];
			//calculate quantities depending on: rik, ijk
			add_vectors(3, 1.0, neigh_r_vec[m], -1.0,particle[i].r_vec, 0., rik_vec);
			double rik=vector_norm(3,rik_vec);
			double fc_ik=compute_fc(rik,R,D);
			double cos_teta_ijk=compute_cos_teta_ijk(rij_vec,rij,rik_vec,rik);
			double gik=compute_gik(gamma,c,d,h,cos_teta_ijk);	
			
			//calculate chi_ij	
			chi_ij+=fc_ik*gik;
		}		
			
		double bij=1./sqrt(1.+chi_ij);
					
		vi+=0.5*fc_ij*(vr_ij-bij*va_ij);
		
	}
	return vi;
}	
	

	
	
	
/************************************************************************************
 *	scaling function: fc(r)
 ************************************************************************************/
inline double compute_fc(const double & r,const double & R,const double & D){
	double fc;
	if(r<=fabs(R-D)) fc=1.;
	else if(r<=(R+D))fc=0.5*(1-sin(M_PI/2/D*(r-R)));
	else fc=0.;
	return fc;
}
/************************************************************************************
 *	gradient of scaling function: grad_ri(fc(rij)) - with respect to ri
 ************************************************************************************/
inline void compute_grad_fc_ri(const double & r,const double & R,const double & D,
					 const vector<double> & r_vec,  vector<double> & grad)
{
	double wsp;
	if(r<=fabs(R-D)) wsp=0.;
	else if(r<=(R+D))wsp=-M_PI/4/D*cos(M_PI/2/D*(r-R));
	else wsp=0.;
	for(int i=0;i<3;i++) grad[i]=wsp*(-r_vec[i]/r);
	return;
}




/************************************************************************************
 *	attractive potential: VA(r)
 ************************************************************************************/
inline double compute_va(const double & r,const double & D0, const double & S,const double & beta, const double & r0){
	double va=S*D0/(S-1.)*exp(-beta*sqrt(2./S)*(r-r0));
	return va;
}
/************************************************************************************
 *	gradient of attractive potential: grad(VA(r)) with respect to ri
 ************************************************************************************/
inline void compute_grad_va_ri(const double & r,const double & D0, const double & S,const double & beta, const double & r0,
				    const vector<double> & r_vec,  vector<double> & grad)
{					
	double wsp=-beta*sqrt(2./S)*compute_va(r,D0,S,beta,r0);
	for(int i=0;i<3;i++) grad[i]=wsp*(-r_vec[i]/r);
	return;
}

/************************************************************************************
 *	repulsive potential: VR(r)
 ************************************************************************************/
inline double compute_vr(const double & r,const double & D0, const double & S,const double & beta, const double & r0){
	double vr=D0/(S-1.)*exp(-beta*sqrt(2.*S)*(r-r0));
	return vr;
}	
/************************************************************************************
 *	gradient of repulsive potential: grad(VR(r)) with respect to ri
 ************************************************************************************/
inline void compute_grad_vr_ri(const double & r,const double & D0, const double & S,const double & beta, const double & r0,
				    const vector<double> & r_vec,  vector<double> & grad)
{
	double wsp=-beta*sqrt(2.*S)*compute_vr(r,D0,S,beta,r0);
	for(int i=0;i<3;i++) grad[i]=wsp*(-r_vec[i]/r);
	return;
}



/************************************************************************************
 *	cos(teta_ijk):
 ************************************************************************************/
inline double compute_cos_teta_ijk(const vector<double> & rij_vec,const double & rij,
					     const vector<double> & rik_vec,const double & rik)
{
	double cos_teta_ijk=0;
	for(int i=0;i<3;i++) cos_teta_ijk+=rij_vec[i]*rik_vec[i];
	return cos_teta_ijk/rij/rik;
}

/************************************************************************************
 *	g_ik:
 ************************************************************************************/
inline double compute_gik(const double & gamma,const double & c,const double & d,const double & h,const double & cos_teta_ijk){
	
	double c2=c*c;
	double d2=d*d;
	double hcos=h+cos_teta_ijk;
	double hcos2=hcos*hcos;
	double gik=gamma*(1+c2/d2-c2/(d2+hcos2));
	return gik;
}	
	
/************************************************************************************
 *	gradients(g_ik): with respect to xi (grad_i), xj (grad_j), xk(grad_k)
 ************************************************************************************/
void compute_grad_gik(const double & gamma,const double & c,const double & d,const double & h,const double & cos_teta_ijk,
			    const vector<double> & rij_vec,const double & rij,
			    const vector<double> & rik_vec,const double & rik,
			    vector<double> & grad_i,
			    vector<double> & grad_j,
			    vector<double> & grad_k)
{
	double c2=c*c;
	double d2=d*d;
	double hcos=h+cos_teta_ijk;
	double hcos2=hcos*hcos;
	double denominator=pow(d2+hcos2,2);
	double wsp=2*gamma*c2*hcos/denominator;
	
	add_vectors(3, wsp/rij/rik, rik_vec, -wsp*cos_teta_ijk/rij/rij, rij_vec, 0., grad_j);
	add_vectors(3, wsp/rij/rik, rij_vec, -wsp*cos_teta_ijk/rik/rik, rik_vec, 0., grad_k);
	add_vectors(3, -1., grad_j, -1., grad_k, 0., grad_i);
	
	return;
}	
	
	
	
	
/*********************************************************************************************************
 * 1: dla czastki o indeksie "i0" szukamy jej sasiadow (indeksow)
 * 2: czastki sa w tej samej komorce lub sasiednich
 * 3: stosujemy periodyczne WB w kierunkach x/y
 * n_neigh - liczba sasiednich czastek
 * neigh_index - tablica z indeksami sasiadow (zwracana)
 * neigh_r_vec - tablica z polozeniem sasiadow z uwzglednieniem periodycznych WB w kierunkach x/y
**********************************************************************************************************/	
	void find_nearest_particles_and_positions(const int & i0, int & n_neigh, vector<int> & neigh_index,vector<vector<double>> & neigh_r_vec){
		
		//numer komorki w ktorej siedzi czastka
		int kx=particle[i0].cell[0]; 
		int ky=particle[i0].cell[1];
		int kz=particle[i0].cell[2];
				
		n_neigh=0;
		for(int i=-1;i<=1;i++){   		//periodyczne WB
			for(int j=-1;j<=1;j++){		//periodyczne WB
				//sztywne WB
				int kmin,kmax;
				kz==0?kmin=0:kmin=-1;		//jesli jestemy w najnizszej komorce to NIE schodzimy nizej
				kz==(mz_cell-1)?kmax=0:kmax=1;//jesli jestemy w najwyzszej komorce to NIE wchodzimy wyzej
				for(int k=kmin;k<=kmax;k++){
					//indeks komorki z lokalnego otoczenia czastki
					int kx2=(kx+i+mx_cell) % mx_cell;	
					int ky2=(ky+j+my_cell) % my_cell;
					int kz2=kz+k;
					
					//szukamy czastek w komorce: linked-cell-method
					int ipos=header[kx2][ky2][kz2];  //czastka o najwyzszym indeksie w komorce
					while(ipos>0){
						if(ipos!=i0){
							neigh_index[n_neigh]=ipos;
							
							//polozenie czastki przesuwamy w zaleznosci czy komorka wychodzi poza obszar pudla (periodyczne WB)
							double x_shift=0.;
							double y_shift=0.;
							if((kx+i)==-1) 	x_shift=-xmax;  //czastke z prawej skrajnej komorki przesuwamy w lewo 
							else if((kx+i)==mx_cell)  x_shift= xmax; //czastke z lewej skrajnej komorki przesuwamy w prawo
							if((ky+j)==-1)	 y_shift=-ymax;  //czastke z prawej skrajnej komorki przesuwamy w lewo 
							else if((ky+j)==my_cell)  y_shift= ymax; //czastke z lewej skrajnej komorki przesuwamy w prawo
							neigh_r_vec[n_neigh][0]=particle[ipos].r_vec[0]+x_shift;
							neigh_r_vec[n_neigh][1]=particle[ipos].r_vec[1]+y_shift;
							neigh_r_vec[n_neigh][2]=particle[ipos].r_vec[2];
							
							n_neigh++;
						}
						ipos=link[ipos]; //siegamy do nastepnej pozycji w komorce
					}
				}
			}
		}
		return;	
	}
	
	
	
	
	
	
/*********************************************************************************************************
 * sztywne WB
**********************************************************************************************************/	
	void find_nearest_particles_and_positions_SZTYWNE_nieuzywane(const int & i0, int & n_neigh, vector<int> & neigh_index,vector<vector<double>> & neigh_r_vec){
		
		//numer komorki w ktorej siedzi czastka
		int kx=particle[i0].cell[0]; 
		int ky=particle[i0].cell[1];
		int kz=particle[i0].cell[2];
				
		
		
		
		n_neigh=0;
		
		//sztywne WB
		int imin,imax;
		kx==0?imin=0:imin=-1;		//jesli jestemy w najnizszej komorce to NIE schodzimy nizej
		kx==(mx_cell-1)?imax=0:imax=1;//jesli jestemy w najwyzszej komorce to NIE wchodzimy wyzej
		
		for(int i=imin;i<=imax;i++){   		//sztywne
			
			//sztywne WB	
			int jmin,jmax;
			ky==0?jmin=0:jmin=-1;		//jesli jestemy w najnizszej komorce to NIE schodzimy nizej
			ky==(my_cell-1)?jmax=0:jmax=1;//jesli jestemy w najwyzszej komorce to NIE wchodzimy wyzej
			
			
			for(int j=jmin;j<=jmax;j++){		//sztywne
				
				//sztywne WB
				int kmin,kmax;
				kz==0?kmin=0:kmin=-1;		//jesli jestemy w najnizszej komorce to NIE schodzimy nizej
				kz==(mz_cell-1)?kmax=0:kmax=1;//jesli jestemy w najwyzszej komorce to NIE wchodzimy wyzej
				
				for(int k=kmin;k<=kmax;k++){
					//indeks komorki z lokalnego otoczenia czastki
					int kx2=(kx+i+mx_cell) % mx_cell;	
					int ky2=(ky+j+my_cell) % my_cell;
					int kz2=kz+k;
					
					//szukamy czastek w komorce: linked-cell-method
					int ipos=header[kx2][ky2][kz2];  //czastka o najwyzszym indeksie w komorce
					while(ipos>0){
						if(ipos!=i0){
							neigh_index[n_neigh]=ipos;
							
							//polozenie czastki przesuwamy w zaleznosci czy komorka wychodzi poza obszar pudla (periodyczne WB)
							double x_shift=0.;
							double y_shift=0.;
							if((kx+i)==-1) 	x_shift=-xmax;  //czastke z prawej skrajnej komorki przesuwamy w lewo 
							else if((kx+i)==mx_cell)  x_shift= xmax; //czastke z lewej skrajnej komorki przesuwamy w prawo
							if((ky+j)==-1)	 y_shift=-ymax;  //czastke z prawej skrajnej komorki przesuwamy w lewo 
							else if((ky+j)==my_cell)  y_shift= ymax; //czastke z lewej skrajnej komorki przesuwamy w prawo
							neigh_r_vec[n_neigh][0]=particle[ipos].r_vec[0]+x_shift;
							neigh_r_vec[n_neigh][1]=particle[ipos].r_vec[1]+y_shift;
							neigh_r_vec[n_neigh][2]=particle[ipos].r_vec[2];
							
							n_neigh++;
						}
						ipos=link[ipos]; //siegamy do nastepnej pozycji w komorce
					}
				}
			}
		}
		return;	
	}	
	
	
	
	
/**************************************************************
 * sortowanie: linked-cell-method
***************************************************************/	
	void sort(){
		//flaga=-1 :blokada (na poczatku brak czastek w komorce)
		for(int i=0;i<mx_cell;i++){
			for(int j=0;j<my_cell;j++){
				for(int k=0;k<mz_cell;k++){
					header[i][j][k]=-1;
				}	
			}	
		}
		for(int i=0;i<N;i++) link[i]=-1;
		
		for(int i=0;i<N;i++){
			double x=particle[i].r_vec[0];
			double y=particle[i].r_vec[1];
			double z=particle[i].r_vec[2];
			int kx=min(int(x/dx_cell),mx_cell-1);
			int ky=min(int(y/dy_cell),my_cell-1);
			int kz=min(int(z/dz_cell),mz_cell-1);
			link[i]=header[kx][ky][kz];
			header[kx][ky][kz]=i;
			
			particle[i].cell[0]=kx;
			particle[i].cell[1]=ky;
			particle[i].cell[2]=kz;
		}		
		return;
	}
	
	
	
/************************************************************************************
 *	vector norm
 ************************************************************************************/
	inline double vector_norm(const int & n, const vector<double> & vec)
	{
		double cn=0.;
		for(int i=0;i<n;i++)cn+=vec[i]*vec[i];
		return sqrt(cn);
	}

/************************************************************************************
 *	add two scaled vectors
 ************************************************************************************/
	inline void add_vectors(const int & n, 
					const double & wsp1, const vector<double> & vec1,
					const double & wsp2, const vector<double> & vec2, 
					const double & wsp3, vector<double> & vec3)
	{
		for(int i=0;i<n;i++) vec3[i]=wsp1*vec1[i]+wsp2*vec2[i]+wsp3*vec3[i];
		return;
	}

/************************************************************************************
 *	zeroing vector components
 ************************************************************************************/
	inline void vector_zero(const int & n,vector<double> & vec){
		for(int i=0;i<n;i++) vec[i]=0.;
		return;
	}	

/************************************************************************************
 *	scale vector
 ************************************************************************************/
	inline void scale_vector(const int & n,const double & wsp, vector<double> & vec){
		for(int i=0;i<n;i++) vec[i]=wsp*vec[i];
		return;
	}
	
	
/**************************************************************
 * generator liczb pseudolosowych o rozkladzie jednorodnym
***************************************************************/	
	double uniform(){
		return (rand()+1.0)/(RAND_MAX+2.0);
	}
	
};
