#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>

/*Spins states SPIN UP or SPIN DOWN*/
#define SPIN_UP 1 
#define SPIN_DOWN -1

/*Lattice side size. Lattice size is NxN*/
#define N 128

/*Coordination sphere number for square lattice  equal 4, for simple cubic lattice equal 6 */
#define Z  4

/*Interaction spin parameter. For model paramagnetic = 0; fero- = 1; antifero- = -1;*/
const double Jij = 1;

/*Displacement array*/
const int di[Z] = {0,0,1,-1};
const int dj[Z] = {1,-1,0,0};

/* Spins array*/
 int s[N][N] = {0}; 

/* Periodic boundary conditions */
int bk(int index)
{
	int bc = index;
	if (index < 0) bc = index + N;
	if (index > N-1) bc = index - N;
	
	return bc;	
}

/*Randomize function (generation seed for pseudo random) */
void randomize(void) {
	srand(time(NULL));
}

/*Random function  [0; 1) */
double rnd(void) {
	return  ((double) rand()) / ((double) RAND_MAX+1);
}

/*Output spins array to xyz format for Ovito */
void out_to_xy(const char *filename);

/*Random state initilization*/
void init() {
	randomize();
	for (int i = 0; i<N; i++)
		for(int j = 0; j<N; j++)
		{
			double r = rnd();
			if (r<0.5) 
			  	s[i][j] = SPIN_DOWN;
			else
				s[i][j] = SPIN_UP;
		}

}

/* average magnetization per  spin */
double get_magnetization(void) {
	double M = 0;
	for (unsigned i = 0; i < N; i++)
		 for (unsigned j =0; j<N; j++)
		  M += s[i][j];

 return M/(N*N);		
}


/*Hamiltonian function with external magnetic field*/
double hamiltonian(int i, int j,  double B){
	int sum = 0;
	for (int nn = 0; nn<Z; nn++)  sum += s[bk(di[nn]+i)][ bk(dj[nn]+j)];
	return -s[i][j]*(Jij*sum + B);	
}

/*Ising model fuction that let us output to .xyz files result of modelling*/
void model_with_plot(int mc_steps,  double T, double B,  int nOut){
  out_to_xy("0.xyz");
  randomize();
  int ns = 1;
  for (int l = 1; l<=mc_steps*N*N; l++)
  {
  	 int i = (int)floor(rnd()*N);
  	 int j = (int)floor(rnd()*N);

  	 double dE = -2*hamiltonian(i,j, B);
  	 if (dE < 0) s[i][j] = -s[i][j];
  	 else 
  	 	if (rnd() < exp(-dE/T) )  s[i][j] = -s[i][j];

  	 if (l % (nOut*N*N) == 0)
  	 {
  	 	ns+=1;
  	 	char fname[255];
  	 	sprintf(fname, "%d.xyz",ns);

  	 	out_to_xy(fname);
  	 }
  	 
  } 
}

 /*Ising model that calculate  the average  magnetization per spin (<M/(N*N)> ) */
double model_avg_magnetization(int mc_steps,  double T, double B){
	double M = 0;
	for (int l = 1; l<=mc_steps*N*N; l++)
  {
  	 int i = (int)floor(rnd()*N);
  	 int j = (int)floor(rnd()*N);

  	 double dE = -2*hamiltonian(i,j, B);
  	 if (dE < 0) s[i][j] = -s[i][j];
  	 else 
  	 	if (rnd() < exp(-dE/T) )  s[i][j] = -s[i][j];
  	 M += get_magnetization();
  } 
  return M / (N*N*mc_steps);
}

/*Function that let us get the hysteresis loop of average magnetization per spin  <M/(N*N)> = F(B)*/
void check_hysteresis(double T, double Bmax, int mc_steps) {
	printf("Checking Hysteresis\n");
	char fname[255];
	sprintf(fname,"hysteresis_T_%.2f.txt",T);
	FILE *file = fopen(fname, "w");
	double B = 0;
	double M = model_avg_magnetization(mc_steps, T, B);
	printf("Step 1 \nfrom 0 to Bmax\n");
	while (B<Bmax){
		M = model_avg_magnetization(mc_steps, T, B);
		//fprintf(file, "%g   %g\n",B, M);
		B+= 0.01;	
	}

	printf("Step 2 \nfrom Bmax to -Bmax\n");
	while (B>-Bmax){
		M = model_avg_magnetization(mc_steps, T, B);
		fprintf(file, "%g   %g\n",B, M);
		B-= 0.01;	
	}

	printf("Step 3 \nfrom -Bmax to Bmax\n");
	while (B<Bmax){
		M = model_avg_magnetization(mc_steps, T, B);
		fprintf(file, "%g   %g\n",B, M);
		B+= 0.01;	
	}

	fclose(file);
}

/*Calculating average energy per spin */
double energy_per_spin(void){

	double sum = 0;
	for(unsigned i = 0; i<N; i++)
		for (unsigned j = 0; j<N; j++)
		sum += hamiltonian(i,j,0);
	return sum / (N*N)/2.0;
}

/*Ising model with calculation average energy of system */
double model_energy(int mc_steps,  double T, double B){
	double Etotal = 0;
  for (int l = 1; l<=mc_steps*N*N; l++)
  {
  	 int i = (int)floor(rnd()*N);
  	 int j = (int)floor(rnd()*N);

  	 double dE = -2*hamiltonian(i,j, B);
  	 if (dE < 0) s[i][j] = -s[i][j];
  	 else 
  	 	if (rnd() < exp(-dE/T) )  s[i][j] = -s[i][j];

  	 Etotal += energy_per_spin();

  } 
  return Etotal / (mc_steps*N*N);
}

/*Energy per spin dependence from temperature */
double energy_dep_temperature(void){
	int mc_steps = 10000;
	double T = 0.1, avE = 0;
	init();
	char fname[255];
	sprintf(fname,"%dx%d_E_T.txt",N,N);
	FILE *file = fopen(fname,"w");
	avE = model_energy(mc_steps, T, 0);
	while (T<=4) {

	  avE = model_energy(mc_steps, T, 0);
	  fprintf(file, "%g   %g\n",T, avE);

	  T += 0.1;
	}
	fclose(file);
}

 /*Ising model that calculate  the average  absolute magnetization per spin (<|M|/(N*N)> ) */
double model_abs_avg_magnetization(int mc_steps,  double T, double B){
	double M = 0;
	for (int l = 1; l<=mc_steps*N*N; l++)
  {
  	 int i = (int)floor(rnd()*N);
  	 int j = (int)floor(rnd()*N);

  	 double dE = -2*hamiltonian(i,j, B);
  	 if (dE < 0) s[i][j] = -s[i][j];
  	 else 
  	 	if (rnd() < exp(-dE/T) )  s[i][j] = -s[i][j];
  	 M += fabs(get_magnetization());
  } 
  return M / (N*N*mc_steps);
}

/* Absolute Average Magnetization per spin dependence from temperature */
void magnetization_dep_temperature(void){
	int mc_steps = 10000;
	double T = 0.1, avM = 0;
	init();
	char fname[255];
	sprintf(fname,"%dx%d_M_T.txt",N,N);
	FILE *file = fopen(fname,"w");
	avM = model_abs_avg_magnetization(mc_steps, T, 0);
	while (T<=4) {

	  avM = model_abs_avg_magnetization(mc_steps, T, 0);
	  fprintf(file, "%g   %g\n",T, avM);

	  T += 0.1;
	}
	fclose(file);
}

double sqr(double x){
	return x*x;
}

/* model (dE)^2 per spin calculation for checking C(T) - heat capasity as function of temperature*/
double model_de2(int mc_steps,  double T, double B){
	double Etotal = 0;
	double E2total =0;
  for (int l = 1; l<=mc_steps*N*N; l++)
  {
  	 int i = (int)floor(rnd()*N);
  	 int j = (int)floor(rnd()*N);

  	 double dE = -2*hamiltonian(i,j, B);
  	 if (dE < 0) s[i][j] = -s[i][j];
  	 else 
  	 	if (rnd() < exp(-dE/T) )  s[i][j] = -s[i][j];

  	 double E = energy_per_spin();	
  	 Etotal += E;
  	 E2total += sqr(E);

  } 
  return ((E2total/(mc_steps*N*N)) -  sqr(Etotal/(mc_steps*N*N)) )/sqr(T);
}

/* heat capasity dependence from temperature C(T) */
void heat_capasity_dep_temperature(void){
	double T = 0.1, dE2 = 0;
	int mc_steps = 10000;
	init();
	char fname[255];
	sprintf(fname,"%dx%d_C_T.txt",N,N);
	FILE *file = fopen(fname,"w");
	dE2 = model_de2(mc_steps, T, 0);
	while (T<=4) {

	  dE2 = model_de2(mc_steps, T, 0);
	  fprintf(file, "%g   %g\n",T, dE2);

	  T += 0.1;
	}
	fclose(file);
}

double model_susceptibility(int mc_steps,  double T, double B) {
	double Mtotal = 0,  M2total = 0;
	for (int l = 1; l<=mc_steps*N*N; l++)
	  {
	  	 int i = (int)floor(rnd()*N);
	  	 int j = (int)floor(rnd()*N);

	  	 double dE = -2*hamiltonian(i,j, B);
	  	 if (dE < 0) s[i][j] = -s[i][j];
	  	 else 
	  	 	if (rnd() < exp(-dE/T) )  s[i][j] = -s[i][j];

	  	 double M   = get_magnetization();	
	  	 M2total += sqr(M);
	  	 Mtotal += fabs(M);
	  } 
 return (M2total/(mc_steps*N*N) -  sqr(Mtotal/(mc_steps*N*N)))/T;	  
}


void susceptibility_dep_temperature(void){
	double T = 0.5, hi = 0;
	int mc_steps = 10000;
	init();
	char fname[255];
	sprintf(fname,"%dx%d_Sc_T.txt",N,N);
	FILE *file = fopen(fname,"w");
	hi = model_susceptibility(mc_steps, T, 0);
	while (T<=4) {

	  hi = model_susceptibility(mc_steps, T, 0);
	  fprintf(file, "%g   %g\n",T, hi);

	  T += 0.1;
	}
	fclose(file);
}

/*Ising model with calculation average energy of system */
void model_magnetization_dep_MC(int mc_steps,  double T, double B){
	double Mtotal = 0;
	char fname[255];
	sprintf(fname,"%dx%d_M_mcs_%.2f.txt",N,N,T);
	FILE *file = fopen(fname,"w");
  for (int l = 1; l<=mc_steps*N*N; l++)
  {
  	 int i = (int)floor(rnd()*N);
  	 int j = (int)floor(rnd()*N);

  	 double dE = -2*hamiltonian(i,j, B);
  	 if (dE < 0) s[i][j] = -s[i][j];
  	 else 
  	 	if (rnd() < exp(-dE/T) )  s[i][j] = -s[i][j];

  	 Mtotal += get_magnetization();
  	 if(l % (N*N) ==0 )
  	 {
  	 	 fprintf(file, "%d  %g\n",l/(N*N), Mtotal/(N*N));
  	 	 Mtotal = 0;
  	 }
  	
  } 
  fclose(file);
}

int main(int argc, char const *argv[])
{
	
	///energy_dep_temperature();
	//magnetization_dep_temperature();
	//heat_capasity_dep_temperature();
	//magnetization_dep_temperature();
	//susceptibility_dep_temperature();
	//energy_dep_temperature();
	init();
	/*for (double T =1; T<=3; T+= 0.5)
	{
		check_hysteresis(T, 2, 5000);
	}*/
	//model_magnetization_dep_MC(1000,2.25,0);
	model_with_plot(100, 4.0, 0, 1);
	return 0;
}

void out_to_xy(const char *filename) {
	FILE *file = fopen(filename, "w");
  	fprintf(file, "%d\n",N*N);
  fprintf(file, "X Y P\n");

  for (int i = 0; i < N; i++)
  	for(int j = 0; j < N; j++ )
  	{
  		fprintf(file, "%d %d %d\n",i,j, s[i][j]);
  	}
  	
  fclose(file);
}