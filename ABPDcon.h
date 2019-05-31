
//N- number of state variables
#define N 11

double time_;
double step;
double state[N];
double deriv[N];
double   arg[N];

#define E_MAX         0.0000001      /* max voltage change (mV/step) */
#define E_MIN         0.000000031    /* min voltage change (mV/step) */
#define STEP          0.02           /* init step size (ms) */

#define START_TIME    0    /* start time in ms */
#define END_TIME    13000/* end time in ms */


//capacitance values
#define CMa 1.0  //nF
#define CMs 1.0 //nF
#define CMd 1.0  //nF

double isoma,idend;
double iext=0.2; //nA0.2
double iext_s=0.0; 
double iext_d=0.0; 


//conductance values
#define gna 15.0//15.0 microsiemens
#define gk 8.0 //8.0 microsiemens
#define gca 0.04//0.04 microsiemens

#define gl 0.0354 //0.0854 microsiemens
#define gl_pd 0.001
#define gl_s 0.001 // microsiemens
#define gl_a 0.001 // microsiemens

#define G_PS 0.065//0.065 microsiemens
#define G_kf 0.07//0.07 microsiemens

double ga = 100;//100.00 microsiemens
double gkca = 0.273537;//0.273537 microsiemens
double gsyn; //microsiemens or 1000 nS
double Gsa=0.0,Gdpd=0.04, Gapd=0.5, Gda=0.0,Gspd=0.05;

//reversal potential
#define vk    -75.0    /* membrane currents. mV */
#define vl    -40.0 //-40.0    /* V_x= x_NERNST_POT 40.0 mV */
#define vna   30.0     /* V_x= x_NERNST_POT mV */
#define vca   140.0 //mV
#define vsyn 0.0 //mV

//values for kinetic variables

#define rho 0.0016
#define lambdan 0.8
#define lambdah 0.8
#define ka 1.0
#define tauz 23
#define kca 0.0078
#define btau 1
#define zb -50.0
#define va -12.0
#define vb -62.0
#define sa -26.0
#define sb 6.0

#define tps 3000
#define tpsmin 100//20
#define psst 0.05
double slope_ca_act=12,ca_inact_var=0.43;

double tm_ica,tm_ikca;
FILE *t_ica,*t_ikca;
double i_ds,i_as;
FILE *rate_ca;
