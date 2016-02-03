#include <iostream>
#include <fstream>
#include <complex>
#include <sstream>
//-----------------------------------
using namespace std;
//-----------------------------------
typedef complex<double> cmplx;
//-----------------------------------
void init( cmplx* const psi0, const double alpha, const double lambda,
           const double dx, const double dt,
           const int Nx, const double xmin);


void writeToFile(const cmplx* const v, const string s, const double dx,
         const int Nx, const double xmin, const double alpha,
         const double lambda, const double omega, const double t);

void step(cmplx* const, cmplx* const, const double, const double,
	  const double, const double, const int, const cmplx, const double k, double* v);
//-----------------------------------
int main(){

	const int Nx = 300;
	const double xmin = -40;
  const double xmax = 40;
	const double Tend = 10*M_PI;
	const double dx = (xmax-xmin)/(Nx-1);
	const double dt =  dx/10 ;
	const cmplx di = cmplx(0,1);
  double t = 0;
	const int Na = 10;
	int Nk = int(Tend / Na / dt + 0.5);

	const double lambda = 10;
  const double omega = 0.2;
  
  cmplx* u1 = new cmplx[Nx];
  cmplx* h;
  const double alpha=pow(0.04,1/4);
  const double D=1;
  double v[Nx];
  const double k = omega*omega;

  stringstream strm;

	cmplx* psi0 = new cmplx[Nx];

	init(psi0, alpha, lambda, dx, dt, Nx, xmin);

	writeToFile(psi0,"psi_0", dx,Nx,xmin, alpha, lambda, omega,t);



	for (int i = 1; i <= Na; i++) {
		for (int j = 1; j <= Nk-1; j++) {

		  step(u1, psi0, dt, dx, xmin, D, Nx, di, k, v);
		  h=psi0;
		  psi0=u1;
		  u1=h;
         t+=dt;
		}
		strm.str("");
		strm << "psi_" << i;
		writeToFile(psi0,strm.str(), dx,Nx,xmin, alpha, lambda, omega,t);
	}
  cout << "t = " << t << endl;
	
	delete[] psi0;
	delete[] u1;
	return 0;
}
//-----------------------------------
void step(cmplx* const f1, cmplx* const f0, const double dt, const double dx,
	  const double xmin, const double D, const int Nx, const cmplx di, const double k, double* v)
{
cmplx* d=new cmplx[Nx];
cmplx* u=new cmplx[Nx];
cmplx* l=new cmplx[Nx];

  double x;
  for(int i=0;i<Nx;i++){
    x=xmin+i;
    v[i]=k*x*x/2;
  }
  for(int i=0;i<Nx;i++) d[i] = 1. + di*(dt/(2*dx*dx)+(dt/2)*v[i]);
  for(int i=0;i<Nx;i++) u[i] = - di*dt/(4*dx*dx);
  for(int i=0;i<Nx;i++) l[i] = - di*dt/(4*dx*dx);
  for(int i=1;i<Nx;i++){
    d[i]-=u[i-1]*l[i]/d[i-1];
    f0[i]-=f0[i-1]*l[i]/d[i-1];
    }
  f1[Nx-1]=f0[Nx-1]/d[Nx-1];
  for(int i=Nx-2;i>=0; i--) f1[i]=(f0[i]-u[i+1]*f1[i+1]);
delete[] d;
delete[] u;
delete[] l;
}
//-----------------------------------
void writeToFile(const cmplx* const v, const string s, const double dx,
                 const int Nx, const double xmin, const double alpha,
                 const double lambda, const double omega, const double t)
{
	ofstream out(s.c_str());
  double x, xi, xil;
  double h1, h2, h3;
  cmplx ana;
	for(int i=0; i<Nx; i++){
		x = xmin + i * dx;
    xi = alpha * x;
    xil = alpha * lambda;
    h1 = -0.5 * pow(xi - xil*cos(omega*t),2 );
    h2 = omega*t/2 + xi * xil * sin(omega*t);
    h3 =  - 0.25 * xil*xil* sin(2*omega*t);
    ana = cmplx( h1 , h2 + h3  );
    ana = sqrt(alpha / sqrt(M_PI) ) * exp(ana);
		out << x << "\t" << norm(v[i]) << "\t" << v[i].real() << "\t" << v[i].imag()
         << "\t" << norm(ana) << "\t" << ana.real() << "\t" << ana.imag() <<  endl;
	}
	out.close();
}
//-----------------------------------

void init( cmplx* const psi0, const double alpha, const double lambda,
           const double dx, const double dt,
           const int Nx, const double xmin)
{
	const double x0 = dx*Nx * 0.5;
	for(int i=0;i<Nx; i++){
		double x = xmin + i*dx ;
		psi0[i] = sqrt(alpha/sqrt(M_PI)) * exp(- pow(alpha*(x-lambda),2)/2 );
	}
}
