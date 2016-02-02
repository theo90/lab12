#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <complex>
#include <sstream>
#include <cmath>
//-----------------------------------
using namespace std;
//-----------------------------------
typedef complex<double> cmplx;
//-----------------------------------
void init( cmplx* const psi0, const double eta, const double sigma, const double dx,
          const int Nx);

void writeToFile(const cmplx* const v, const string s, const double dx,
                 const int Nx, const double xmin, const double t);
void lstep(cmplx* v0, cmplx* v1, const double dx, const double dt, const int Nx);
void rstep(cmplx* v0, const double dt, cmplx* v1, const int Nx);
//-----------------------------------
int main(){

	const int Nx = 4000;
	const double L = 800;
	const double xmin = 0;
	const double Tend = 50;
	const double dx = L / (Nx - 1);
	const double dt = dx  / 10;
	const int Na = 10;
	int Nk = int(Tend / Na / dt + 0.5);
	double t=0;
	const double eta = 0.2;

	stringstream strm;

	cmplx* psi0 = new cmplx[Nx];
	cmplx* psi1 = new cmplx[Nx];
	cmplx* h;
	init(psi0, eta, dx, dt,Nx);

	writeToFile(psi0,"psi_0.txt", dx,Nx,xmin,t);

	
	for (int i = 1; i <= Na; i++) {

		for (int j = 1; j <= Nk-1; j++) {
			lstep(psi0,psi1,dx,dt/2,Nx);
			rstep(psi1,dt,psi0,Nx);
			lstep(psi0,psi1,dx,dt/2,Nx);
			h = psi0;
			psi0 = psi1;
			psi1 = h;
			t +=dt;
		}
		
		strm.str("");
		strm << "psi_" << i<<".txt";
		writeToFile(psi0,strm.str(), dx,Nx,xmin,t);
	}

	return 0;
}
void lstep(cmplx* v0, cmplx* v1, const double dx, const double dt, const int Nx)
{
	cmplx* d=new cmplx[Nx];
	cmplx* u=new cmplx[Nx];
	cmplx* l=new cmplx[Nx];

	for(int i=0; i<Nx; i++)
	{
		d[i]=cmplx(1., -2*dt/(dx*dx));
		u[i]=cmplx(0.,dt/(dx*dx));
		l[i]=cmplx(0.,dt/(dx*dx));
	}
	for(int i=1;i<Nx;i++)
	{
		d[i]-=l[i]*u[i-1]/d[i-1];
		v0[i]-=l[i]/d[i-1]*v0[i-1];
	}
	v1[Nx-1]=v0[Nx-1]/d[Nx-1];
	for(int i=Nx-2;i>-1;i--)
		v1[i]=(v0[i]-u[i]*v1[i+1])/d[i];

	delete[] d;
	delete[] u;
	delete[] l;
}

void rstep(cmplx* v0, const double dt, cmplx* v1, const int Nx)
{
	for (int i=0; i<Nx; i++)
		v1[i]=v0[i]*exp(cmplx(0., -norm(v0[i])*dt));
}
//-----------------------------------
void writeToFile(const cmplx* const v, const string s, const double dx,
                 const int Nx, const double xmin, const double t)
{
	ofstream out(s.c_str());
	for(int i=0; i<Nx; i++){
		double x = xmin + i * dx + t*5;
		out << x << "\t" << norm(v[i]) << "\t" << v[i].real() << "\t" << v[i].imag() << endl;
	}
	out.close();
}
//-----------------------------------

void init( cmplx* const psi0, const double eta,  const double dx, const double dt,
          const int Nx)
{
	const double x0 = dx*Nx * 0.5;
	const double f = sqrt(2) * eta;
	for(int i=0;i<Nx; i++){
		double x = i*dx - x0;
		psi0[i] = 2*f/cosh(eta * x);
	}
}

