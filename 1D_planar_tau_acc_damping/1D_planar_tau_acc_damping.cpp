#include "omp.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cmath>
#include <sys/time.h>
using namespace std;

#define pi 3.1415926
#define yr 31536000.0
#define Msun 1.98847e33
#define me 9.1e-28
#define mp 1.673e-24
#define mH 1.674e-24
#define c 3.0e10
#define e 4.8e-10
#define mec2 8.19e-7
#define mpc2 1.5e-3
#define r0 e * e / mec2
#define kB 1.38e-16
#define sigmaT 6.65e-25
#define pc 3.08e18
#define h 6.626e-27
#define hbar h / 2.0 / pi
#define eV 1.6e-12 // erg
#define GeV 1.6e-3
#define TeV 1.6
#define PeV 1.6e3

int main(){

    double Tage = 3.0e6*yr;
    double z0 = 96.0*pc;
    double u1 = 2.8e8;
    double u2 = u1/4.0;
    double Mdot = 1.5e-4*Msun/yr;
    double etaB = 1.0;

    double Lw = 0.5*Mdot*pow(u1,2.0);
    double V = 4.0*pi/3.0*pow(z0,3.0);
    double S = etaB*Lw/V;
    double rho = Mdot*Tage/V;
    double a = 0.3/sqrt(4.0*pi);

    double Lc = 1.0*pc;
    double k0 = 2.0*pi/Lc;
    double B = sqrt(4.0*pi*3.0/2.0*pow(sqrt(rho)*S/a/k0,2.0/3.0));
    double vA = B/sqrt(4.0*pi*rho);
    double Sqrt_Delta = sqrt(pow(u2,2.0)-4.0/9.0*pow(vA,2.0));

    cout << "density: " << rho/mp << endl;
    cout << "magnetic field: " << B << endl;
    cout << "Alven speed: " << vA << endl;
    cout << "velocity ratio: " << u2/vA << ", " << u2/Sqrt_Delta << endl;

    // create meshgrid in momentum space

    int nj = 7001;
    int i, j, k;
    double p[nj];
    for(j=0; j<nj; ++j){
        p[j] = pow(10, j*7.0/(nj-1))*GeV/c;
    }

    // calculate the diffusion coefficient

    double Dz_i[nj];
    double rL, kres, Wp;
    for(j=0; j<nj; ++j){
        rL = p[j]*c/e/B;
        kres = 2.0*pi/rL;
        Wp = pow(kres,-5.0/3.0)*pow(sqrt(rho)*S/a,2.0/3.0);
        Dz_i[j] = 1.0/3.0*c*rL/(kres*Wp/(B*B/4.0/pi));
    }

    cout << "diffusion coefficient: " << Dz_i[0] << ", " << Dz_i[3000] << ", " << Dz_i[6000] << endl;

    // timescale ratio

    double tau_sa[nj], tau_esc[nj];
    for(j=0; j<nj; ++j){
        tau_sa[j] = 9.0*Dz_i[j]/pow(vA,2.0);
        tau_esc[j] = pow(z0,2.0)/Dz_i[j];
    }

    cout << "timescale ratio: " << tau_sa[6000]/tau_esc[6000] << endl;

    // calculate the spectrum at the shock

    double Int_1, Int_2;

    double f0_1[nj], f0_2[nj];
    for(j=0; j<nj; ++j){
        Int_1 = 0.0;
        Int_2 = 0.0;
        for(i=0; i<j; ++i){
            Int_1 += -3.0*u2/(u2-u1)/(1.0-exp(u2*z0/Dz_i[i]))/p[i]*(p[i+1]-p[i]);
            Int_2 += -(3.0/2.0*u2/(u2-u1)+3.0*Sqrt_Delta/(u2-u1)*(1.0/(1.0-exp(Sqrt_Delta*z0/Dz_i[i]))-1.0/2.0))/p[i]*(p[i+1]-p[i]);
        }
        f0_1[j] = pow(p[j],3.0*u1/(u2-u1))*exp(Int_1);
        f0_2[j] = pow(p[j],3.0*u1/(u2-u1))*exp(Int_2);
    }

    // calculate normalization
 
    // double Pram = rho*u1*u1;
    double Pram = Mdot*u1/(4.0*pi*pow(15.0*pc,2.0));
    double v[nj]; 
    for(j=0; j<nj; ++j){
        v[j] = p[j]*c/sqrt(pow(mp*c,2.0)+pow(p[j],2.0));
    }

    double Ppar_1 = 0.0, Ppar_2 = 0.0;
    for(j=0; j<nj-1; ++j){ 
        Ppar_1 += 1.0/3.0 * 4.0*pi*pow(p[j], 2.0) *  f0_1[j] * p[j] * v[j] * (p[j+1]-p[j]);
        Ppar_2 += 1.0/3.0 * 4.0*pi*pow(p[j], 2.0) *  f0_2[j] * p[j] * v[j] * (p[j+1]-p[j]);
    }

    double xi = 0.1;
    double norm_1 = xi * Pram / Ppar_1;
    double norm_2 = xi * Pram / Ppar_2;

    for(j=0; j<nj; ++j){
        f0_1[j] = norm_1 * f0_1[j];
        f0_2[j] = norm_2 * f0_2[j];
    }

    // calculate spectrum away from the shock

    double z = 0.5*z0;

    double B1[nj], B2[nj];
    for(j=0; j<nj; ++j){
        B1[j] = (u2+Sqrt_Delta)/2/Dz_i[j];
        B2[j] = (u2-Sqrt_Delta)/2/Dz_i[j];
    }

    double f_1[nj], f_2[nj];
    for(j=0; j<nj; ++j){
        f_1[j] = f0_1[j]/(1.0-exp(-u2*z0/Dz_i[j]))*(1.0-exp(u2*(z-z0)/Dz_i[j]));
        f_2[j] = f0_2[j]/(1.0-exp((B1[j]-B2[j])*z0))*(exp(B1[j]*z)-exp((B1[j]-B2[j])*z0)*exp(B2[j]*z));
    }

    // update the diffusion coefficient

    double Dz_f[nj], Int_3[nj], Gamma, dk, kp;
    #pragma omp parallel for private(i, k, rL, kres, dk, kp, Gamma, Wp) shared(Dz_f, f0_2, p) num_threads(30)
    for(j=0; j<nj; ++j){
        rL = p[j]*c/e/B;
        kres = 2.0*pi/rL;
        Int_3[j] = 0.0;
        if(kres>k0){
            dk = (kres-k0)/1000.0;
            kp = k0;
            for(k=0; k<1000; ++k){
                Gamma = 0.0;
                for(i=0; i<nj; ++i){
                    if(p[i]>e*B/kp/c){
                        Gamma += f0_2[i]*p[i]*(p[i+1]-p[i]);
                    }
                }
                Gamma = 8.0*pow(pi,3.0)*pow(e,2.0)*pow(vA/c,2.0)/kp * Gamma;
                Int_3[j] += Gamma/pow(kp,5.0/3.0)*dk;
                kp = kp + dk;
            }
            Int_3[j] = 1.0/3.0*pow(rho/S/a/a,1.0/3.0) * Int_3[j];
        }
        Wp = pow(kres,-5.0/3.0)*pow(sqrt(rho)*S/a,2.0/3.0) * pow(1.0-Int_3[j],2.0);
        Dz_f[j] = 1.0/3.0*c*rL/(kres*Wp/(B*B/4.0/pi));
    }

    cout << "updated diffusion coefficient: " << Dz_f[0] << ", " << Dz_f[3000] << ", " << Dz_f[6000] << endl;

    // output to file

    ofstream fout;

    fout.open("p.csv");
    for(j=0; j<nj; ++j)
        fout << p[j] << endl;
    fout.close();

    fout.open("Dz_i.csv");
    for(j=0; j<nj; ++j)
        fout << Dz_i[j] << endl;
    fout.close();

    fout.open("Dz_f.csv");
    for(j=0; j<nj; ++j)
        fout << Dz_f[j] << endl;
    fout.close();

    fout.open("tau_sa.csv");
    for(j=0; j<nj; ++j)
        fout << tau_sa[j] << endl;
    fout.close();

    fout.open("tau_esc.csv");
    for(j=0; j<nj; ++j)
        fout << tau_esc[j] << endl;
    fout.close();

    fout.open("f0_1.csv");
    for(j=0; j<nj; ++j)
        fout << f0_1[j] << endl;
    fout.close();

    fout.open("f0_2.csv");
    for(j=0; j<nj; ++j)
        fout << f0_2[j] << endl;
    fout.close();

    fout.open("f_1.csv");
    for(j=0; j<nj; ++j)
        fout << f_1[j] << endl;
    fout.close();

    fout.open("f_2.csv");
    for(j=0; j<nj; ++j)
        fout << f_2[j] << endl;
    fout.close();

}