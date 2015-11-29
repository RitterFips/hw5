#include <cmath>
#include <iostream>
#include <cstdlib>
#include <string>
#include <fstream>

using namespace std;
//Unterfunktionen
void write(double* const, double* const, const int, const double, const string, const string);
void forward(double* const, const int, const double);
void backward(double* const, const int, const double);

//Hauptfunktion
int main(){
  const int n = 10;
  const double tmax = 20*M_PI; // obere integrationsgrenze
  const double dt = M_PI/n;// Schrittgröße
  const int N = tmax/dt; // Anzahl der Schritte
  double* Mf = new double[2*N]; // Array um die Werte der forward-Methode zu speichern
  double* Mb = new double[2*N]; // Array um die Werte der backward-Methode zu speichern
  const string filename1 = "forward.txt";//Name der Zieldatei
  const string filename2 = "backward.txt";//Name der Zieldatei
  Mf[0] = 1; //Randbedingungen
  Mf[N] = 0;
  Mb[0] = 1;
  Mb[N] = 0;

  //Unterfunktionen aufrufen
  forward(Mf,N,dt);
  backward(Mb,N,dt);
  write(Mf, Mb, N, dt, filename1, filename2);
  //dynamische Arrays löschen
  delete[] Mf;
  delete[] Mb;
  return 0;
}

//Unterfunktionen
void forward(double* const Mf, const int N, const double dt){//forward-Methode
  for(int i = 1; i < N; i++){
    Mf[i] = Mf[i-1] + dt*Mf[N+i-1];
    Mf[N+i] = Mf[N+i-1] - dt*Mf[i-1];
  }
}

void backward(double* const Mb, const int N, const double dt){//Backward-Methode
  for(int i = 1; i < N; i++){
    Mb[i] = (Mb[i-1] + dt*Mb[N+i-1])/(1+dt*dt);
    Mb[N+i] = (Mb[N+i-1] - dt*Mb[i-1])/(1+dt*dt);
  }
}

void write(double* const Mf, double* const Mb, const int N, const double dt, const string fname, const string fname2){//Werte in  Files schreiben
 ofstream forw(fname.c_str());
 for(int i = 0; i < N; i++){
  forw << i*dt << "\t" << Mf[i] << "\t" << Mf[N+i] << endl; 
 }
 forw.close();
 ofstream out(fname2.c_str());
 for(int i = 0; i < N; i++){
  out << i*dt << "\t" << Mb[i] << "\t" << Mb[N+i] << endl; 
 }
 out.close();
}