// fft séquentielle 1D
// Implementation d'une transformee de Fourier rapide
#include <iostream>
#include <complex>
#include <vector>
#include <bitset>
#include <cstdlib>
#include <ctime>

using namespace std;
const double pi      = 3.14159265358979323846264;

// On réordonne les éléments de data en renversant
// l'ordre des bits des indices 
void bitInv(vector<complex<double> >& data) { /*à compléter*/ }

// Une etape de la FFT sequentielle
void stepSeq(vector<complex<double> >& data,
             complex<double>& w,int d) { /*à compléter*/ }

// FFT séquentielle
void fft(vector<complex<double> >& data) {
   bitInv(data);
   complex<double> w = polar(1.0,-pi);
   int n = data.size();
   for (int d=1;d<n;d*=2) stepSeq(data,w,d);
}

// Initialisation aléatoire d'un vecteur de nombres complexes
void randInit(vector<complex<double> >& data,double min, double max) {
   for (int k=0;k<data.size();k++)
      data[k] = complex<double>(min+(max-min)*drand48(),min+(max-min)*drand48());
}

// Affichage d'un vecteur de nombres complexes
void print(vector<complex<double> > data) {
   for (int k=0;k<data.size();k++)
      cout << data[k].real() << " " << data[k].imag() << endl;
}

int main(int argc,char ** argv) {
   int seed = atoi(argv[1]);
   srand48(seed);
   int n = atoi(argv[2]);
   vector<complex<double> > data(n);
   randInit(data,0.0,100.0);
   print(data);
   fft(data);
   print(data);
}
