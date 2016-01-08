// fft séquentielle 1D
// Implementation d'une transformee de Fourier rapide
#include <iostream>
#include <complex>
#include <vector>
#include <bitset>
#include <cstdlib>
#include <ctime>
#include <cmath>

using namespace std;
const double pi      = 3.14159265358979323846264;

// inversion des données en fonction des nouveaux indices
void bitInv(vector<complex<double> >& data) {
    complex<double> tmp;
    int resI;
    int size = log2(data.size());

    bool index[data.size()] = {false};

    // Parsing only the half array
    for (int i = 0; i < data.size(); i++) {
        if (!index[i]) {
            resI = 0;
            for (int j = 0; j < size; j++) {
                if (i & (0x1 << j)) {
                   resI |= 0x1 << (size - 1 - j);
                }
            }

            tmp = data[i];
            data[i] = data[resI];
            data[resI] = tmp;
            index[i] = true;
            index[resI] = true;
        }
    }
}



// Une etape de la FFT sequentielle
void stepSeq(vector<complex<double> >& data,
             complex<double>& w,int d) {
    complex<double> wk, impair;
    for (int k = 0; k < data.size(); k++) {
        if (k % (2*d) == 0)
            wk = 1;
        if ((k & (0x1 << (int)log2(d))) == 0) {
            impair = wk * data[k + d];
            data[k + d] = data[k] - impair;
            data[k] = data[k] + impair;
            wk = w * wk;
        }
    }
    w = sqrt(w);
}

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
void print(vector<complex<double> > data, char* label) {
<<<<<<< HEAD
    cout << label << "=[";
=======
    cout << label << " = [";
>>>>>>> 105867f708e9b65bb963c65d358a02fbf0f9e261
   for (int k=0;k<data.size();k++)
      cout << data[k].real() << "+(" << data[k].imag() << "i);" << endl;
   cout << "]" << endl;
}

int main(int argc,char ** argv) {
   int seed = atoi(argv[1]);
   srand48(seed);
   int n = atoi(argv[2]);
   vector<complex<double> > data(n);
   randInit(data,0.0,100.0);
   print(data, "A");
   fft(data);
<<<<<<< HEAD
   //cout << "After FFT" << endl;
   print(data, "B");
   //cout << "Done" << endl;
=======
   cout << "After FFT" << endl;
   print(data, "B");
   cout << "Done" << endl;
>>>>>>> 105867f708e9b65bb963c65d358a02fbf0f9e261
}
