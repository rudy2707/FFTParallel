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

// On réordonne les éléments de data en renversant
// l'ordre des bits des indices
//void bitInv(vector<complex<double> >& data) {
//    complex<double> tmp;
//    int resI;
//    int size = log2(data.size());
//
//    // Parsing only the half array
//    for (int i = 1; i < data.size()/2; i++) {
//        resI = 0;
//        for (int j = 0; j < size; j++) {
//            if (i & (0x1 << j)) {
//               resI |= 0x1 << (size - 1 - j);
//            }
//        }
//        tmp = data[i];
//        data[i] = data[resI];
//        data[resI] = tmp;
//    }
//}


int inversionBits(int nbBit, int indice){
	int m = 0;
	int i = indice;
	for(int b = 0; b < nbBit; b++){
		m = m|(i&1);
		m = m<<1;
		i = i>>1;
	}
	m = m >> 1;
	return m;
}

bool IndiceDansTableau(int* tab, int indice, int taille){
	for(int i = 0; i < taille; i++){
		if(tab[i] == indice)
			return true;
	}
	return false;
}

// inversion des données en fonction des nouveaux indices
void bitInv(vector<complex<double> >& data) {
	complex<double> tmp;
	int n = data.size();
	int* indices = new int[n];
	for(int i = 0; i < n; i++)
		indices[i] = i;

	for(int k = 0; k < n; k++){
		int m = inversionBits(log2(n),k);
		if(IndiceDansTableau(indices,m,n) ) { // test si l'échange a déjà été fait
			//cout << "k " << k << " i " << m << endl;
			//cout << "dataM " << data[m].real() << " dataK " << data[k].real() << endl;
			tmp = data[m];
			data[m] = data[k];
			data[k] = tmp;
			indices[m] = -1;
			indices[k] = -1;
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
    cout << label << "=[";
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
   //cout << "After FFT" << endl;
   print(data, "B");
   //cout << "Done" << endl;
}
