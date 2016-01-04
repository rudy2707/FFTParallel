// fft parallèle 1D
// Implementation d'une transformee de Fourier rapide
#include <iostream>
#include <complex>
#include <algorithm>
#include <vector>
#include <bitset>
#include <cstdlib>
#include <ctime>
#include <cstring>
#include <cmath>
#include <mpi.h>

using namespace std;

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

void bitReversedPar(vector<complex<double>>& data) {

	int nbPE,myPE;
    MPI_Comm_size(MPI_COMM_WORLD,&nbPE);
    MPI_Comm_rank(MPI_COMM_WORLD,&myPE);
    
    cout << ("id en des : ") << myPE;

    // Converti son id en binaire
    string MyPEbit = DestoBin(myPE);
    
    cout << ("id en bit : ") << MyPEbit;
    
    //inverse les bits
    string MyPEbitinverse[MyPEbit.length()];
    
    int j = 0;
    for(i = MyPEbit.length; i >= MyPEbit.length; i---){
		    MyPEbitinverse[j]=MyPEbit[i];
		    j++;		    
    }
    
    cout << ("id inverse en bit : ") << MyPEbitinverse;
    
    // converti les bits en decimal
    int Id_Proc_binome;
    for(i = 0; i < MyPEbitinverse.length; i++){
		    Id_Proc_binome += (2^i)*MyPEbitinverse[i];		    
    } 
    
    cout << ("id inverse en dec : ") << Id_Proc_binome;
}

string DecToBin(int number){
    string result = "";

    do
    {
        if ( (number & 1) == 0 )
            result += "0";
        else
            result += "1";

        number >>= 1;
    } while ( number );

    reverse(result.begin(), result.end());
    return result;
}

void swapPar(complex<double> loc, int proc) {
    // TODO 
}

// Une etape de la FFT parallele
void stepPar(vector<complex<double> >& data, complex<double>& w, int d) {
    int nbPE,myPE;
    MPI_Comm_size(MPI_COMM_WORLD,&nbPE);
    MPI_Comm_rank(MPI_COMM_WORLD,&myPE);

    int nloc = nbPE / data.size();
    int wk = 1;
    
    for (int i = 1; i < ((myPE * nloc) % size); i++) {
        wk = w*wk;
    }

    int indGlobal;
    complex<double> oldElement;
    for (int k = 0; k < nloc-1; k++) {
        indGlobal = k + myPE * nloc;    // Indice global dans la partie paire

        if ((indGlobal & (0x1 << (int)log2(size))) == 0) {
            oldElement = data[k];
            swapPar(data[k], myPE ^ (size/nloc));
            data[k] = data[k] + oldElement;
        }
        else {
            data[k] = wk * data[k];
            oldElement = data[k];
            swapPar(data[k], myPE ^ (size/nloc));
            data[k] = data[k] - oldElement;
        }

        wk *= w;
    }

    w = sqrt(w);
}

// FFT parallèle 1D
void fftPar(vector<complex<double> >& data) { 
    int nbPE,myPE;
    MPI_Comm_size(MPI_COMM_WORLD,&nbPE);
    MPI_Comm_rank(MPI_COMM_WORLD,&myPE);

    complex<double> w = -1;
    data = bitReversedPar(data);

    int size = 0;
    for (int etape = 0; etape < log2(nbPE); etape++) {
        size = pow(2, etape);
        
        // Partie sans communication
        if (size < (nbPE/data.size()) {
            stepSeq(data, w, size);
        }

        // Partie avec communication
        else {
            stepPar(data, w, size);
        }
    }
}

// Initialisation aléatoire d'un vecteur de nombres complexes
void randInit(vector<complex<double> >& data,double min, double max) {
   for (int k=0;k<data.size();k++)
      data[k] = complex<double>(min+(max-min)*drand48(),min+(max-min)*drand48());
}

// Affichage d'un vecteur de nombres complexes
void printAll(vector<complex<double> > data,string label) {
   int nbPE,myPE;
   MPI_Comm_size(MPI_COMM_WORLD,&nbPE);
   MPI_Comm_rank(MPI_COMM_WORLD,&myPE);
   int buf_size = data.size()*2;
   double* send_buf = new double[buf_size];
   for (int k=0;k<buf_size/2;k++) {
      send_buf[2*k] = data[k].real();
      send_buf[2*k+1] = data[k].imag();
   }
   double* recv_buf = new double[nbPE*buf_size];
   MPI_Gather(send_buf,buf_size,MPI_DOUBLE,recv_buf,buf_size,MPI_DOUBLE,0,MPI_COMM_WORLD);
   if (myPE == 0) {
      cout << label;
      for (int k=0;k<nbPE*buf_size;k++) {
         cout << recv_buf[k] << " ";
         if ((k+1)%2 == 0) cout << endl;
      }
   }
   delete recv_buf; 
}

int main(int argc,char ** argv) {
   MPI_Init(&argc,&argv);
   int myPE;
   MPI_Comm_rank(MPI_COMM_WORLD,&myPE);
   int seed = atoi(argv[1]);
   srand48(seed+myPE);
   int nloc = atoi(argv[2]);
   vector<complex<double> > data(nloc);
   randInit(data,0.0,100.0);
   printAll(data,"%A\n");
   //fftPar(data);
   printAll(data,"%B\n");
   MPI_Finalize();
   return 0;
}
