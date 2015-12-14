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
#include <mpi.h>

using namespace std;

// FFT parallèle 1D
void fftPar(vector<complex<double> >& data) { /*à compléter*/}

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
   fftPar(data);
   printAll(data,"%B\n");
   MPI_Finalize();
   return 0;
}
