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
#include <bitset>
#include <string>

using namespace std;
const double pi      = 3.14159265358979323846264;

void printAll(vector<complex<double> > data,string label);

// Une etape de la FFT sequentielle
void stepSeq(vector<complex<double> >& data,
             complex<double>& w,int d) {
	int nbPE,myPE;
    MPI_Comm_size(MPI_COMM_WORLD,&nbPE);
    MPI_Comm_rank(MPI_COMM_WORLD,&myPE);

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

void DecToBinReverse(int id, char* id_bin){
	int i = 0;

    do
    {
        if ( (id & 1) == 0 )
            id_bin[i] = '0';
        else
            id_bin[i] = '1';

        id >>= 1;
        i++;
    } while ( id );
}

void bitReversedPar(vector<complex<double> >& data) {
	int nbPE,myPE;
    MPI_Comm_size(MPI_COMM_WORLD,&nbPE);
    MPI_Comm_rank(MPI_COMM_WORLD,&myPE);

    for (int j = 0;j < data.size();j++){

		//Calcul de l'id réel (en global sur tous les processeurs)
		int idreeldata = j+data.size()*myPE;

		//Calcul sur combien de bits coder les addresses
		int nbbits = log2(nbPE*data.size());	//fois data.size car plusieur data par proc

		//alloue pour l id en binaire
		char* id_bin = (char *)calloc((nbbits+1),sizeof(char));

		//initialiste l id
		for(int i = 0; i < nbbits; i++)
			id_bin[i]='0';

		//convertie l id en binaire et l inverse
		DecToBinReverse(idreeldata,id_bin);

		int idbinome = 0;

		//converti l id binaire inverse en decimal
		int x = nbbits-1;

		for(int i = 0; i < nbbits; i++){
            char tmp[2];
            tmp[0] = id_bin[i];
            tmp[1] = 0;

            int value = atoi(tmp);
			idbinome += pow(2,x) * value;
			x--;
		}

		//Calcul du processeur contenant le binome
		int idprocbinome = idbinome/data.size();

		//stock nombre complex dans un buffer pour l'envoie
		int buf_size = 2;
		double* send_buf = new double[buf_size];
		double* recv_buf = new double[buf_size];

  		send_buf[0] = data[j].real();
  		send_buf[1] = data[j].imag();

		//si plus petit que son binome attend et ensuite envoie (pour eviter collision)
		if(myPE > idprocbinome){
			//recois de son binome
			MPI_Recv(recv_buf, buf_size, MPI_DOUBLE, idprocbinome, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			//envoie a son binome
			MPI_Send(send_buf, buf_size, MPI_DOUBLE, idprocbinome , 0, MPI_COMM_WORLD);
		}
        else{
			//envoie a son binome
			MPI_Send(send_buf, buf_size, MPI_DOUBLE, idprocbinome , 0, MPI_COMM_WORLD);

			//recois de son binome
			MPI_Recv(recv_buf, buf_size, MPI_DOUBLE, idprocbinome, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}

		//met à jour data local avec celle de son binome
		complex<double> tmpLoc(recv_buf[0], recv_buf[1]);
		data[j] = tmpLoc;

        free(id_bin);
	}
}

void swapPar(complex<double>& loc, int proc) {
    // Get process id
    int myPE;
    MPI_Comm_rank(MPI_COMM_WORLD,&myPE);

    // Init of buffer to send and receive the complex number to send
    double* send_buf = new double[2];
    double* recv_buf = new double[2];

    // Put in send buffer the complex number
    send_buf[0] = loc.real();
    send_buf[1] = loc.imag();

    // Determine an order for send and receive
    if (proc > myPE) {
        MPI_Send(send_buf, 2, MPI_DOUBLE, proc, 0, MPI_COMM_WORLD);
        MPI_Recv(recv_buf, 2, MPI_DOUBLE, proc, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
        MPI_Recv(recv_buf, 2, MPI_DOUBLE, proc, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(send_buf, 2, MPI_DOUBLE, proc, 0, MPI_COMM_WORLD);
    }

    // Update the loc data
    complex<double> tmpLoc(recv_buf[0], recv_buf[1]);
    loc = tmpLoc;
}

// Une etape de la FFT parallele
void stepPar(vector<complex<double> >& data, complex<double>& w, int d) {
    int nbPE,myPE;
    MPI_Comm_size(MPI_COMM_WORLD,&nbPE);
    MPI_Comm_rank(MPI_COMM_WORLD,&myPE);

    int nloc = data.size();
    complex<double> wk = 1;

    for (int i = 1; i <= ((myPE * nloc) % d); i++) {
        wk = w*wk;
    }

    int indGlobal;
    complex<double> oldElement;

    for (int k = 0; k < nloc; k++) {
        indGlobal = k + myPE * nloc;    // Indice global dans la partie paire

        if ((indGlobal & (0x1 << (int)log2(d))) == 0) {
            oldElement = data[k];
            swapPar(data[k], myPE ^ (d/nloc));
            data[k] = data[k] + oldElement;
        }
        else {
            data[k] = wk * data[k];
            oldElement = data[k];
            swapPar(data[k], myPE ^ (d/nloc));
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

    complex<double> w = polar(1.0, -pi);

    bitReversedPar(data);

    int size = 0;
    for (int etape = 0; etape < log2(nbPE*data.size()); etape++) {
        size = pow(2, etape);

        // Partie sans communication
        if (size < data.size()) {
            stepSeq(data, w, size);
        }

        // Partie avec communication
        else {
            stepPar(data, w, size);
        }
    }
}

// iFFT parallèle 1D
void ifftPar(vector<complex<double> >& data) {
    int nbPE,myPE;
    MPI_Comm_size(MPI_COMM_WORLD,&nbPE);
    MPI_Comm_rank(MPI_COMM_WORLD,&myPE);

    complex<double> w = polar(1.0, pi);

    bitReversedPar(data);

    int size = 0;
    for (int etape = 0; etape < log2(nbPE*data.size()); etape++) {
        size = pow(2, etape);

        // Partie sans communication
        if (size < data.size()) {
            stepSeq(data, w, size);
        }

        // Partie avec communication
        else {
            stepPar(data, w, size);
        }
    }

    // Division by N (total size of data)
    for (int i = 0; i < data.size(); i++) {
        data[i] *= (double)1 / (double)(data.size() * nbPE);
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
      //cout << label << " = [";
      //for (int k=0;k<nbPE*buf_size;k++) {
      //   if ((k+1)%2 == 0) cout << "+(" << recv_buf[k] << "i);";
      //   else cout << recv_buf[k] << "";
      //}
      //cout << "];" << endl;
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
   MPI_Barrier(MPI_COMM_WORLD);
   printAll(data,"A");
   fftPar(data);
   MPI_Barrier(MPI_COMM_WORLD);
   printAll(data,"B");
   MPI_Barrier(MPI_COMM_WORLD);
   ifftPar(data);
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();
   return 0;
}
