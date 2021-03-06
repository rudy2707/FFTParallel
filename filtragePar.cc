/*
 * =====================================================================================
 *
 *       Filename:  filtragePar.cc
 *
 *    Description:  Read a PGM file and apply the FFT on it in order
 *                  to filter it.
 *
 *        Version:  1.0
 *        Created:  01/07/2016 03:30:40 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Axel Fahy (),
 *   Organization:  HES-SO hepia section ITI
 *
 * =====================================================================================
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <string>
#include <stdlib.h>
#include <mpi.h>

using namespace std;
const double pi      = 3.14159265358979323846264;

#define NAME_FILTERED "_filtre"

vector<complex<double> > bufferToComplex(double* buffer, int size);
double* complexToBuffer(vector<complex<double> > vec);
void processFFT(vector<complex<double> > &data, int imgSize, bool isRow, bool isFFT);

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

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    int nbPE,myPE;
    MPI_Comm_size(MPI_COMM_WORLD, &nbPE);
    MPI_Comm_rank(MPI_COMM_WORLD, &myPE);

    // Check the number of parameters
    if (argc < 2) {
        // Help message
        cerr << "Usage: " << argv[0] << " file.pgm filter_percent" << endl;
        return 1;
    }
    // Get the filename
    string filename = argv[1];

    // Get the filter percent
    double percent = atof(argv[2]);

    string header[4];

    int imgSize = 0;

    // Store the image into a vector of complex<double>
    vector<complex<double> > data;

    // Process 0 read the file
    if (myPE == 0) {
        // Create an input file stream
        ifstream in(filename.c_str(), ios::in);

        cout << "Reading file : " << filename << endl;

        // Get the header (4 lines)
        for (int i = 0; i < 4; i++) {
            in >> header[i];
        }

        // Get the size of the file
        imgSize = atoi(header[1].c_str());

        // Variable to store the current value when parsing the file
        int number;

        while (in >> number) {
            // Add the number into the array
            // with the complex part set to 0
            data.push_back(complex<double>(number, 0.0));
        }

        // Close the filestream
        in.close();
    }

    // Broadcast the size of the image to every process
    MPI_Bcast((void*)&imgSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    int sizeLocal = imgSize / nbPE;

    // Apply FFT on image's rows
    processFFT(data, imgSize, true, true);

    MPI_Barrier(MPI_COMM_WORLD);

    // Apply FFT on image's columns
    processFFT(data, imgSize, false, true);

    MPI_Barrier(MPI_COMM_WORLD);

    // Filter frequencies by setting corner to 0
    int nbFilter = imgSize * percent / 100;
    if (myPE == 0) {
        for (int i = 0; i < imgSize; i++) {
            for (int j = 0; j < imgSize; j++) {
                data[i + j * imgSize] = 0;
                if (j == nbFilter) {
                    j += imgSize - nbFilter * 2;
                }
                if (i == nbFilter) {
                    i += imgSize - nbFilter * 2;
                }
            }
        }
    }

    // Apply iFFT on image's rows
    processFFT(data, imgSize, true, false);

    MPI_Barrier(MPI_COMM_WORLD);

    // Apply iFFT on image's columns
    processFFT(data, imgSize, false, false);

    MPI_Barrier(MPI_COMM_WORLD);

    // Only root process writes the filtered image
    if (myPE == 0) {
        ofstream out;
        string ext = filename.substr(filename.find("."));
        string base = filename.substr(0, filename.find("."));
        string fout = base + NAME_FILTERED + ext;
        cout << "Filename out : " << fout << endl;

        // Write the filter file
        out.open(fout.c_str());
        for (int i = 0; i < 4; i++) {
            out << header[i];
            out << "\n";
        }
        for (int i = 0; i < data.size(); i++) {
            out << (int)(abs(real(data[i])));
            //out << real(data[i]);
            out << "\n";
        }
        out.close();
    }

    MPI_Finalize();
    return 0;
}

/**
 * @brief Process the FFT, split the row / column and send them to all the processes.
 *
 * This is the process number 0 which dispatch data to other processes.
 *
 * If isRow is true, apply on rows, else on column
 * If isFFT is true, apply FFT, else apply iFFT
 *
 * @param data      All data
 * @param imgSize   Size of image (one dimension)
 * @param isRow     If we are processing the rows
 * @param isFFT     If we are doing the fft of the ifft
 */
void processFFT(vector<complex<double> > &data, int imgSize, bool isRow, bool isFFT) {
    int nbPE,myPE;
    MPI_Comm_size(MPI_COMM_WORLD, &nbPE);
    MPI_Comm_rank(MPI_COMM_WORLD, &myPE);

    int sizeLocal = imgSize / nbPE;

    double* recv_buf;
    double* send_buf;

    if (myPE == 0) {
        // Send parts of row to other process and then process its own part
        // Apply FFT on each row and then on each column
        // Run over each row and apply the FFT
        for (int i = 0; i < imgSize; i++) {
            vector<complex<double> > rowOrCol;  // Curent row or column depending on the operation we are running
            for (int j = 0; j < imgSize; j++) {
                if (isRow)
                    rowOrCol.push_back(data[i * imgSize + j]);
                else
                    rowOrCol.push_back(data[j * imgSize + i]);
            }
            // Split the row by the number of process and send the data to them
            // Send all the parts and then process its own
            for (int k = 1; k < nbPE; k++) {
                vector<complex<double> > dataLocal(rowOrCol.begin() + (k * sizeLocal), rowOrCol.begin() + (k * sizeLocal) + sizeLocal);
                double* send_buf = complexToBuffer(dataLocal);
                MPI_Send(send_buf, sizeLocal * 2, MPI_DOUBLE, k, 0, MPI_COMM_WORLD);
                delete send_buf;
            }

            // Process its own part
            vector<complex<double> > dataLocal(rowOrCol.begin(), rowOrCol.begin() + sizeLocal);

            if (isFFT) {
                fftPar(dataLocal);
            }
            else {
                ifftPar(dataLocal);
            }

            for (int j = 0; j < sizeLocal; j++) {
                if (isRow) {
                    data[i * imgSize + j] = dataLocal[j];
                }
                else {
                    data[j * imgSize + i] = dataLocal[j];
                }
            }

		    recv_buf = new double[sizeLocal * 2];
            // Receive data from other process
            for (int k = 1; k < nbPE; k++) {
                MPI_Recv(recv_buf, sizeLocal * 2, MPI_DOUBLE, k, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                // Modification of computed value in data vector
                vector<complex<double> > dataPartFiltered = bufferToComplex(recv_buf, sizeLocal * 2);
                for (int j = 0; j < sizeLocal; j++) {
                    if (isRow)
                        data[i * imgSize + (k * sizeLocal) + j] = dataPartFiltered[j];
                    else
                        data[imgSize * k * sizeLocal + i + j * imgSize] = dataPartFiltered[j];
                }
            }
            delete recv_buf;
        }
    }
    else {
        // Reception of data and FFT process
		double* recv_buf = new double[sizeLocal * 2];
        // Process each rows
        for (int i = 0; i < imgSize; i++) {
            // Recreate the vector of complex
            MPI_Recv(recv_buf, sizeLocal * 2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            vector<complex<double> > dataLocal = bufferToComplex(recv_buf, sizeLocal * 2);

            if (isFFT)
                fftPar(dataLocal);
            else
                ifftPar(dataLocal);

            double* send_buf = complexToBuffer(dataLocal);
            MPI_Send(send_buf, sizeLocal * 2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            delete send_buf;
        }
        delete recv_buf;
    }
}

/**
 * @brief Transformation of vector of complex into an array of double
 *
 * The array of double as twice as much value as the vector because we need to
 * separate the imag and real part.
 *
 * The memory is allocated in the function and needs to be free outside.
 *
 * @param vec Vector to transform
 *
 * @return An array of double with values of vector
 */
double* complexToBuffer(vector<complex<double> > vec) {
    // Buffer with the double of the size for the real and the imaginary part
    double* buffer = new double[vec.size() * 2];
    for (int i = 0; i < vec.size(); i++) {
        buffer[2 * i] = vec[i].real();
        buffer[2 * i + 1] = vec[i].imag();
    }
    return buffer;
}

/**
 * @brief Transformation of an array of double into a vector of complex number
 *
 * @param buffer Array of double to transform
 * @param size   Size of the array (double of the vector size)
 *
 * @return
 */
vector<complex<double> > bufferToComplex(double* buffer, int size) {
    vector<complex<double> > vec;
    for (int i = 0; i < size / 2; i++) {
        complex<double> tmpComplex(buffer[2 * i], buffer[2 * i + 1]);
        vec.push_back(tmpComplex);
    }
    return vec;
}





