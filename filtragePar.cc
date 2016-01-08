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

#include "fftPar1D.h"

using namespace std;

#define NAME_FILTERED "_filtre"

vector<complex<double> > bufferToComplex(double* buffer, int size);
double* complexToBuffer(vector<complex<double> > vec);

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
        //imgSize = header[1]; // WTF ??? TODO
        imgSize = 256;

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
    cout << "[Process " << myPE << "] Img size : " << imgSize << endl;
    int sizeLocal = imgSize / nbPE;
/*
    if (myPE == 0) {
        // Send parts of row to other process and then process its own part
        // Apply FFT on each row and then on each column
        // Run over each row and apply the FFT
        for (int i = 0; i < imgSize; i++) {
            vector<complex<double> > row;
            for (int j = 0; j < imgSize; j++) {
                row.push_back(data[i * imgSize + j]);
            }
            // Split the row by the number of process and send the data to them
            // Send all the parts and then process its own
            for (int k = 1; k < nbPE; k++) {
                vector<complex<double> > rowLocal(row.begin() + (k * sizeLocal), row.begin() + (k * sizeLocal) + sizeLocal );
                double* send_buf = complexToBuffer(rowLocal);
                MPI_Send(send_buf, sizeLocal * 2, MPI_DOUBLE, k, 0, MPI_COMM_WORLD);
                free(send_buf);
            }

            // Process its own part
            vector<complex<double> > rowLocal(row.begin(), row.begin() + sizeLocal);
            fftPar(rowLocal);
            // Gather data from each process and replace data with value computed in row
            MPI_Barrier(MPI_COMM_WORLD);
		    double* recv_buf = new double[imgSize * 2];
            MPI_Gather(send_buf, buf_size, MPI_DOUBLE, recv_buf, buf_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            for (int j = 0; j < imgSize; j++) {
                complex<double> tmpComplex(recv_buf[2 * i], recv_buf[2 * i + 1]);
                rowLocal[i] = tmpComplex;
            }
            for (int j = 0; j < imgSize; j++) {
                data[i * imgSize + j] = row[j];
            }
        }

    }
    else {
        // Reception of data and FFT process
		double* recv_buf = new double[imgSize * 2];
        // Process each rows
        for (int i = 0; i < imgSize; i++) {
            // Recreate the vector of complex
            MPI_Recv(recv_buf, sizeLocal * 2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            vector<complex<double> > rowLocal = bufferToComplex(recv_buf, imgSize * 2);
            fftPar(rowLocal);
            double* send_buf = complexToBuffer(rowLocal);
            MPI_Send(send_buf, sizeLocal * 2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            free(send_buf);
        }
        free(recv_buf);
    }

    // Run over each column and apply the FFT
    for (int j = 0; j < imgSize; j++) {
        vector<complex<double> > column;
        for (int i = 0; i < imgSize; i++) {
            column.push_back(data[i * imgSize + j]);
        }
        //fftPar(column);
        // Replace data with value computed in column
        for (int i = 0; i < imgSize; i++) {
            data[i * imgSize + j] = column[j];
        }
    }
*/
    // prendre la valeur absolue des reel pour reobtenir l'image (% 255 ???)
    // Only root process wrties the filtered image
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
            out << real(data[i]);
            out << "\n";
        }
        out.close();
    }

    MPI_Finalize();
    return 0;
}

double* complexToBuffer(vector<complex<double> > vec) {
    // Buffer with the double of the size for the real and the imaginary part
    double* buffer = new double[vec.size() * 2];
    for (int i = 0; i < vec.size(); i++) {
        buffer[2 * i] = vec[i].real();
        buffer[2 * i + 1] = vec[i].imag();
    }
    return buffer;
}

vector<complex<double> > bufferToComplex(double* buffer, int size) {
    vector<complex<double> > vec;
    for (int i = 0; i < size / 2; i++) {
        complex<double> tmpComplex(buffer[2 * i], buffer[2 * i + 1]);
        vec[i] = tmpComplex;
    }
    return vec;
}





