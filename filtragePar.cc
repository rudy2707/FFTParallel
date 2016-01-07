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

#include "fftPar1D.h"

using namespace std;

#define NAME_FILTERED "_filtre"

int main(int argc, char* argv[])
{
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

    // Store the image into a vector of complex<double>
    vector<complex<double> > data;

    // Create an input file stream
	ifstream in(filename.c_str(), ios::in);

    // Get the header (4 lines)
    string header[4];

    for (int i = 0; i < 4; i++) {
        in >> header[i];
    }

    // Variable to store the current value when parsing the file
    int number;

    while (in >> number) {
        // Add the number into the array
        // with the complex part set to 0
        data.push_back(complex<double>(number, 0.0));
    }

    // Close the filestream
    in.close();

    // Apply FFT on each row and then on each column
    // Run over each row and apply the FFT
    for (int i = 0; i < 256; i++) {
        vector<complex<double> > row;
        for (int j = 0; j < 256; j++) {
            row.push_back(data[i * 256 + j]);
        }
        //fftPar(row);
        // Replace data with value computed in row
        for (int j = 0; j < 256; j++) {
            data[i * 256 + j] = row[j];
        }
    }

    // Run over each column and apply the FFT
    for (int j = 0; j < 256; j++) {
        vector<complex<double> > column;
        for (int i = 0; i < 256; i++) {
            column.push_back(data[i * 256 + j]);
        }
        //fftPar(column);
        // Replace data with value computed in column
        for (int i = 0; i < 256; i++) {
            data[i * 256 + j] = column[j];
        }
    }

    // prendre la valeur absolue des reel pour reobtenir l'image (% 255 ???)
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

    return 0;
}





