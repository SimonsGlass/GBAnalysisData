#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
using namespace std;




// Interlace two files: add softness field to LAMMPS dump files

int main(int argc, char* argv[]) {

        // Verify number of command line arguments correct
        if ( argc != 3 ) { // argc should be 2 for correct execution
                cout<<"Usage: "<< argv[0] <<" <filenameDump> <filenameGBIDField> \n";
                return 0;
        }

        // Set up constants
        cout << "Please check nndist, nheader, initial file length guess." << endl;
        int nheader = 9; // header length of "small" files
        long idmaxallowed = 4000000; // max length of arrays required
        int* softnessarr = new int[idmaxallowed];
        int cluster, softness;
        long particleID;
        char newsuffix[20]; sprintf(newsuffix, ".GBID"); // output files will have this suffix
        //string oldfields("ITEM: ATOMS id x y z phop soft Cluster StructureType");
        string oldfields("ITEM: ATOMS id x y z phop soft");
        //string newfields("ITEM: ATOMS id x y z phop soft Cluster StructureType GBID");
        string newfields("ITEM: ATOMS id x y z phop soft GBID");

        string readinline;
        double posx, posy, posz = -1.0, phop, soft;
        int StructureType;

        // Open the file with the softness array and load it into an array, softnessarr
        ifstream soft_file; std::string sfilename(argv[2]); soft_file.open( sfilename.c_str() );
        if ( !soft_file.is_open() ) {
                cout<<"Could not open softness file\n";
                return 0;
        }
        particleID =1;
        while (!soft_file.eof()) {
                soft_file >> softness;
                // The array is indexed by actual particle ID (starting from 1)
                softnessarr[particleID++] = softness;
                if (particleID < 5) printf("softness value read is %d \n",softness);
        }
        soft_file.close();


        // Open the dump and output files
        ifstream dump_file; std::string dfilename(argv[1]); dump_file.open( dfilename.c_str() );
        if ( !dump_file.is_open() ) {
                cout<<"Could not open LAMMPS file\n";
                return 0;
        }
        string outfilename = dfilename; outfilename.append(newsuffix);
        ofstream outfile; outfile.open(outfilename.c_str()); outfile.precision(9);
        if ( !outfile.is_open() ) {
                cout<<"Could not open new dump file " << outfile << endl;
                return 0;
        }


        // Copy the LAMMPS header over
        for (int jline=0; jline<nheader; jline++) {
                getline(dump_file,readinline);
                if (jline == 0) cout << "First header line: " << readinline << endl;
                if (0 == readinline.compare(oldfields)) readinline = newfields;
                outfile << readinline << endl;
        }
        cout << "Last line of header copied into outfile:" << endl;
        cout << readinline << endl;


        // Copy the data fields over but adding the softness of the given particleID
        while (!dump_file.eof()) {
                posz=-2e8;
                dump_file >> particleID >> posx >> posy >> posz >> phop >> soft;// >> cluster >> StructureType;
                if (posz > -1e8)
                    outfile << particleID << " " << posx << " " << posy << " " << posz << " " 
                            << phop << " " << soft << " " // << cluster << " " << StructureType << " "
                            << softnessarr[particleID] << endl;
        }

        dump_file.close(); outfile.close();

        return 0;
}
