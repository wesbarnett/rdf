
#include <cstdlib>
#include <iostream>
#include <string>
#include <cstring>
#include "Trajectory.h"
#include "Utils.h"

int main(int argc, char *argv[]) {

    const double start = 0.0;
    double r2;
    double binvol;
    double boxvol;
    double r;
    int frame;
    int i;
    int ig;
    int j;
    int nOW;
    matrix box;
    ofstream oFS;
    rvec atomi;
    rvec atomj;
    rvec dx;

    if (argc != 9) {
        cout << "Usage: " << endl;
        cout << "  " << argv[0] << " xtcfile ndxfile outfile group1 group2 excldist binwidth endofrdf" << endl;
        return -1;
    }

    const string xtcfile = argv[1];
    const string ndxfile = argv[2];
    const string outfile = argv[3];
    const string grp1 = argv[4];
    const string grp2 = argv[5];
    const double rexcl = atof(argv[6]);
    const double rexcl2 = rexcl * rexcl;
    const double binwidth = atof(argv[7]);
    const double end = atof(argv[8]);
    const double end2 = end * end;

    cout << "Trajectory file:      " << xtcfile << endl;
    cout << "Index file:           " << ndxfile << endl;
    cout << "Output file:          " << outfile << endl;
    cout << "Group 1:              " << grp1 << endl;
    cout << "Group 2:              " << grp2 << endl;
    cout << "Exclusion distance:   " << rexcl << endl;
    cout << "Bin width:            " << binwidth << endl;
    cout << "Location of last bin: " << end << endl;

    const int nBins = (end-start)/binwidth + 1;
    double g[nBins];

    Trajectory traj(xtcfile, ndxfile);

    const int nFrames = traj.GetNFrames();
    const int nGrp1 = traj.GetNAtoms(grp1);
    const int nGrp2 = traj.GetNAtoms(grp2);

    // Constant volume
    traj.GetBox(0,box);

    #pragma omp parallel for private(frame,i,j,atomi,atomj,dx,r2,ig)
    for (frame = 0; frame < nFrames; frame++) {

        if (frame % 100 == 0 ) cout << "Calculating frame: " << frame << endl;;

        for (i = 0; i < nGrp1; i++) {

            traj.GetXYZ(frame,grp1,i,atomi);

            for (j = 0; j < nGrp2; j++) {

                traj.GetXYZ(frame,grp1,j,atomj);
                dx[X] = atomi[X] - atomj[X];
                dx[Y] = atomi[Y] - atomj[Y];
                dx[Z] = atomi[Z] - atomj[Z];
                pbc(dx,box);
                r2 = dot(dx,dx);
                if (r2 > rexcl2 && r2 < end2) {
                    ig = ceil(sqrt(r2)/binwidth);
                    g[ig] += 1.0;
                }

            }

        }

    }

    boxvol = volume(box);

    for (i = 0; i < nBins; i++) {
        r = (double) i + 0.5;
        binvol = pow(r,3) - pow((r-1.0),3);
        binvol *= 4.0/3.0 * M_PI * pow(binwidth,3);
        g[i] *= boxvol / (nGrp1 * nGrp2 * binvol * nFrames);
    }

    oFS.open(outfile.c_str());
    oFS << "# Trajectory file:      " << xtcfile << endl;
    oFS << "# Index file:           " << ndxfile << endl;
    oFS << "# Output file:          " << outfile << endl;
    oFS << "# Group 1:              " << grp1 << endl;
    oFS << "# Group 2:              " << grp2 << endl;
    oFS << "# Exclusion distance:   " << rexcl << endl;
    oFS << "# Bin width:            " << binwidth << endl;
    oFS << "# Location of last bin: " << end << endl << endl;
    oFS << fixed << setprecision(6);
    for (i = 0; i < nBins; i++) {
        oFS << ((double) i) * binwidth + start << "  " << g[i] << endl;
    }
    oFS.close();

    return 0;
}
