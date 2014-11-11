
#include <cstdlib>
#include <iostream>
#include "Trajectory.h"
#include "Utils.h"

int main(int argc, char *argv[]) {

    const double start = 0.0;
    double *g;
    double r2;
    int frame;
    int i;
    int ig;
    int j;
    int nFrames;
    int nGrp1;
    int nGrp2;
    int nOW;
    int total;
    matrix box;
    ofstream oFS;
    rvec atomi;
    rvec atomj;
    rvec dx;

    if (!argc == 8) {
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

    cout << "Trajectory file:      " << xtcfile << endl;
    cout << "Index file:           " << ndxfile << endl;
    cout << "Output file:          " << outfile << endl;
    cout << "Group 1:              " << grp1 << endl;
    cout << "Group 2:              " << grp2 << endl;
    cout << "Exclusion distance:   " << rexcl << endl;
    cout << "Bin width:            " << binwidth << endl;
    cout << "Location of last bin: " << end << endl << endl;

    const int nBins = (end-start)/binwidth + 1;
    g = new double[nBins];

    Trajectory traj(xtcfile, ndxfile);

    nFrames = traj.GetNFrames();
    nGrp1 = traj.GetNAtoms(grp1);
    nGrp2 = traj.GetNAtoms(grp2);
    total = 0;

    for (frame = 0; frame < nFrames; frame++) {

        traj.GetBox(frame,box);

        for (i = 0; i < nGrp1; i++) {

            traj.GetXYZ(frame,grp1,i,atomi);

            for (j = 0; j < nGrp2; j++) {

                traj.GetXYZ(frame,grp1,j,atomj);
                dx[X] = atomi[X] - atomj[X];
                dx[Y] = atomi[Y] - atomj[Y];
                dx[Z] = atomi[Z] - atomj[Z];
                pbc(dx,box);
                r2 = dot(dx,dx);
                if (r2 > rexcl2) {
                    ig = sqrt(r2/binwidth);
                    g[ig] += 1.0;
                    total++;
                }

            }

        }

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
    for (i = 0; i < nBins; i++) {
        oFS << ((double) i) * binwidth + start << "  " << g[i]/total << endl;
    }
    oFS.close();

    return 0;
}
