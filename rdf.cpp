
#include <cstdlib>
#include <iostream>
#include <string>
#include <cstring>
#include "Trajectory.h"
#include "Utils.h"

void doRdf(Trajectory &traj, string grp, double rexcl2, double end2, double binwidth, vector <double> &g) {

    const int nFrames = traj.GetNFrames();
    const int nGrp = traj.GetNAtoms(grp);
    double boxvol;
    double r2;
    int frame;
    int i;
    int ig;
    int j;
    matrix box;
    rvec atomi;
    rvec atomj;
    rvec dx;

    #pragma omp parallel for schedule(guided) private(frame,i,j,atomi,atomj,dx,r2,ig)
    for (frame = 0; frame < nFrames; frame++) {

        traj.GetBox(frame,box);
        boxvol = volume(box);

        for (i = 0; i < nGrp-1; i++) {

            traj.GetXYZ(frame,grp,i,atomi);

            for (j = i+1; j < nGrp; j++) {

                traj.GetXYZ(frame,grp,j,atomj);
                dx[X] = atomi[X] - atomj[X];
                dx[Y] = atomi[Y] - atomj[Y];
                dx[Z] = atomi[Z] - atomj[Z];
                pbc(dx,box);
                r2 = dot(dx,dx);
                if (r2 > rexcl2 && r2 < end2) {
                    ig = ceil(sqrt(r2)/binwidth);
                    g.at(ig) += boxvol;
                }

            }

        }

    }

    return;

}

void doRdf(Trajectory &traj, string grp1, string grp2, double rexcl2, double end2, double binwidth, vector <double> &g) {

    const int nFrames = traj.GetNFrames();
    const int nGrp1 = traj.GetNAtoms(grp1);
    const int nGrp2 = traj.GetNAtoms(grp2);
    double boxvol;
    double r2;
    int frame;
    int i;
    int ig;
    int j;
    matrix box;
    rvec atomi;
    rvec atomj;
    rvec dx;

    #pragma omp parallel for schedule(guided) private(frame,i,j,atomi,atomj,dx,r2,ig)
    for (frame = 0; frame < nFrames; frame++) {

        traj.GetBox(frame,box);
        boxvol = volume(box);

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
                    g.at(ig) += boxvol;
                }

            }

        }

    }

    return;
}

void normalize(Trajectory &traj, string grp1, string grp2, int nBins, double binwidth, vector <double> &g) {

    const double f = 4.0/3.0 * M_PI;
    const int nFrames = traj.GetNFrames();
    const int nGrp1 = traj.GetNAtoms(grp1);
    const int nGrp2 = traj.GetNAtoms(grp2);
    double binvol;
    double r;
    int i;

    for (i = 0; i < nBins; i++) {
        r = (double) i;
        binvol = pow(r,3) - pow((r-1.0),3);
        binvol *= f * pow(binwidth,3);
        g.at(i) /= ((nGrp1-1) * nGrp2 * binvol * nFrames);
    }

}

int main(int argc, char *argv[]) {

    const double start = 0.0;
    double binwidth;
    double end;
    double end2;
    double rexcl;
    double rexcl2;
    ifstream iFS;
    ofstream oFS;
    string configfile;
    string grp1;
    string grp2;
    string ndxfile;
    string outfile;
    string xtcfile;

    if (argc != 2) {
        cout << "Usage: " << endl;
        cout << "  " << argv[0] << " configfile" << endl;
        return -1;
    }

    configfile = argv[1];

    cout << "Reading " << configfile << "." << endl;
    iFS.open(configfile.c_str());
    iFS >> xtcfile;
    iFS >> ndxfile;
    iFS >> outfile;
    iFS >> grp1;
    iFS >> grp2;
    iFS >> rexcl;
    iFS >> binwidth;
    iFS >> end;
    iFS.close();

    cout << "Trajectory file:      " << xtcfile << endl;
    cout << "Index file:           " << ndxfile << endl;
    cout << "Output file:          " << outfile << endl;
    cout << "Group 1:              " << grp1 << endl;
    cout << "Group 2:              " << grp2 << endl;
    cout << "Exclusion distance:   " << rexcl << endl;
    cout << "Bin width:            " << binwidth << endl;
    cout << "Location of last bin: " << end << endl;
    rexcl2 = rexcl*rexcl;
    end2 = end * end;

    const int nBins = (end-start)/binwidth + 1;
    vector <double> g(nBins,0.0);

    Trajectory traj(xtcfile, ndxfile);

    const int nFrames = traj.GetNFrames();
    const int nGrp1 = traj.GetNAtoms(grp1);
    const int nGrp2 = traj.GetNAtoms(grp2);

    if (grp1 == grp2) {
        doRdf(traj,grp1,rexcl2,end2,binwidth,g);
    } else {
        doRdf(traj,grp1,grp2,rexcl2,end2,binwidth,g);
    }

    normalize(traj, grp1, grp2, nBins, binwidth, g);

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
    for (int i = 0; i < nBins; i++) {
        oFS << ((double) i) * binwidth + start << "  " << g[i] << endl;
    }
    oFS.close();

    return 0;
}
