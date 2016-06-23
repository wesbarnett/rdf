/* 
 * James W. Barnett
 *
 * Calculates the radial distribution function. Requires libgmxcpp, which is
 * custom library written for reading in Gromacs files using OOP. Also needs a
 * configuration file to read in. See below on what is needed in such a file.
 *
 */


#include <chrono>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <string>
#include "gmxcpp/Trajectory.h"
#include "gmxcpp/Utils.h"
#include "omp.h"

const double FourThirdsPi = 4.0/3.0 * M_PI;
void doRdf(Trajectory &traj, string grp, double rexcl2, double end2, double binwidth, vector <double> &g);
void doRdf(Trajectory &traj, string grp1, string grp2, double rexcl2, double end2, double binwidth, vector <double> &g);
void doBinning(coordinates atomi, coordinates atomj, triclinicbox box, double boxvol, double binwidth, double rexcl2, double end2, vector <double> &g);
void normalize(Trajectory &traj, string grp1, string grp2, double binwidth, vector <double> &g);

using namespace std;

int main(int argc, char *argv[]) 
{


    if (argc != 2) 
    {
        cout << "Usage: " << endl;
        cout << "  " << argv[0] << " configfile" << endl;
        return -1;
    }

    chrono::time_point<std::chrono::system_clock> start_tp, end_tp;
    start_tp = std::chrono::system_clock::now(); 
    time_t start_time = chrono::system_clock::to_time_t(start_tp);

    const double start = 0.0;
    double binwidth;
    double end;
    double end2;
    double rexcl;
    double rexcl2;
    ifstream iFS;
    int nthreads;
    ofstream oFS;
    string configfile;
    string grp1;
    string grp2;
    string ndxfile;
    string outfile;
    string xtcfile;

    configfile = argv[1];

    #pragma omp parallel
    #pragma omp master
    nthreads = omp_get_num_threads();

    cout << "Reading " << configfile << "...";
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
    cout << "done." << endl;

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

    const int nBins = (end - start)/binwidth + 1;
    vector <double> g(nBins, 0.0);

    Trajectory traj(xtcfile, ndxfile);

    if (grp1 == grp2) 
    {
        doRdf(traj, grp1, rexcl2, end2, binwidth, g);
    } 
    else 
    {
        doRdf(traj, grp1, grp2, rexcl2, end2, binwidth, g);
    }
    normalize(traj, grp1, grp2, binwidth, g);

    end_tp = chrono::system_clock::now(); 
    chrono::duration<double> elapsed_seconds = end_tp-start_tp;
    time_t end_time = chrono::system_clock::to_time_t(end_tp);

    cout << "Writing output to " << outfile << "...";
    oFS.open(outfile.c_str());
    oFS << "# -------------------------------------" << endl;
    oFS << "# - RADIAL DISTRIBUTION - WES BARNETT -" << endl;
    oFS << "# -------------------------------------" << endl;
    oFS << "# Started computation:  " << ctime(&start_time);
    oFS << "# Finished computation: " << ctime(&end_time);
    oFS << "# OMP Threads:          " << nthreads << endl;
    oFS << "# Trajectory file:      " << xtcfile << endl;
    oFS << "# Index file:           " << ndxfile << endl;
    oFS << "# Output file:          " << outfile << endl;
    oFS << "# Group 1:              " << grp1 << endl;
    oFS << "# Group 2:              " << grp2 << endl;
    oFS << "# Exclusion distance:   " << rexcl << endl;
    oFS << "# Bin width:            " << binwidth << endl;
    oFS << "# Location of last bin: " << end << endl << endl;
    oFS << fixed << setprecision(6);
    for (int i = 0; i < nBins; i++) 
    {
        oFS << ((double) i) * binwidth + start << "  " << g[i] << endl;
    }
    oFS.close();
    cout << "done." << endl;

    return 0;
}

/* 
 * Calculates the radial distribution function when the two groups are the same
 * (e.g, OW-OW). It cycles through all frames stored in the Trajectory object.
 * In each frame it cycles through each pair of atoms, using a trick to make
 * sure we don't double count. Because normalization is done the same way as
 * when the groups are different, we add 2 * the boxvol each time. We add the
 * box volume instead of just counting because the box volume may not be
 * constant. If the box could be assured to be constant we would just put the
 * box volume in the normalization function, but instead we put it here. We also
 * exclude atoms closer than rexcl; this is so we don't count atoms on the same
 * molecule.
 */
void doRdf(Trajectory &traj, string grp, double rexcl2, double end2, double binwidth, vector <double> &g) 
{

    const int nFrames = traj.GetNFrames();
    const int nGrp = traj.GetNAtoms(grp);

    #pragma omp parallel
    {

        coordinates atomi;
        coordinates atomj;
        double boxvolX2;
        triclinicbox box;
        vector <double> g_thread(g.size(), 0.0);

        #pragma omp for schedule(static)
        for (int frame = 0; frame < nFrames; frame++) 
        {

            if (frame % 10 == 0) cout << "Frame: " << frame << endl;
            box = traj.GetBox(frame);
            boxvolX2 = 2.0 * volume(box);

            for (int i = 0; i < nGrp-1; i++) 
            {

                atomi = traj.GetXYZ(frame, grp, i);

                for (int j = i+1; j < nGrp; j++) 
                {

                    atomj = traj.GetXYZ(frame, grp, j);
                    doBinning(atomi, atomj, box, boxvolX2, binwidth, rexcl2, end2, g_thread);

                }

            }

        }

        #pragma omp critical
        {
            for (unsigned int i = 0; i < g.size(); i++)
            {
                g[i] += g_thread[i];
            }

        }

    }

    return;

}

/* 
 * Calculates the radial distribution function for two groups that are
 * different. For more details see the function above.
 */
void doRdf(Trajectory &traj, string grp1, string grp2, double rexcl2, double end2, double binwidth, vector <double> &g) 
{

    const int nFrames = traj.GetNFrames();
    const int nGrp1 = traj.GetNAtoms(grp1);
    const int nGrp2 = traj.GetNAtoms(grp2);

    #pragma omp parallel
    {
        coordinates atomi;
        coordinates atomj;
        double boxvol;
        triclinicbox box;
        vector <double> g_thread(g.size(), 0.0);

        #pragma omp for schedule(static)
        for (int frame = 0; frame < nFrames; frame++) 
        {

            if (frame % 10 == 0)
            {
                cout << "Frame: " << frame << endl;
            }

            box = traj.GetBox(frame);
            boxvol = volume(box);

            for (int i = 0; i < nGrp1; i++) 
            {

                atomi = traj.GetXYZ(frame, grp1, i);

                for (int j = 0; j < nGrp2; j++) 
                {

                    atomj = traj.GetXYZ(frame, grp2, j);
                    doBinning(atomi, atomj, box, boxvol, binwidth, rexcl2, end2, g_thread);
                }

            }

        }

        #pragma omp critical
        {
            for (unsigned int i = 0; i < g.size(); i++)
            {
                g[i] += g_thread[i];
            }

        }

    }
    

    return;
}


/* 
 * Does the actual binning - gets the squared distance, finds out if it is
 * within the noexcluded range and then adds the box volume to the bin
 * specified.
 */
void doBinning(coordinates atomi, coordinates atomj, triclinicbox box, double boxvol, double binwidth, double rexcl2, double end2, vector <double> &g) 
{

    double r2 = distance2(atomi,atomj,box);
    if (r2 > rexcl2 && r2 < end2) 
    {
        g[floor(sqrt(r2)/binwidth)] += boxvol;
    }

    return;

}

/* 
 * Normalizes the radial distribution function. Each bin has an associated
 * volume, which we divide by. The volume is spherical (but not a sphere, since
 * it is missing the inner part). The normalization factor also includes the
 * number of atoms in each group as well as the number of frames.
 */
void normalize(Trajectory &traj, string grp1, string grp2, double binwidth, vector <double> &g) 
{

    const int nFrames = traj.GetNFrames();
    const int nGrp1 = traj.GetNAtoms(grp1);
    const int nGrp2 = traj.GetNAtoms(grp2);

    double normFactor = (nGrp1-1) * nGrp2 * nFrames * FourThirdsPi * pow(binwidth,3);
    for (unsigned int i = 0; i < g.size(); i++) 
    {
        double r = (double) i;
        double binvol = pow(r+1.0,3) - pow(r,3);
        g[i] /= (binvol * normFactor);
    }

    return;

}

