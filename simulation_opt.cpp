#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "geometry.hpp"
#include "solid.hpp"
#include "epot_efield.hpp"
#include "epot_bicgstabsolver.hpp"
#include "meshvectorfield.hpp"
#include "meshscalarfield.hpp"
#include "particledatabase.hpp"
#include "ibsimu.hpp"
#include "error.hpp"

using namespace std;

class SolidTube : public Solid {
    double x_start, x_end, r_inner, r_outer;
public:
    SolidTube(const Vec3D &p1, const Vec3D &p2, double r_in, double r_out) {
        x_start = p1[0];
        x_end   = p2[0];
        r_inner = r_in;
        r_outer = r_out;
    }
    
    virtual bool inside(const Vec3D &pos) const {
        return (pos[0] >= x_start && pos[0] <= x_end && 
                pos[1] >= r_inner && pos[1] <= r_outer);
    }

    virtual void debug_print(std::ostream &os) const {}
    virtual void save(std::ostream &os) const {}
};

int main( int argc, char **argv ) {
    if( argc < 5 ) {
        cerr << "Nutzung: " << argv[0] << " V_L1 V_L2 V_L3 V_L4 [traj_file]" << endl;
        return 1;
    }

    double V_L1 = atof(argv[1]);
    double V_L2 = atof(argv[2]);
    double V_L3 = atof(argv[3]);
    double V_L4 = atof(argv[4]);
    double V_Repeller = 30.0;
    double V_L5       = 0.0;

    try {
        ibsimu.set_message_threshold( MSG_VERBOSE, 0 );

        double h = 0.0005;
        Int3D mesh_size( 1000, 80, 1 );
        Vec3D origin( 0.0, 0.0, 0.0 );

        Geometry geom( MODE_CYL, mesh_size, origin, h );

        geom.set_boundary( 1, Bound(BOUND_NEUMANN, 0.0) );
        geom.set_boundary( 2, Bound(BOUND_NEUMANN, 0.0) );
        geom.set_boundary( 3, Bound(BOUND_NEUMANN, 0.0) );
        geom.set_boundary( 4, Bound(BOUND_DIRICHLET, 0.0) );
        geom.set_boundary( 5, Bound(BOUND_NEUMANN, 0.0) );
        geom.set_boundary( 6, Bound(BOUND_NEUMANN, 0.0) );

        geom.set_solid( 7,  new SolidTube(Vec3D(0.013, 0, 0), Vec3D(0.018, 0, 0), 0.0050, 0.035) );
        geom.set_solid( 8,  new SolidTube(Vec3D(0.018, 0, 0), Vec3D(0.020, 0, 0), 0.0075, 0.035) );
        geom.set_solid( 9,  new SolidTube(Vec3D(0.030, 0, 0), Vec3D(0.094, 0, 0), 0.0125, 0.035) );
        geom.set_solid( 10, new SolidTube(Vec3D(0.096, 0, 0), Vec3D(0.099, 0, 0), 0.0075, 0.035) );
        geom.set_solid( 11, new SolidTube(Vec3D(0.101, 0, 0), Vec3D(0.341, 0, 0), 0.0130, 0.035) );
        geom.set_solid( 12, new SolidTube(Vec3D(0.361, 0, 0), Vec3D(0.363, 0, 0), 0.0150, 0.035) );

        geom.set_boundary( 7,  Bound(BOUND_DIRICHLET, V_Repeller) );
        geom.set_boundary( 8,  Bound(BOUND_DIRICHLET, V_L1) );
        geom.set_boundary( 9,  Bound(BOUND_DIRICHLET, V_L2) );
        geom.set_boundary( 10, Bound(BOUND_DIRICHLET, V_L3) );
        geom.set_boundary( 11, Bound(BOUND_DIRICHLET, V_L4) );
        geom.set_boundary( 12, Bound(BOUND_DIRICHLET, V_L5) );

        geom.build_mesh();

        MeshScalarField scharge( geom );
        EpotBiCGSTABSolver solver( geom );
        EpotField epot( geom );

        solver.solve( epot, scharge );
        EpotEfield efield( epot );

        ParticleDataBaseCyl pdb( geom );
        pdb.add_2d_beam_with_energy( 100, 1.0e-6, 1.0, 40.0, 100.0, 0.0, 0.0, 
                                     0.001, 0.0, 0.001, 0.003 );

        MeshVectorField bfield;
        pdb.iterate_trajectories( scharge, efield, bfield );

        double sum_y_sq = 0.0;
        int treffer = 0;
        for( size_t i = 0; i < pdb.size(); ++i ) {
            const ParticleCyl &p = pdb.particle(i);
            if( p.traj_size() > 0 ) {
                size_t last_idx = p.traj_size() - 1;
                double end_x = p.traj(last_idx)[1];
                double end_y = p.traj(last_idx)[3];

                if( end_x >= 0.380 ) {
                    sum_y_sq += end_y * end_y;
                    treffer++;
                }
            }
        }

        if( treffer < 30 ) {
            cout << 999.0 << endl;
            return 0;
        }

        double rms_r_mm = sqrt(sum_y_sq / treffer) * 1000.0;

        // Schreiben der Daten bei Übergabe eines 5. Argumentes (Optimum)
        if( argc >= 6 ) {
            string filename = argv[5];
            ofstream fout( filename.c_str() );
            for( size_t i = 0; i < pdb.size(); ++i ) {
                const ParticleCyl &p = pdb.particle(i);
                for( size_t j = 0; j < p.traj_size(); ++j ) {
                    fout << p.traj(j)[1] << " " << p.traj(j)[3] << "\n";
                }
                fout << "\n";
            }
            fout.close();

            // Korrigierter Aufruf: epot verlangt Raumkoordinaten Vec3D(x, y, 0.0)
            ofstream fpot( "epot.dat" );
            for( int i = 0; i < mesh_size[0]; ++i ) {
                for( int j = 0; j < mesh_size[1]; ++j ) {
                    double x = origin[0] + i * h;
                    double y = origin[1] + j * h;
                    fpot << x << " " << y << " " << epot( Vec3D(x, y, 0.0) ) << "\n";
                }
                fpot << "\n";
            }
            fpot.close();
        }

        cout << rms_r_mm << endl;

    } catch( Error &e ) {
        cout << 999.0 << endl;
        return 1;
    }

    return 0;
}
