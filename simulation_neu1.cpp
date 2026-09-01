#include <iostream>
#include <fstream>
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

    virtual void debug_print(std::ostream &os) const {
        os << "SolidTube(x_start=" << x_start << ", x_end=" << x_end 
           << ", r_in=" << r_inner << ", r_out=" << r_outer << ")";
    }

    virtual void save(std::ostream &os) const {}
};

int main( int argc, char **argv ) {
    try {
        ibsimu.set_message_threshold( MSG_VERBOSE, 0 ); // Ausgaben für Sweep reduzieren

        double h = 0.0005; // Gitterabstand (0.5 mm)
        Int3D mesh_size( 1000, 80, 1 );
        Vec3D origin( 0.0, 0.0, 0.0 );

        double V_Repeller = 30.0;
        double V_L1       = -500.0;
        double V_L2       = -1000.0;
        double V_L4       = -500.0;
        double V_L5       = 0.0;

        double V_L3_start   = 0.0;
        double V_L3_min     = -3000.0;
        double schrittweite = 20.0; // 20 V Schritte

        double best_V_L3 = 0.0;
        int max_treffer_global = -1;

        ofstream fsweep("sweep_data.dat");
        fsweep << "# V_L3 Treffer\n";

        cout << "Starte vollstaendigen Sweep von " << V_L3_start << " V bis " << V_L3_min << " V..." << endl;

        for( double V_L3 = V_L3_start; V_L3 >= V_L3_min; V_L3 -= schrittweite ) {

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

            int treffer = 0;
            for( size_t i = 0; i < pdb.size(); ++i ) {
                const ParticleCyl &p = pdb.particle(i);
                if( p.traj_size() > 0 ) {
                    size_t last_idx = p.traj_size() - 1;
                    double end_x = p.traj(last_idx)[1];
                    double end_y = p.traj(last_idx)[3];

                    // Verschärfter Fokus-Radius: 2 mm anstatt 5 mm
                    if( end_x >= 0.380 && end_y <= 0.002 ) {
                        treffer++;
                    }
                }
            }

            fsweep << V_L3 << " " << treffer << "\n";
            cout << "V_L3 = " << V_L3 << " V | Treffer (r <= 2mm): " << treffer << endl;

            if( treffer > max_treffer_global ) {
                max_treffer_global = treffer;
                best_V_L3 = V_L3;

                // Speichere Trajektorien nur fuer den bisher besten Fokus
                ofstream fout( "trajectories_best.dat" );
                for( size_t i = 0; i < pdb.size(); ++i ) {
                    const ParticleCyl &p = pdb.particle(i);
                    for( size_t j = 0; j < p.traj_size(); ++j ) {
                        fout << p.traj(j)[1] << " " << p.traj(j)[3] << "\n";
                    }
                    fout << "\n";
                }
                fout.close();
            }
        }
        fsweep.close();

        cout << "\n==========================================" << endl;
        cout << "SWEEP ABGESCHLOSSEN!" << endl;
        cout << "Optimum V_L3: " << best_V_L3 << " V mit " << max_treffer_global << " Treffern im 2mm-Fokus." << endl;
        cout << "==========================================\n" << endl;

    } catch( Error &e ) {
        cerr << "IBSimu Fehler: " << e.get_error_message() << endl;
        return 1;
    }

    return 0;
}
