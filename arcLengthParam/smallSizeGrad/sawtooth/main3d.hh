#include <math.h>
#include <stdlib.h>
#include <functional>
#include "../../../src/Eigen/Dense"
#include "../../../src/lbm.hh"


// Quadrature (simple Simpson's rule for demo)
double simpsonIntegral(std::function<double(double)> f, double a, double b, int n = 200) {
    if (n % 2 != 0) n++; // Simpson's requires even number of intervals
    double h = (b - a) / n;
    double s = f(a) + f(b);
    for (int i = 1; i < n; i++) {
        double x = a + i*h;
        s += (i % 2 == 0 ? 2 : 4) * f(x);
    }
    return s * h / 3.0;
}

// System F(x)
Eigen::Vector3d F(const Eigen::Vector3d& x, double alpha, double S1, double S2) {
    double a = x(0), c = x(1), theta = x(2);
    Eigen::Vector3d f;

    // f1
    f(0) = std::atan(c/a) - alpha;

    // f2: integral from pi/2 to 0 -> negative of 0->pi/2
    auto integrand = [&](double z) {
        double val = 1.0 + ( (a*a)/(c*c) - 1.0 ) * std::sin(z)*std::sin(z);
        return std::sqrt(val);
    };
    double I = simpsonIntegral(integrand, 0.0, M_PI/2);
    f(1) = -c * (-I) - S1; // = c*I - S1

    // f3
    f(2) = a*theta - S2;

    return f;
}

// Numerical Jacobian via finite differences
Eigen::Matrix3d numericalJacobian(
    const Eigen::Vector3d& x,
    std::function<Eigen::Vector3d(const Eigen::Vector3d&)> Ffun)
{
    double eps = 1e-6;
    Eigen::Matrix3d J;
    Eigen::Vector3d f0 = Ffun(x);
    for (int j = 0; j < 3; j++) {
        Eigen::Vector3d xh = x;
        xh(j) += eps;
        Eigen::Vector3d fh = Ffun(xh);
        J.col(j) = (fh - f0) / eps;
    }
    return J;
}

int lx;             // Simulation size in x direction
int ly;             // Simulation size in y direction
int lz;               // Simulation size in z direction
int timesteps ;     // Total number of iterations
int saveInterval;  // How frequently to save order parameter, velocity etc
double radius;     // Droplet radius
double theta;      // Contact angle
double A;         // Parameter for free energy functional (depends on surface tension and interface width, you can
                          // Ignore for now)
double kappa; // Parameter for free energy functional (depends on surface tension and interface width, you can ignore for now)
double surfacetension;  // Surface tension in lattice units
double interfacewidth;      // Width of the diffuse interface between fluids
double posx;        // Droplet position in the x direction
double posy;        // Droplet position in the y direction
double posz;        // Droplet position in the z direction, uncomment for full 3D
double dens1;         // Droplet density
double dens2;         // Air density
double bodyforcex;    // Magnitude of body force acting on the droplet in the x direction
std::string datadir = "data/"; // Data directory
double postfraction; // Approximate fraction of area in the x direction taken up by posts (approximate because we have a finite number of lattice points)
double postfractionZ; // Approximate fraction of area in the x direction taken up by posts
int equilibriumtimesteps; // Number of timesteps until body force is applied
int postheight;      // Height of the posts
int nbPost;          // Number of posts in the x direction
int nbPostZ;          // Number of posts in the z direction
double tau1 = 1;          // Relaxation time of droplet
double tau2 = 1;          // Relaxation time of droplet

double w;
double a1;
double a2;
double b1;
double b2;
double c1;
double c2;
double a1p;
double a2p;
double b1p;
double b2p;
double c1p;
double c2p;
double angular_width;
double tip_to_tip_dist;
double row_to_row_dist;

double tilt_angle;
double S1;
double S2;

int Num_ratchets_per_row;
int Num_rows_of_ratchets;

double left_liquid_level;
double right_liquid_level;

double inter_ratchet_dist;
int total_num_ratchets;

double reservoir_length;
double reservoir_height;

int channel_width;
double channel_height;
double channel_length;
double y0_val = channel_height;

double constraining_plate_height;
double constraining_plate_length;

int offset;

int x_pos_first_row_front;
int x_pos_first_row_behind;
int x_pos_second_row;
int z_pos_first_row;

int apertureHeight;
int apertureWidth;
int reservoirHeight;

int x_c;
int z_c;
int apWidth2;

double pi = 3.141592653589;

// Class that will handle input parameters
InputParameters params;

// Lattice classes, contain information about the size of the simulation and the parallelisation.
#ifdef MPIPARALLEL //If you have -DM
using Lattice = LatticePropertiesRuntime<ParallelX<2> /*Class to handle parallelisation*/, 3 /*Number of physical dimensions*/>;
#else
using Lattice = LatticePropertiesRuntime<NoParallel, 3>;
#endif

// Function to initialise the parameters from an input file
void initParams(std::string inputfile) {
    //params.addParameter<int>(lx, "lx");
    //params.addParameter<int>(ly, "ly");
    //params.addParameter<int>(lz, "lz");
    params.addParameter<int>(timesteps, "timesteps");
    params.addParameter<int>(saveInterval, "saveInterval");
    params.addParameter<double>(radius, "radius");
    params.addParameter<double>(theta, "theta");
    params.addParameter<double>(A, "A");
    params.addParameter<double>(kappa, "kappa");
    params.addParameter<double>(surfacetension, "surfacetension");
    params.addParameter<double>(interfacewidth, "interfacewidth");
    params.addParameter<double>(posx, "posx");
    params.addParameter<double>(posy, "posy");
    params.addParameter<double>(posz, "posz"); // Uncomment for 3D
    params.addParameter<double>(dens1, "dens1");
    params.addParameter<double>(dens2, "dens2");
    params.addParameter<double>(tau1, "tau1");
    params.addParameter<double>(tau2, "tau2");
    params.addParameter<std::string>(datadir, "datadir");
    params.addParameter<double>(postfraction, "postfraction");
    params.addParameter<double>(postfractionZ, "postfractionZ");
    params.addParameter<double>(bodyforcex, "bodyforcex");
    params.addParameter<int>(postheight, "postheight");
    params.addParameter<int>(nbPost, "nbPost");
    params.addParameter<int>(nbPostZ, "nbPostZ");
    params.addParameter<int>(equilibriumtimesteps, "equilibriumtimesteps");
    params.addParameter<double>(w, "w");
    params.addParameter<double>(tilt_angle, "tilt_angle");
    params.addParameter<double>(S1, "S1");
    params.addParameter<double>(S2, "S2");
    params.addParameter<double>(tip_to_tip_dist, "tip_to_tip_dist");
    params.addParameter<double>(row_to_row_dist, "row_to_row_dist");
    params.addParameter<int>(Num_ratchets_per_row, "Num_ratchets_per_row");
    params.addParameter<int>(Num_rows_of_ratchets, "Num_rows_of_ratchets");

    /*
    If you want to add a parameter here, follow the format above
    params.addParameter< *parameter type* >( *parameter name in this file*, *parameter name in the input file* );
    */

    // Read the input file and initialise the parameters
    params.readInput(inputfile);

    tilt_angle = tilt_angle * 3.14159265 / 180.;

    auto Ffun = [&](const Eigen::Vector3d& x) {
        return F(x, tilt_angle, S1, S2);
    };

    Eigen::Vector3d x(1.0, 1.0, 1.0); // initial guess: (a,c,theta)

    for (int iter = 0; iter < 20; ++iter) {
        Eigen::Vector3d f = Ffun(x);
        if (f.norm() < 1e-10) break;
        Eigen::Matrix3d J = numericalJacobian(x, Ffun);
        Eigen::Vector3d dx = J.fullPivLu().solve(f);
        x = x - dx;
        std::cout << "Iter " << iter
                  << ": x = " << x.transpose()
                  << " |f| = " << f.norm() << "\n";
    }


    std::cout << "\n\n angular_width = " << x[2] << std::endl;
    a1 = x[0];
    b1 = a1;
    c1 = x[1];
    angular_width = x[2]*180/3.14159;

    a2 = a1 + w;
    b2=  b1 + w;
    c2 = c1 + w;

    inter_ratchet_dist = tip_to_tip_dist;

    total_num_ratchets = Num_ratchets_per_row * Num_rows_of_ratchets;

    offset = 15;
    channel_length = ( Num_ratchets_per_row - 3 )*( std::max(a2,c2) + inter_ratchet_dist );
    lx = (int)1*channel_length;
    //lx = 250;
    std::cout << "lx = " << lx << std::endl;

    ly = (int)7*b2;
    //ly = ((ly + 99) / 100) * 100;
    std::cout << "ly = " << ly << std::endl;

    lz = std::max(a2, b2) + 2.0*row_to_row_dist;
    //lz = ((lz + 99) / 100) * 100;

    std::cout << "lz = " << lz << std::endl;

    x_c = (int)(lx/2);
    z_c = (int)(lz/2);
    

    apertureHeight = 5;
    apertureWidth = 2*apertureHeight;
    apWidth2 = (int)(apertureWidth/2);

    reservoirHeight = int(0.3 * ly);

    y0_val = reservoirHeight + apertureHeight;

    x_pos_first_row_front = a2;
    x_pos_second_row = x_pos_first_row_front + 20;
    x_pos_first_row_behind = x_c - a2/2;
    z_pos_first_row = z_c - (int)(row_to_row_dist/2);

    // Initialise free energy parameters from the surface tension and interface width
    A = 12 * surfacetension / interfacewidth;
    kappa = pow(interfacewidth, 2) / 8.0 * A;

    // Initialise the lattice class with the simulation size
    Lattice::init(lx, ly, lz);

}


// double stepFunction(int xx, int yy, int zz, double y0_val)
// {

//     // The step function works as follows: if you're in the solid, return 0. If you're in the liquid, return 1. 
//     double x = (double)xx;
//     double y = (double)yy;
//     double z = (double)zz;

//     std::vector<std::vector<std::array<bool, 4>>> ratchet_condition_arr(
//     Num_rows_of_ratchets, std::vector<std::array<bool,4>>(Num_ratchets_per_row));

//     std::vector<std::vector<std::array<double,2>>> center_coords_arr(
//         Num_rows_of_ratchets, std::vector<std::array<double,2>>(Num_ratchets_per_row));

//     for(int i = 0; i <  Num_rows_of_ratchets; i++)
//     {

//         for (int j = 0; j < Num_ratchets_per_row; j++)
//         {
//             center_coords_arr[i][j][0] = x_pos_first_row + j*tip_to_tip_dist;
//             center_coords_arr[i][j][1] = z_pos_first_row + i*row_to_row_dist;
            

//             ratchet_condition_arr[i][j][0] = pow(c1, 2)*( pow(x-center_coords_arr[i][j][0],2) + pow(z-center_coords_arr[i][j][1],2) ) + pow(a1,2)*pow(y-y0_val,2) > pow(a1*c1,2);
//             ratchet_condition_arr[i][j][1] = pow(c2, 2)*( pow(x-center_coords_arr[i][j][0],2) + pow(z-center_coords_arr[i][j][1],2) ) + pow(a2,2)*pow(y-y0_val,2) <= pow(a2*c2,2);
//             ratchet_condition_arr[i][j][2] = x - center_coords_arr[i][j][0] > 0;
//             ratchet_condition_arr[i][j][3] = atan( (x-center_coords_arr[i][j][0]) / abs(z-center_coords_arr[i][j][1]) ) > (pi/2 - angular_width*pi/180);
//         }
//     }

//     for (int i = 0; i < Num_rows_of_ratchets; i++)
//     {
//         for (int j = 0; j < Num_ratchets_per_row; j++)
//         if( ( ratchet_condition_arr[i][j][0] && ratchet_condition_arr[i][j][1] && ratchet_condition_arr[i][j][2] && ratchet_condition_arr[i][j][3] ) )
//         {

//             return 0;
//         }
//     }

//     if( y <= 2)
//     {
//         return 0;
//     }

//     return 1;
    

// }

// Function used to define the solid geometry
int initBoundary(const int k) {
    // x coordinate of the lattice node k
    int xx = computeXGlobal<Lattice>(k);
    // y coordinate of the lattice node k
    int yy = computeY(ly, lz, k);
    // z coordinate of the lattice node k
    int zz = computeZ(ly, lz, k);

    double x = (double)xx;
    double y = (double)yy;
    double z = (double)zz;

    std::vector<std::vector<std::array<bool, 5>>> ratchet_condition_arr(
    Num_rows_of_ratchets, std::vector<std::array<bool,5>>(Num_ratchets_per_row));

    std::vector<std::vector<std::array<double,2>>> center_coords_arr(
        Num_rows_of_ratchets, std::vector<std::array<double,2>>(Num_ratchets_per_row));



    // These two for loops set up the logical conditions describing the ratchet geometry.
    // unfortunately, these loops have to be invoked for every point in the domain

    double size_fac;
    int anchor_pt;
    for(int i = 0; i <  Num_rows_of_ratchets; i++)
    {

        for (int j = 0; j < Num_ratchets_per_row; j++)
        {
            size_fac =  1 + 0.6 * double( ( j % 8 ) ) / 9 ;
            //std::cout << "size_factor = " << size_fac << std::endl; 
            a1p = size_fac*a1;
            a2p = size_fac*a2;
            b1p = size_fac*b1;
            b2p = size_fac*b2;
            c1p = size_fac*c1;
            c2p = size_fac*c2;

            if(i == 0)
            {
                int anchor_pt = x_pos_first_row_front + (j)*(tip_to_tip_dist);
                center_coords_arr[i][j][0] = anchor_pt;
                //std::cout << "x_c = " << center_coords_arr[i][j][0] << std::endl;
            }
            
            if(i == 1)
            {
                int anchor_pt = x_pos_second_row + (j)*(tip_to_tip_dist);
                center_coords_arr[i][j][0] = anchor_pt;
            }

            
          


            center_coords_arr[i][j][1] = z_pos_first_row + i*row_to_row_dist;

            ratchet_condition_arr[i][j][0] = pow(c1p*b1p, 2)*pow(-x+center_coords_arr[i][j][0],2) + pow(a1p*c1p,2)*pow(y-y0_val,2) + pow(a1p*b1p,2)*pow(z-center_coords_arr[i][j][1],2) > pow(a1p*b1p*c1p,2);
            ratchet_condition_arr[i][j][1] = pow(c2p*b2p, 2)*pow(-x+center_coords_arr[i][j][0],2) + pow(a2p*c2p,2)*pow(y-y0_val,2) + pow(a2p*b2p,2)*pow(z-center_coords_arr[i][j][1],2) <= pow(a2p*b2p*c2p,2);
            ratchet_condition_arr[i][j][2] = x - center_coords_arr[i][j][0] > 0;
            ratchet_condition_arr[i][j][3] = y - y0_val > 0;
            ratchet_condition_arr[i][j][4] = atan( (x-center_coords_arr[i][j][0]) / abs(z-center_coords_arr[i][j][1]) ) > (pi/2 - angular_width*pi/180);
        }
    }

    // This loop checks if any of the ratchet conditions are satisfied for a given point in the domain. If any of the criteria are satisfied, 
    // we return 1, indicating that the given point is inside a ratchet. 
    for (int i = 0; i < Num_rows_of_ratchets; i++)
    {
        for (int j = 0; j < Num_ratchets_per_row; j++)
        if( ratchet_condition_arr[i][j][0] && ratchet_condition_arr[i][j][1] && ratchet_condition_arr[i][j][2] && ratchet_condition_arr[i][j][3] && ratchet_condition_arr[i][j][4] )
        {
            return 1;
        }
    }


    bool cond1 = ( xx < x_c + apWidth2 ) && ( xx > x_c - apWidth2 ) && ( zz < z_c + apWidth2 ) && ( zz > z_c - apWidth2 );

    if( ( yy >= reservoirHeight ) && ( yy <= reservoirHeight + apertureHeight - 2 ) && !cond1 )
    {
        return 2;
    }

    if( ( yy >= reservoirHeight + apertureHeight -1 ) && ( yy <= reservoirHeight + apertureHeight ) && !cond1 )
    {
        return 1;

    }
    

    return 0;
}

//Initialises the fluid with a bulk value of 1 for the droplet and a bulk value of 0 for air
double initFluid(const int k) {
    // x coordinate of the lattice node k
    int xx = computeXGlobal<Lattice>(k);
    // y coordinate of the lattice node k
    int yy = computeY(ly, lz, k);
    // z coordinate of the lattice node k
    int zz = computeZ(ly, lz, k); // Uncomment for 3D

    double right_reservoir = lx - (reservoir_length - reservoir_length/2);

    // Radial distance from the centre of the droplet at (posx, posy)
    //double rr2 = (xx - posx) * (xx - posx) + (yy - posy) * (yy - posy);
    double rr2 = (xx - posx) * (xx - posx) + (yy - posy) * (yy - posy) + (zz - posz) * (zz - posz); // Switch these for 3D

    // Smooth droplet

    int x_c = (int)(lx/2);
    int z_c = (int)(lz/2);
    int apWidth2 = (int)(apertureWidth/2);

    bool cond1 = ( xx < x_c + apWidth2 ) && ( xx > x_c - apWidth2 ) && ( zz < z_c + apWidth2 ) && ( zz > z_c - apWidth2 );
    
    if( ( yy <= reservoirHeight ) )
    {
        return 1;
    }

    if( ( yy <= reservoirHeight + apertureHeight + 1 ) && cond1 )
    {
        return 1;
    }

    int posx = x_c;
    int posy = reservoirHeight + apertureHeight;
    int posz = z_c;
    double radius = 0.85*2*b2;

    double dist = std::max( abs(yy - posy), abs(xx - posx) );



    if( ( yy > reservoirHeight + apertureHeight ) ) //&& ( zz >= z_c - apWidth2 ) && ( zz <= z_c + apWidth2 ) )
    {
        //return 0.5 - 0.5*tanh(2*(yy - ( reservoirHeight + apertureHeight + 3 ) )/4.0);
        if(initBoundary(k) == 1)
        {
            return 0;
        }
        else
        {
            return (0.5 - 0.5 * tanh(2 * (dist - radius) / interfacewidth ));
        }
    }

    return 0;
}

//
///////// Simulation details
//
using densitygradients =
    GradientsMultiStencil<Density<>, CentralXYZBounceBack, CentralQBounceBack>;
using velocitygradients =
    GradientsDirectional<Velocity<>, CentralXYZBounceBackDirectional>;

template<class TLattice>
using DefaultTraitPressureWellBalancedNew = typename DefaultTrait<TLattice,2>::template SetBoundary<BounceBack, FreeSlip>
                                                               :: template SetProcessor<std::tuple<MirrorBoundary<Density<>>,MirrorBoundary<Velocity<>,Lattice::NDIM>>,std::tuple<densitygradients,velocitygradients,ViscousStressCalculator>>
                                                               :: template SetForce< PressureForceWellBalanced<WellBalancedForce<NoTauDependence>,ChemicalForceBinaryMu>, PressureForceWellBalancedTau<SimpleForcingQ<SplitGuoPrefactor>,ChemicalForceBinaryMu>, BodyForce<> >;


// Function to create the pressure LBM model (solves navier stokes)
auto initPressure() {

    // PressureLee is the model class, and we give the lattice and traits as template parameters.
    FlowFieldPressureWellBalanced3<Lattice, typename DefaultTraitPressureWellBalancedNew<Lattice>::template SetCollisionOperator<MRT>> pressure;

    // Boundary ids to apply the LBM model
    pressure.setCollideID({0});
    //pressure.template getBoundary<EquilibriumPL>().setNodeID(11);
    // Set relaxation times of each component
    pressure.setTaus(tau1,tau2);
    pressure.setDensities(dens1,dens2);

    pressure.template getBoundary<BounceBack>().setNodeID({1,2});

    // Apply the mirror boundary condition on all nodes with id 2
    pressure.template getBoundary<FreeSlip>().setNodeID(4);


    // Set magnitude of bodyforce and component to apply it to, will be zero until AfterEquilibration is called
    pressure.template getForce<BodyForce<>>().setForce({0, 0}, // 0 in the x direction, 0 in the y direction
                                                        0); // Act on 0 component (droplet)

    pressure.template getProcessor<densitygradients>().setBoundaryID({1,2});
    pressure.template getProcessor<velocitygradients>().setBoundaryID({1,2});
    pressure.template getProcessor<MirrorBoundary<Density<>>>().setNodeID(4);
    pressure.template getProcessor<MirrorBoundary<Velocity<>,Lattice::NDIM>>().setNodeID(4);                         

    // Return the model so it can be used in main.cc
    return pressure;
}
using orderparamgradients = GradientsMultiStencil<OrderParameter<>, CentralXYZBounceBack, CentralQBounceBack, LaplacianCentralWetting>;


template<int N, int Nmax, class TLattice>
using DefaultTraitWellBalancedCH2 = typename DefaultTrait<TLattice,Nmax>::template SetBoundary<BounceBack, FreeSlip>
                                                                ::template SetProcessor<std::tuple<LinearWetting,MirrorBoundary<OrderParameter<>>>,std::tuple<orderparamgradients, ChemicalPotentialCalculatorBinaryLee>,std::tuple<MirrorBoundary<ChemicalPotential<>>>> // Note the inclusion of a class to calculate the chemical potential
                                                               :: template SetForce< NoSource<N> >;

// Function to create the order parameter LBM model (solves cahn hilliard)
auto initBinary() {

    // BinaryLee is the model class, and we give the lattice and traits as template parameters.
    //BinaryWellBalanced<Lattice, typename DefaultTraitBinaryWellBalanced<Lattice>::template SetStencil<D2Q5>> binary;
    WellBalancedCH<0,2,Lattice, DefaultTraitWellBalancedCH2<0,2,Lattice>::template AddProcessor<std::tuple<GradientsMultiStencil<ChemicalPotential<0>,CentralQBounceBack,CentralXYZBounceBack>>>> binary;

    // Boundary ids to apply the LBM model
    binary.setCollideID({0,11});

    // Chemical potential calculation needs the free energy parameters
    binary.template getProcessor<ChemicalPotentialCalculatorBinaryLee>().setA(A);
    binary.template getProcessor<ChemicalPotentialCalculatorBinaryLee>().setKappa(kappa);
    
    binary.template getProcessor<ChemicalPotentialCalculatorBinaryLee>().setOmega(0.01);

    binary.setMij({-1,0});    
    
    binary.template getBoundary<BounceBack>().setNodeID({1,2});

    // Apply the mirror boundary condition on all nodes with id 2
    binary.template getBoundary<FreeSlip>().setNodeID(4);

    binary.template getProcessor<orderparamgradients>().setBoundaryID({1,2});
    binary.template getProcessor<GradientsMultiStencil<ChemicalPotential<0>,CentralQBounceBack,CentralXYZBounceBack>>().setBoundaryID({1,2});

    binary.template getProcessor<MirrorBoundary<OrderParameter<>>>().setNodeID(4);
    binary.template getProcessor<MirrorBoundary<ChemicalPotential<>>>().setNodeID(4);

    double wettingprefactor1; 
    double wettingprefactor2;
    
    wettingprefactor1 = -2 * cos(30 * M_PI / 180.0) * sqrt(2 * A / kappa);
    wettingprefactor2 = -2 * cos(90 * M_PI / 180.0) * sqrt(2 * A / kappa);
    
    // wettingprefactor1 corresponds to boundary label 1, ie where the ratchets are 
    // wettingprefactor2 corresponds to boundary label 2, ie the reservoir, which should /always/ have \theta = 90.
    
    binary.template getProcessor<orderparamgradients>().setWettingPrefactor({1,2}, {wettingprefactor1, wettingprefactor2});
    // Stabilisation parameter to stop the order parameter going below 0


    return binary;
}


//Will apply the body force after the equilibrium timesteps
template <typename T>
void AfterEquilibration(int eqsteps, T& model) {
    if (eqsteps == equilibriumtimesteps) model.template getForce<BodyForce<>>().setForce({bodyforcex, 0}, // bodyforcex in the x direction, 0 in the y direction
                                                                                         0); // Act on 0 component (droplet)
}
