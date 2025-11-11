#include <math.h>
#include <stdlib.h>

#include "../src/lbm.hh"

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
double angular_width;
double tip_to_tip_dist;
double row_to_row_dist;

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

int x_pos_first_row;
int z_pos_first_row;

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
    params.addParameter<double>(a1, "a1");
    params.addParameter<double>(b1, "b1");
    params.addParameter<double>(c1, "c1");
    params.addParameter<double>(angular_width, "angular_width");
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

    a2 = a1 + w;
    b2=  b1 + w;
    c2 = c1 + w;

    inter_ratchet_dist = tip_to_tip_dist - std::max(a2, b2);

    total_num_ratchets = Num_ratchets_per_row * Num_rows_of_ratchets;

    offset = 10;
    channel_length = ( Num_ratchets_per_row + 1 )*( std::max(a2,c2) + inter_ratchet_dist) + offset;
    reservoir_length = channel_length;
    lx = (int)3*reservoir_length;

    channel_height = 3*b2;
    reservoir_height = channel_height;
    ly = (int)3*channel_height;

    channel_width = ( Num_rows_of_ratchets )*( std::max(a2, b2) + row_to_row_dist );
    lz = channel_width;

    y0_val = channel_height;

    constraining_plate_height = reservoir_height + b2+10;
    constraining_plate_length = reservoir_length + 15;


    x_pos_first_row = reservoir_length + offset;
    z_pos_first_row = channel_width/2;

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
    for(int i = 0; i <  Num_rows_of_ratchets; i++)
    {

        for (int j = 0; j < Num_ratchets_per_row; j++)
        {
            center_coords_arr[i][j][0] = x_pos_first_row + j*tip_to_tip_dist;
            center_coords_arr[i][j][1] = z_pos_first_row + i*row_to_row_dist;
            

            ratchet_condition_arr[i][j][0] = pow(c1*b1, 2)*pow(-x+center_coords_arr[i][j][0],2) + pow(a1*c1,2)*pow(y-y0_val,2) + pow(a1*b1,2)*pow(z-center_coords_arr[i][j][1],2) > pow(a1*b1*c1,2);
            ratchet_condition_arr[i][j][1] = pow(c2*b2, 2)*pow(-x+center_coords_arr[i][j][0],2) + pow(a2*c2,2)*pow(y-y0_val,2) + pow(a2*b2,2)*pow(z-center_coords_arr[i][j][1],2) <= pow(a2*b2*c2,2);
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

    // Need bounceback conditions for reservoir (all walls need bounceback)
    // Need periodic boundary conditions in the x-direction for the reservoir
    
    if( ( x <= constraining_plate_length ) && ( y >= constraining_plate_height - 1 ) && ( y <= constraining_plate_height + 1 ) )
    {
        return 1; // Gives us the constraining plate on the left side of the domain
    }

    if( ( x >= lx - (reservoir_length - reservoir_length/2) ) && ( y >= constraining_plate_height - 1 ) && ( y <= constraining_plate_height + 1 ) )
    {
        return 2; // Gives us the constraining plate on the right side of the domain
    }

    // if(  ( x >= lx - (reservoir_length - reservoir_length/2 ) - 1 ) && ( x <= lx - (reservoir_length - reservoir_length/2 ) + 1 ) && ( y <= constraining_plate_height + 1)  )
    // {
    //     return 1; // Gives us a solid (vertical) wall going from the bottom of the domain to the constraining plate on the right side of the domain. 
    // }
    
    if( ( ( x <= reservoir_length - 1 ) && ( y <=2 ) ) || ( ( x <= reservoir_length - 1) && ( y >= ly - 2) ) )
    {
        return 1; // Gives us solid (horizontal) walls on the top and bottom of the left reservoir
    }

    if( ( ( x >= lx - (reservoir_length )  ) && ( y <= 2 ) ) || ( ( x >= lx - (reservoir_length-2) ) && ( y >= ly - 2) ) )
    {
        return 2; // Gives us solid (horizontal) walls on the top and bottom of the right reservoir
    }

    if( ( x >= reservoir_length - 1) && ( x <= reservoir_length + 1) && ( y >= 0 ) && ( y <= reservoir_height + 1 ) )
    {
        return 1; // Gives us a solid (vertical) wall between going from the bottom of the domain to the channel for the left part of the domain
    }

      if( ( x >= reservoir_length - 1) && ( x <= reservoir_length + 1) && ( y >= 2*(reservoir_height) -1 ) && ( y <= ly ) )
    {
        return 1; // Gives us a solid (vertical) wall between going from the top of the channel to the top of the domain for the left part of the domain
    }


    if( ( x >= lx - (reservoir_length+1) ) && ( x <= lx - (reservoir_length -1) ) && ( y >= 0 ) && ( y <= reservoir_height + 1 ) )
    {
        return 1; // Gives us a solid (vertical) wall between going from the bottom of the domain to the channel for the right part of the domain
    }

    if( ( x >= lx - (reservoir_length+1) ) && ( x <= lx - (reservoir_length-1) ) && ( y >= 2*(reservoir_height) -1 ) && ( y <= ly ) )
    {
        return 1; // Gives us a solid (vertical) wall between going from the top of the channel to the top of the domain for the right part of the domain
    }


    if( ( x > reservoir_length) && ( x <= reservoir_length + channel_length ) && ( y >= reservoir_height - 1) && ( y <= reservoir_height + 1 ) )
    {
        //std::cout << "here" << std::endl; 
        return 1; // Givves us solid walls on the bottom of the channel
    }

    if( ( x > reservoir_length) && ( x <= reservoir_length + channel_length ) && ( y >= reservoir_height + channel_height - 1) && ( y <= reservoir_height + channel_height + 1 ) )
    {
        // std::cout << "here" << std::endl; 
        return 1; // Gives us solid walls on the top of the channel
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
    
    if( ( xx < reservoir_length) && ( yy < constraining_plate_height) )
    {
        return 1;
    }

    if( (xx >= reservoir_length) && (xx <= reservoir_length + channel_length) && (yy > channel_height) && (yy <= constraining_plate_height) )
    {

        return 0.5 - 0.5*tanh(2*(xx - constraining_plate_length)/4.0);
    }

    if( ( xx >= right_reservoir - 10 ) && ( yy < constraining_plate_height -1) )
    {
        return 0.5 + 0.5*tanh(2*(xx - right_reservoir)/4.0);
    }

    // Sharp droplet
    // If (rr2 < radius * radius&&yy>postheight)
    //     return 1;
    // Else
    //     return 0;

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
    WellBalancedCH<0,Lattice, DefaultTraitWellBalancedCH2<0,Lattice>::template AddProcessor<std::tuple<GradientsMultiStencil<ChemicalPotential<0>,CentralQBounceBack,CentralXYZBounceBack>>>> binary;

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
