/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Klothakis, A. and Nikolos, I., “Modeling of Rarefied Hypersonic Flows Using the Massively
 *        Parallel DSMC Kernel “SPARTA”,” in 8th GRACM International Congress on Computational Mechanics,
 *        Volos, Greece, July 2015.
 *      Plimpton, S. and Gallis, M., SPARTA Users Manual, Sandia National Laboratories, United States
 *        Department of Energy, July 2017.
 *      Liechty, D., “Aeroheating Analysis for the Mars Reconnaissance Orbiter with Comparison to Flight Data,”
 *        Journal of Spacecraft and Rockets, vol. 44, no. 6, pp. 1226–1231, 2007.
 *
 */

#ifndef TUDAT_RAREFIED_FLOW_ANALYSIS_H
#define TUDAT_RAREFIED_FLOW_ANALYSIS_H

#include <map>
#include <string>
#include <vector>

#include <boost/array.hpp>
#include <boost/multi_array.hpp>
#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/Aerodynamics/aerodynamicCoefficientGenerator.h"
#include "Tudat/Astrodynamics/Aerodynamics/tabulatedAtmosphere.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"

#include "Tudat/Basics/basicTypedefs.h"

namespace tudat
{

namespace aerodynamics
{

//! Enumaration of key values for map of atmospheric conditions.
enum AtmosphericConditionVariables
{
    density_index = 0,
    pressure_index = 1,
    temperature_index = 2,
    speed_of_sound_index = 3,
    number_density_index = 4
};

//! Returns default values of altitude for use in RarefiedFlowAnalysis.
/*!
 *  Returns default values of altitude for use in RarefiedFlowAnalysis.
 */
std::vector< double > getDefaultRarefiedFlowAltitudePoints(
        const std::string& targetPlanet = "Earth" );

//! Returns default values of Mach number for use in RarefiedFlowAnalysis.
/*!
 *  Returns default values of Mach number for use in RarefiedFlowAnalysis.
 */
std::vector< double > getDefaultRarefiedFlowMachPoints(
        const std::string& machRegime = "Full" );

//! Returns default values of angle of attack for use in RarefiedFlowAnalysis.
/*!
 *  Returns default values of angle of attack for use in RarefiedFlowAnalysis.
 */
std::vector< double > getDefaultRarefiedFlowAngleOfAttackPoints(
        const std::string& angleOfAttackRegime = "Reduced" );

//! Class for aerodynamic analysis in rarefied flow using the SPARTA DSMC method.
/*!
 * Class for aerodynamic analysis in rarefied flow using the SPARTA DSMC method.
 * This method uses a Monte Carlo simulation, to determine the pressure and shear forces
 * acting on each element of the vehicle. These values are output by default every 200 time
 * steps, and are used to compute the average pressure and friction coefficients on each surface
 * element, which are then translated to aerodynamic coefficients for the whole surface. One can
 * find a description of the SPARTA software as Reference [1,2], where the second reference is the
 * official user manual. The user should also pay careful attention to the requirements for the geometry
 * of the vehicle to be analyzed.
 */
class RarefiedFlowAnalysis: public AerodynamicCoefficientGenerator< 3, 6 >
{
public:

    //! Default constructor.
    /*!
     *  Default constructor of class, specified vehicle geometry, discretization properties,
     *  independent variable ranges, reference values and local inclination methods that are
     *  to be used
     *  \param SPARTAExecutable Path to executable for SPARTA simulation.
     *  \param dataPointsOfIndependentVariables Vector of vectors, with each subvector containing
     *          the data points of each of the independent variables for the coefficient generation.
     *          The physical meaning of each of the three independent variables is: 0 = altitude,
     *          1 = Mach number, 1 = angle of attack.
     *          Each of the subvectors must be sorted in ascending order.
     *  \param simulationGases String of gases making up the atmosphere of the planet, to be used for
     *          the simulation.
     *  \param atmosphereModel Pointer to the atmosphere model of the planet. Should provide information on
     *          density, pressure, temperature, gas constant and specific heat ratio.
     *  \param geometryFileUser Path to the file describing the geometry of the vehicle, where
     *          the surface elements are discretized as triangles.
     *  \param referenceArea Reference area used to non-dimensionalize aerodynamic forces
     *          and moments.
     *  \param referenceLength Reference length used to non-dimensionalize aerodynamic moments.
     *  \param referenceAxis Index of main axis of the vehicle (i.e., axis opposite in direction to
     *          incoming flow, when angles of attack and sideslip are zero).
     *  \param momentReferencePoint Reference point wrt which aerodynamic moments are calculated.
     *  \param gridSpacing Grid size for simulation environment, used to define the size of each cell, and
     *          the number of cells in the environment.
     *  \param simulatedParticlesPerCell Number of simulated particles per cell.
     *  \param wallTemperature Temperature of the surface of the vehicle (default value is 300 K [3]).
     *  \param accommodationCoefficient Accommodation coefficient of the surface of the vehicle. This
     *          value indicates the degree of diffusivity during molecular-surface collisions (default value
     *          is 1.0, i.e., diffuse reflection).
     */
    RarefiedFlowAnalysis(
            const std::string& SPARTAExecutable,
            const std::vector< std::vector< double > >& dataPointsOfIndependentVariables,
            boost::shared_ptr< TabulatedAtmosphere > atmosphereModel,
            const std::string& simulationGases,
            const std::string& geometryFileUser,
            const double referenceArea,
            const double referenceLength,
            const int referenceAxis,
            const Eigen::Vector3d& momentReferencePoint,
            const double gridSpacing,
            const double simulatedParticlesPerCell,
            const double wallTemperature = 300.0,
            const double accommodationCoefficient = 1.0 );

    //! Default destructor.
    /*!
     *  Default destructor.
     */
    ~RarefiedFlowAnalysis( ) { }

    //! Get aerodynamic coefficients.
    /*!
     *  Returns aerodynamic coefficients.
     *  The physical meaning of each of the three independent variables is: 0 = altitude,
     *  1 = Mach number, 1 = angle of attack.
     * \param independentVariables Array of values of independent variable
     *          indices in dataPointsOfIndependentVariables_.
     * \return vector of coefficients at specified independent variable indices.
     */
    Eigen::Vector6d getAerodynamicCoefficientsDataPoint(
            const boost::array< int, 3 > independentVariables );

private:

    //! Open and read geometry file for SPARTA simulation.
    /*!
     * Open and read geometry file for SPARTA simulation, to extract information on the
     * vehicle and surface elements dimensions and properties.
     * \param geometryFileUser Path to the file describing the geometry of the vehicle, where
     *          the surface elements are discretized as triangles.
     */
    void analyzeGeometryFile( const std::string& geometryFileUser );

    //! Get conditions for simulation.
    /*!
     * Get conditions for simulation, including simulation environment boundaries, velocity of incoming flow,
     * time step of simulation (set as 10 % of the time needed for a particle to traverse the simulation
     * environment), and ratio of real-to-simulated particles.
     */
    void getSimulationConditions( );

    //! Generate aerodynamic database.
    /*!
     * Generates aerodynamic database, by running the SPARTA simulation via command line (standard
     * library function std::system), and reads output of simulation to extract values of aerodynamic
     * coefficients, as a function of altitude, Mach number and angle of attack.
     */
    void generateCoefficients( );

    //! Path to SPARTA executable.
    /*!
     * Path to SPARTA executable. Note that SPARTA is an external software and needs to be compiled before
     * it can be used in Tudat. See the instructions in the manual [2].
     */
    std::string SPARTAExecutable_;

    //! String of gases making up the atmosphere of the target planet.
    std::string simulationGases_;

    //! Reference axis for the aerodynamic analysis.
    int referenceAxis_;

    //! Grid size for simulation environment.
    /*!
     *  Grid size for simulation environment. Used to define the size of each cell, and the number of cells in
     *  the environment.
     */
    double gridSpacing_;

    //! Number of simulated particles per cell.
    double simulatedParticlesPerCell_;

    //! Temperature of surface of vehicle.
    double wallTemperature_;

    //! Accommodation coefficient of surface of vehicle.
    double accommodationCoefficient_;

    //! List of points making up the vehicle geometry.
    /*!
     * List of points making up the vehicle geometry. Column size is 3, since only triangular mesh are supported
     * by SPARTA.
     */
    Eigen::Matrix< double, Eigen::Dynamic, 3 > shapePoints_;

    //! List of triangle vertices making up the vehicle geometry.
    /*!
     * List of triangle vertices making up the vehicle geometry. Each element refers to a point from the shapePoints_
     * list. Column size is 3, since only triangular mesh are supported by SPARTA.
     */
    Eigen::Matrix< int, Eigen::Dynamic, 3 > shapeTriangles_;

    //! Total number of points making up the vehicle geometry.
    int numberOfPoints_;

    //! Total number of triangles making up the vehicle geometry.
    int numberOfTriangles_;

    //! Maximum dimensions of the vehicle in x, y and z.
    Eigen::Vector3d maximumDimensions_;

    //! Minimum dimensions of the vehicle in x, y and z.
    Eigen::Vector3d minimumDimensions_;

    //! List of normal vector to each surface element.
    Eigen::Matrix< double, 3, Eigen::Dynamic > elementSurfaceNormal_;

    //! List of surface area of each surface element.
    Eigen::Matrix< double, 1, Eigen::Dynamic > elementSurfaceArea_;

    //! List of moment arm for each surface element.
    /*!
     * List of moment arm for each surface element. Note that the moment arm is defined as the distance
     * of the centroid of each surface element, to a reference point, i.e., momentReferencePoint_.
     */
    Eigen::Matrix< double, 3, Eigen::Dynamic > elementMomentArm_;

    //! List of cross-sectional area of vehicle in x, y and z.
    Eigen::Vector3d shapeCrossSectionalArea_;

    //! Atmospheric conditions at altitudes specified by user in dataPointsOfIndependentVariables_.
    /*!
     * Atmospheric conditions at altitudes specified by user in dataPointsOfIndependentVariables_. The key values
     * are given by the enumeration AtmosphericConditionVariables, at the beginning of this file.
     */
    std::map< AtmosphericConditionVariables, std::vector< double > > atmosphericConditions_;

    //! List defining the simulation environment boundaries.
    /*!
     * List defining the simulation environment boundaries. They are computed from the maximum and minimum
     * dimensions of the vehicle, i.e., maximumDimensions_ and minimumDimensions_, by adding extra space in each
     * direction (50 % of the vehicle size).
     */
    Eigen::Vector6d simulationBoundaries_;

    //! List of grid spacing for each dimension.
    /*!
     * List of grid spacing for each dimension. They are determined by dividing the simulation dimensions (given
     * by simulationBoundaries_) by the grid spacing, i.e., gridSpacing_.
     */
    Eigen::Vector3d simulationGrid_;

    //! List of freestream velocities for each altitude and Mach number.
    Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > freeStreamVelocities_;

    //! List of time steps for the simulation for each altitude and Mach number.
    Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > simulationTimeStep_;

    //! List of ratios of real-to-simulated particles for each altitude.
    Eigen::Matrix< double, Eigen::Dynamic, 1 > ratioOfRealToSimulatedParticles_;

    //! String containing the template of the input file for the SPARTA simulation.
    /*!
     * String containing the template of the input file for the SPARTA simulation. This template contains
     * format specifiers (such as %.3f) where a simulation variable is to be placed. This is done by filling
     * the values input by the user and determined in previous parts of the code.
     */
    std::string inputTemplate_;
};

//! Typedef for shared-pointer to RarefiedFlowAnalysis object.
typedef boost::shared_ptr< RarefiedFlowAnalysis > RarefiedFlowAnalysisPointer;

} // namespace aerodynamics

} // namespace tudat

#endif // TUDAT_RAREFIED_FLOW_ANALYSIS_H
