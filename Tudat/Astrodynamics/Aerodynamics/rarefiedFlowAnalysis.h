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
 *      Gentry, A., Smyth, D., and Oliver, W. The Mark IV Supersonic-Hypersonic Arbitrary Body
 *        Program, Volume II - Program Formulation, Douglas Aircraft Company, 1973.
 *      Anderson Jr. , J.D, Hypersonic and High-Temperature Gas Dynamics, 2nd edition,
 *        AIAA Education Series, 2006
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
#include "Tudat/Astrodynamics/Aerodynamics/standardAtmosphere.h"
#include "Tudat/InputOutput/basicInputOutput.h"
#include "Tudat/InputOutput/matrixTextFileReader.h"
#include "Tudat/InputOutput/SPARTADataReader.h"
#include "Tudat/Basics/basicTypedefs.h"

namespace tudat
{
namespace aerodynamics
{

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

//! Returns default values of angle of sideslip for use in RarefiedFlowAnalysis.
/*!
 *  Returns default values of angle of sideslip for use in RarefiedFlowAnalysis.
 */
std::vector< double > getDefaultRarefiedFlowAngleOfSideslipPoints( );

//! Class for inviscid hypersonic aerodynamic analysis using local inclination methods.
/*!
 * Class for inviscid hypersonic aerodynamic analysis using local inclination
 * methods. These methods assume that the local pressure on the vehicle is only
 * dependent on the local inclination angle w.r.t. the freestream flow and
 * freestream conditions, such as Mach number and ratio of specific heats.
 * All aerodynamic coefficients can be calculated using the generateCoefficients function, or on an
 * as needed basis by using the getAerodynamicCoefficientsDataPoint function. Note that during the
 * panel inclination determination process, a geometry with outward surface-normals is assumed.
 * The resulting coefficients are expressed in the same reference frame as that of the input
 * geometry.
 */
class RarefiedFlowAnalysis: public AerodynamicCoefficientGenerator< 3, 6 >
{
public:

    //! Default constructor.
    /*!
     *  Default constructor of class, specified vehicle geometry, discretization properties,
     *  independent variable ranges, reference values and local inclination methods that are
     *  to be used
     *  \param dataPointsOfIndependentVariables Vector of vector, with each subvector containing
     *  the data points of each of the independent variables for the coefficient generation.
     *  The physical meaning of each of the three independent variables is: 0 = mach numner,
     *  1 = angle of attack, 2 = angle of sideslip.
     *  Each of the subvectors must be sorted in ascending order.
     *  \param inputVehicleSurface Vehicle surface geometry for which the coefficients are to be
     *  determined.
     *  \param numberOfLines Number of discretization points in the first independent surface
     *  variable of each of the subparts of inputVehicleSurface.
     *  \param numberOfPoints Number of discretization points in the second independent surface
     *  variable of each of the subparts of inputVehicleSurface.
     *  \param invertOrders Booleans to denote whether the surface normals of the panels of
     *  each discretized inputVehicleSurface subpart are to be inverted
     *  (i.e. inward-facing->outward facing or vice versa)
     *  \param selectedMethods Array of selected local inclination methods, the first index
     *  represents compression/expansion, the second index denotes the vehicle part index.
     *  The value for each separate method can be found in the updateCompressionPressures and
     *  updateExpansionPressures functions, respectively.
     *  \param referenceArea Reference area used to non-dimensionalize aerodynamic forces
     *  and moments.
     *  \param referenceLength Reference length used to non-dimensionalize aerodynamic moments.
     *  \param momentReferencePoint Reference point wrt which aerodynamic moments are calculated.
     */
    RarefiedFlowAnalysis(
            const std::string& SPARTAExecutable,
            const std::vector< std::vector< double > >& dataPointsOfIndependentVariables,
            const std::string& simulationGases,
            const StandardAtmosphere& atmosphereModel,
            const std::string& geometryFileUser,
            const int& referenceAxis,
            const Eigen::Vector3d& momentReferencePoint,
            const double wallTemperature = 300.0,
            const double accomodationCoefficient = 1.0 );

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~RarefiedFlowAnalysis( ) { }

    //! Get aerodynamic coefficients.
    /*!
     *  Returns aerodynamic coefficients.
     *  The physical meaning of each of the three independent variables is: 0 = mach numner,
     *  1 = angle of attack, 2 = angle of sideslip.
     * \param independentVariables Array of values of independent variable
     *          indices in dataPointsOfIndependentVariables_.
     * \return vector of coefficients at specified independent variable indices.
     */
    Eigen::Vector6d getAerodynamicCoefficientsDataPoint(
            const boost::array< int, 3 > independentVariables );

private:

    void analyzeGeometryFile( const std::string& geometryFileUser );

    void getSimulationConditions( );

    void runSPARTASimulation( const std::string& SPARTAExecutable );

    std::string outputDirectory_ = "results";
    std::string outputPath_ = input_output::getSPARTADataPath( ) + outputDirectory_;
    std::string inputFile_ = input_output::getSPARTADataPath( )  + "in.sparta";
    std::string inputFileTemplate_ = input_output::getSPARTADataPath( )  + "SPARTAInputTemplate.txt";
    std::string geometryFileInternal_ = input_output::getSPARTADataPath( ) + "data.shape";
    double gridSpacing_ = 0.25;
    double simulatedParticlesPerCell_ = 15;
    std::string simulationGases_;
    int referenceAxis_;
    double wallTemperature_;
    double accomodationCoefficient_;
    Eigen::Matrix< double, Eigen::Dynamic, 3 > shapePoints_;
    Eigen::Matrix< int, Eigen::Dynamic, 3 > shapeTriangles_;
    int numberOfPoints_;
    int numberOfTriangles_;
    Eigen::Vector3d maximumDimensions_;
    Eigen::Vector3d minimumDimensions_;
    Eigen::Matrix< double, 3, Eigen::Dynamic > elementSurfaceNormal_;
    Eigen::Matrix< double, 1, Eigen::Dynamic > elementSurfaceArea_;
    Eigen::Matrix< double, 3, Eigen::Dynamic > elementMomentArm_;
    Eigen::Vector3d shapeCrossSectionalArea_;
    std::vector< std::vector< double > > atmosphericConditions_;
    Eigen::Vector6d simulationBoundaries_;
    Eigen::Vector3d simulationGrid_;
    Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > freeStreamVelocities_;
    Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > simulationTimeStep_;
    Eigen::Matrix< double, Eigen::Dynamic, 1 > ratioOfRealToSimulatedParticles_;
    std::string inputTemplate_;
    boost::multi_array< Eigen::Vector6d, 3 > aerodynamicCoefficients_;
};

//! Typedef for shared-pointer to RarefiedFlowAnalysis object.
typedef boost::shared_ptr< RarefiedFlowAnalysis >
RarefiedFlowAnalysisPointer;

} // namespace aerodynamics
} // namespace tudat

#endif // TUDAT_RAREFIED_FLOW_ANALYSIS_H
