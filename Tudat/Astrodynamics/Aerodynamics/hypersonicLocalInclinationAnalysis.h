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

#ifndef TUDAT_HYPERSONIC_LOCAL_INCLINATION_ANALYSIS_H
#define TUDAT_HYPERSONIC_LOCAL_INCLINATION_ANALYSIS_H

#include <map>
#include <string>
#include <vector>

#include <boost/array.hpp>
#include <boost/multi_array.hpp>
#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/Aerodynamics/aerodynamicCoefficientGenerator.h"
#include "Tudat/Basics/basicTypedefs.h"
#include "Tudat/Mathematics/GeometricShapes/lawgsPartGeometry.h"

namespace tudat
{
namespace aerodynamics
{

//! Returns default values of mach number for use in HypersonicLocalInclinationAnalysis.
/*!
 *  Returns default values of mach number for use in HypersonicLocalInclinationAnalysis.
 */
std::vector< double > getDefaultHypersonicLocalInclinationMachPoints(
        const std::string& machRegime );

//! Returns default values of angle of attack for use in HypersonicLocalInclinationAnalysis.
/*!
 *  Returns default values of angle of attack for use in HypersonicLocalInclinationAnalysis.
 */
std::vector< double > getDefaultHypersonicLocalInclinationAngleOfAttackPoints( );

//! Returns default values of angle of sideslip for use in HypersonicLocalInclinationAnalysis.
/*!
 *  Returns default values of angle of sideslip for use in HypersonicLocalInclinationAnalysis.
 */
std::vector< double > getDefaultHypersonicLocalInclinationAngleOfSideslipPoints( );

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
class HypersonicLocalInclinationAnalysis: public AerodynamicCoefficientGenerator< 3, 6 >
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
    HypersonicLocalInclinationAnalysis(
            const std::vector< std::vector< double > >& dataPointsOfIndependentVariables,
            const boost::shared_ptr< SurfaceGeometry > inputVehicleSurface,
            const std::vector< int >& numberOfLines,
            const std::vector< int >& numberOfPoints,
            const std::vector< bool >& invertOrders,
            const std::vector< std::vector< int > >& selectedMethods,
            const double referenceArea,
            const double referenceLength,
            const Eigen::Vector3d& momentReferencePoint );

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~HypersonicLocalInclinationAnalysis( ) { }

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

    //! Determine inclination angles of panels on a given part.
    /*!
     * Determines panel inclinations for all panels on all parts for given attitude.
     * Outward pointing surface-normals are assumed!
     * \param angleOfAttack Angle of attack at which to determine inclination angles.
     * \param angleOfSideslip Angle of sideslip at which to determine inclination angles.
     */
    void determineInclinations( const double angleOfAttack,
                                const double angleOfSideslip );

    //! Get the number of vehicle parts.
    /*!
     *  Returns the number of vehicle parts.
     *  \return The number of vehicle parts.
     */
    int getNumberOfVehicleParts( ) const
    {
        return vehicleParts_.size( );
    }

    //! Get a vehicle part.
    /*!
     * Returns a vehicle part.
     * \param vehicleIndex Index in vehicleParts_ to be retrieved.
     * \return Requested vehicle part.
     */
     boost::shared_ptr< geometric_shapes::LawgsPartGeometry > getVehiclePart(
             const int vehicleIndex ) const
     {
         return vehicleParts_[ vehicleIndex ];
     }

    //! Overload ostream to print class information.
    /*!
     * Overloads ostream to print class information, prints the number of lawgs geometry parts and
     * names.
     * \param stream Stream object.
     * \param hypersonicLocalInclinationAnalysis Hypersonic local inclination analysis.
     * \return Stream object.
     */
    friend std::ostream& operator << ( std::ostream& stream,
                                     HypersonicLocalInclinationAnalysis&
                                     hypersonicLocalInclinationAnalysis );

private:

    //! Generate aerodynamic database.
    /*!
     * Generates aerodynamic database. Settings of geometry,
     * reference quantities, database point settings and analysis methods
     *  should have been set previously.
     */
    void generateCoefficients( );

    //! Generate aerodynamic coefficients at a single set of independent variables.
    /*!
     * Generates aerodynamic coefficients at a single set of independent variables.
     * Determines values and sets corresponding entry in vehicleCoefficients_ array.
     * \param independentVariableIndices Array of indices from lists of Mach number,
     *          angle of attack and angle of sideslip points at which to perform analysis.
     */
    void determineVehicleCoefficients( const boost::array< int, 3 > independentVariableIndices );

    //! Determine aerodynamic coefficients for a single LaWGS part.
    /*!
     * Determines aerodynamic coefficients for a single LaWGS part,
     * calls determinepressureCoefficient_ function for given vehicle part.
     * \param partNumber Index from vehicleParts_ array for which to determine coefficients.
     * \param independentVariableIndices Array of indices of independent variables.
     * \return Force and moment coefficients for requested vehicle part.
     */
    Eigen::Vector6d determinePartCoefficients(
            const int partNumber, const boost::array< int, 3 > independentVariableIndices );

    //! Determine pressure coefficients on a given part.
    /*!
     * Determines pressure coefficients on a single vehicle part.
     * Calls the updateExpansionPressures and updateCompressionPressures for given vehicle part.
     * \param partNumber Index from vehicleParts_ array for which to determine coefficients.
     * \param independentVariableIndices Array of indices of independent variables.
     */
    void determinePressureCoefficients( const int partNumber,
                                        const boost::array< int, 3 > independentVariableIndices );

    //! Determine force coefficients of a part.
    /*!
     * Sums the pressure coefficients of given part and determines force coefficients from it by
     * non-dimensionalization with reference area.
     * \param partNumber Index from vehicleParts_ array for which determine coefficients.
     * \return Force coefficients for requested vehicle part.
     */
    Eigen::Vector3d calculateForceCoefficients( const int partNumber );

    //! Determine moment coefficients of a part.
    /*!
     * Determines the moment coefficients of a given part by summing the contributions of all
     * panels on the part. Moment arms are taken from panel centroid to momentReferencePoint. Non-
     * dimensionalization is performed by product of referenceLength and referenceArea.
     * \param partNumber Index from vehicleParts_ array for which to determine coefficients.
     * \return Moment coefficients for requested vehicle part.
     */
    Eigen::Vector3d calculateMomentCoefficients( const int partNumber );

    //! Determine the compression pressure coefficients of a given part.
    /*!
     * Sets the values of pressureCoefficient_ on given part and at given Mach number for which
     * inclination > 0.
     * \param machNumber Mach number at which to perform analysis.
     * \param partNumber of part from vehicleParts_ which is to be analyzed.
     */
    void updateCompressionPressures( const double machNumber, const int partNumber );

    //! Determine the expansion pressure coefficients of a given part.
    /*!
     * Determine the values of pressureCoefficient_ on given part and at given Mach number for
     * which inclination <= 0.
     * \param machNumber Mach number at which to perform analysis.
     * \param partNumber of part from vehicleParts_ which is to be analyzed.
     */
    void updateExpansionPressures( const double machNumber, const int partNumber );

    //! Array of vehicle parts.
    /*!
     * Array of vehicle parts.
     */
    std::vector< boost::shared_ptr< geometric_shapes::LawgsPartGeometry > > vehicleParts_;

    //! Multi-array as which indicates which coefficients have been calculated already.
    /*!
     * Multi-array as which indicates which coefficients have been calculated already. Indices of
     * entries coincide with indices of aerodynamicCoefficients_.
     */
    boost::multi_array< bool, 3 > isCoefficientGenerated_;

    //! Three-dimensional array of panel inclination angles.
    /*!
     * Three-dimensional array of panel inclination angles at current values of
     * independent variables. Indices indicate part-line-point.
     */
    std::vector< std::vector< std::vector< double > > > inclination_;

    //! Map of angle of attack and -sideslip pair and associated panel inclinations.
    /*!
     * Map of angle of attack and -sideslip pair and associated panel inclinations.
     */
    std::map< std::pair< double, double >, std::vector< std::vector< std::vector< double > > > >
            previouslyComputedInclinations_;

    //! Three-dimensional array of panel pressure coefficients.
    /*!
     * Three-dimensional array of panel pressure coefficients at current values
     * of independent variables. Indices indicate part-line-point.
     */
    std::vector< std::vector< std::vector< double > > > pressureCoefficient_;

    //! Stagnation pressure coefficient.
    /*!
     * Stagnation pressure coefficient for flow which has passed through a
     * normal shock wave at current Mach number.
     */
    double stagnationPressureCoefficient;

    //! Ratio of specific heats.
    /*!
     * Ratio of specific heat at constant pressure to specific heat at constant pressure.
     */
    double ratioOfSpecificHeats;
\
    //! Array of selected methods.
    /*!
     * Array of selected methods, first index represents compression/expansion,
     * second index represents vehicle part.
     */
    std::vector< std::vector< int > > selectedMethods_;
};

//! Typedef for shared-pointer to HypersonicLocalInclinationAnalysis object.
typedef boost::shared_ptr< HypersonicLocalInclinationAnalysis >
HypersonicLocalInclinationAnalysisPointer;

} // namespace aerodynamics
} // namespace tudat

#endif // TUDAT_HYPERSONIC_LOCAL_INCLINATION_ANALYSIS_H
