/*!   \file hypersonicLocalInclinationAnalysis.h
 *    This file contains the definition of the hypersonic local inclination analysis class.
 *
 *    Path              : /Astrodynamics/ForceModels/
 *    Version           : 5
 *    Check status      : Checked
 *
 *    Author            : D. Dirkx
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : d.dirkx@tudelft.nl
 *
 *    Checker           : B. Romgens
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : B.Romgens@student.tudelft.nl
 *
 *    Date created      : 25 November, 2010
 *    Last modified     : 4 February, 2011
 *
 *    References
 *      Gentry, A., Smyth, D., and Oliver, W. The Mark IV Supersonic-Hypersonic
 *        Arbitrary Body Program, Volume II - Program Formulation, Douglas
 *        Aircraft Company, 1973.
 *      Anderson Jr. , J.D, Hypersonic and High-Temperature Gas Dynamics, 2nd
 *        edition, AIAA Education Series, 2006
 *
 *    Notes
 *
 *    Copyright (c) 2010-2011 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      102511    D. Dirkx          First version of file.
 *      110501    D. Dirkx          Added more comments.
 *      112701    D. Dirkx          Finalized for code check.
 *      110131    B. Romgens        Minor modifications during code check.
 *      110204    D. Dirkx          Finalized code.
 */

#ifndef HYPERSONICLOCALINCLINATIONANALYSIS_H
#define HYPERSONICLOCALINCLINATIONANALYSIS_H

// Include statements.
#include <iostream>
#include "Astrodynamics/Bodies/vehicle.h"
#include "Astrodynamics/ForceModels/aerodynamicCoefficientGenerator.h"
#include "Mathematics/GeometricShapes/lawgsPartGeometry.h"

//! Class for inviscid hypersonic aerodynamic analysis using local inclination
//! methods.
/*! Class for inviscid hypersonic aerodynamic analysis using local inclination
 *  methods. These methods assume that the local pressure on the vehicle is only
 *  dependent on the local inclination angle w.r.t. the freestream flow and
 *  freestream conditions, such as Mach number and ratio of specific heats. The
 *  setVehicle function needs to be used to set the Vehicle which is to be analyzed.
 *  This Vehicle must have a VehicleExternalModel, containing a SurfaceGeometry.
 *  Values for the Mach number, angle of attack and angle of sideslip can be set
 *  to default values (see setDefaultMachPoints( ), etc.) or manually by first
 *  using the setNumberOfMachPoints function, and then the setMachPoint function
 *  for each of the values of the independent variables. All aerodynamic
 *  coefficients can be calculated using the generateDatabase function, or on an
 *  as needed basis by using the getAerodynamicCoefficients function.
 */
class HypersonicLocalInclinationAnalysis: public AerodynamicCoefficientGenerator
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    HypersonicLocalInclinationAnalysis( ) : vehicleParts_( NULL ), numberOfVehicleParts_( -0 ),
        inclination_( NULL ), pressureCoefficient_( NULL ), stagnationPressureCoefficient( -0.0 ),
        ratioOfSpecificHeats( 1.4 ),  vehicleName_( " " ), machRegime_( "Full" ),
        selectedMethods_( NULL )
    {
        machIndex_ = 0; angleOfAttackIndex_ = 1; angleOfSideslipIndex_ = 2;
        setNumberOfIndependentVariables( 3 );
    }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    virtual ~HypersonicLocalInclinationAnalysis( );

    //! Constructor from a Vehicle.
    /*!
     * This constructor sets the geometry which is used for the analysis from
     * a Vehicle object. Vehicle must have an external model containing a
     * surface geometry. If surface geometry is a CompositeSurfaceGeometry,
     * it must not contain any CompositeSurfaceGeometries itself.
     * \param vehicle Vehicle which is to be analyzed.
     * \param numberOfLines Array of size equal to number of SingleSurfaceGeometries to set the
     *          number of lines in the resulting LaWGS parts.
     * \param numberOfPoints Array of size equal to number of SingleSurfaceGeometries to set the
     *          number of points per line in the resulting LaWGS parts.
     * \param invertOrders Array of size equal to number of SingleSurfaceGeometries to set whether
     *          to invert the panel orientation of the resulting LaWGS parts.
     */
    void setVehicle( Vehicle& vehicle, int* numberOfLines, int* numberOfPoints,
                     bool* invertOrders );

    //! Get aerodynamic coefficients.
    /*!
     * Returns aerodynamic coefficients.
     * \param independentVariables Array of values of independent variable
     *          indices in dataPointsOfIndependentVariables_.
     */
    VectorXd getAerodynamicCoefficients( int* independentVariables );

    //! Set local inclination methods for all parts (expansion and compression).
    /*!
     * Sets the aerodynamic analysis methods.
     * \param selectedMethods two-dimensional array (2 by numberOfVehicleParts_)
     * containing the identifiers.
     * for the methods.
     * For the first index:
     * 0 = compression
     * 1 = expansion
     * For the second index, if compression:
     * 0 = Newtonian
     * 1 = Modified Newtonian
     * 2 = Modified Newtonian-Prandtl Meyer (currently disabled)
     * 3 = Tangent Wedge (currently disabled)
     * 4 = Tangent Wedge Empirical (only for ratioOfSpecificHeats = 1.4)
     * 5 = Tangent Cone Empirical (only for ratioOfSpecificHeats = 1.4)
     * 6 = Modified Dahlem - Buck
     * 7 = Van Dyke Unified
     * 8 = Smyth Delta WIng
     * 9 = Hankey Flat Surface
     * For second index, if expansion:
     * 0 = Vacuum
     * 1 = Newtonian
     * 2 = Modified Newtonian - Prandtl Meyer (currently disabled)
     * 3 = Prandtl Meyer expansion from freestream (only for ratioOfSpecificHeats = 1.4)
     * 4 = High Mach Base Pressure
     * 5 = Van Dyke Unified
     * 6 = ACM empirical
     */
    void setSelectedMethods( int** selectedMethods );

    //! Set an analysis method on a single vehicle part.
    /*!
     * Sets an analysis method on a single vehicle part.
     * \param method Identifier of method which is set.
     * \param type Type of method which is set:
     * 0 = High hypersonic compression.
     * 1 = High hypersonic expansion.
     * 2 = Low hypersonic compression.
     * 3 = Low hypersonic expansion.
     * \param part Vehicle part on which to apply method.
     */
    void setSelectedMethod( const int& method, const int& type, const int& part )
    { selectedMethods_[ type ][ part ] = method; }

    //! Generate aerodynamic database.
    /*!
     * Generates aerodynamic database. Settings of geometry,
     * reference quantities, database point settings and analysis methods
     *  should have been set previously.
     */
    void generateDatabase( );

    //! Determine inclination angles of panels on a given part.
    /*!
     * Determines panel inclinations for all panels on a given part for given attitude.
     * \param partNumber Index from vehicleParts_ array for which to determine coefficients.
     * \param angleOfAttack Angle of attack at which to determine inclination angles.
     * \param angleOfSideslip Angle of sideslip at which to determine inclination angles.
     */
    virtual void determineInclination( const int& partNumber, const double& angleOfAttack,
                                       const double& angleOfSideslip );

    //! Get the number of vehicle parts.
    /*!
     * Returns the number of vehicle parts.
     */
    int getNumberOfVehicleParts( ) { return numberOfVehicleParts_; }

    //! Allocate aerodynamic coefficient array and NULL independent variables.
    /*!
     * Checks if all independent variables have been set, if not default values are set, depending
     * on which machRegime_ is selected. Subsequently, the vehicleCoefficients_ array is allocated.
     */
    void allocateVehicleCoefficients( );

    //! Allocate pressure coefficient, inclination, and method arrays.
    /*!
     * Allocate pressure coefficient, inclination, and method arrays.
     */
    void allocateArrays( );

    //! Get a vehicle part.
    /*!
     * Returns a vehicle part.
     * \param vehicleIndex Index in vehicleParts_ to be retrieved.
     * \return Requested vehicle part.
     */
     LawgsPartGeometry getVehiclePart( const int& vehicleIndex )
     { return vehicleParts_[ vehicleIndex ]; }

    //! Set mach regime.
    /*!
     * Sets mach regime, see machRegime_.
     * \param machRegime Name of Mach regime which is set.
     */
    void setMachRegime( const std::string& machRegime ) { machRegime_ = machRegime; }

    //! Get mach regime.
    /*!
     * Returns mach regime, see machRegime_.
     * \return Mach regime.
     */
    std::string getMachRegime( ) { return machRegime_; }

    //! Set the vehicle name.
    /*!
     * Sets the vehicle name.
     * \param vehicleName vehicle name.
     */
    void setVehicleName( const std::string& vehicleName ) { vehicleName_ = vehicleName; }

    //! Get the vehicle name.
    /*!
     * Returns the vehicle name.
     * \return Vehicle name.
     */
    std::string getVehicleName( ) { return vehicleName_; }

    //! Overload ostream to print class information.
    /*!
     * Overloads ostream to print class information, prints the number of lawgs geometry parts and
     * names.
     * \param stream Stream object.
     * \param hypersonicLocalInclinationAnalysis Hypersonic local inclination analysis.
     * \return Stream object.
     */
    friend std::ostream& operator<<( std::ostream& stream,
                                     HypersonicLocalInclinationAnalysis&
                                     hypersonicLocalInclinationAnalysis );

private:

    //! Generate aerodynamic coefficients at a single set of independent variables.
    /*!
     * Generates aerodynamic coefficients at a single set of independent variables.
     * Determines values and sets corresponding entry in vehicleCoefficients_ array.
     * \param independentVariableIndices Array of indices from lists of Mach number,
     *          angle of attack and angle of sideslip points at which to perform analysis.
     */
    void determineVehicleCoefficients( int* independentVariableIndices );

    //! Determine aerodynamic coefficients for a single LaWGS part.
    /*!
     * Determines aerodynamic coefficients for a single LaWGS part,
     * calls determinepressureCoefficient_ function for given vehicle part.
     * \param partNumber Index from vehicleParts_ array for which to determine coefficients.
     * \param independentVariableIndices Array of indices of independent variables.
     */
    VectorXd determinePartCoefficients( const int& partNumber, int* independentVariableIndices );

    //! Determine pressure coefficients on a given part.
    /*!
     * Determines pressure coefficients on a single vehicle part.
     * Calls the updateExpansionPressures and updateCompressionPressures for given vehicle part.
     * \param partNumber Index from vehicleParts_ array for which to determine coefficients.
     * \param independentVariableIndices Array of indices of independent variables.
     */
    void determinePressureCoefficients( const int& partNumber, int* independentVariableIndices );

    //! Determine force coefficients of a part.
    /*!
     * Sums the pressure coefficients of given part and determines force coefficients from it by
     * non-dimensionalization with reference area.
     * \param partNumber Index from vehicleParts_ array for which determine coefficients.
     */
    VectorXd calculateForceCoefficients( const int& partNumber );

    //! Determine moment coefficients of a part.
    /*!
     * Determines the moment coefficients of a given part by summing the contributions of all
     * panels on the part. Moment arms are taken from panel centroid to momentReferencePoint. Non-
     * dimensionalization is performed by product of referenceLength and referenceArea.
     * \param partNumber Index from vehicleParts_ array for which to determine coefficients.
     */
    VectorXd calculateMomentCoefficients( const int& partNumber );

    //! Determine the compression pressure coefficients of a given part.
    /*!
     * Sets the values of pressureCoefficient_ on given part and at given Mach number for which
     * inclination > 0.
     * \param machNumber Mach number at which to perform analysis.
     * \param partNumber of part from vehicleParts_ which is to be analyzed.
     */
    void updateCompressionPressures( const double& machNumber, const int& partNumber );

    //! Determine the expansion pressure coefficients of a given part.
    /*!
     * Determine the values of pressureCoefficient_ on given part and at given Mach number for
     * which inclination <= 0.
     * \param machNumber Mach number at which to perform analysis.
     * \param partNumber of part from vehicleParts_ which is to be analyzed.
     */
    void updateExpansionPressures( const double& machNumber, const int& partNumber );

    //! Set the default analysis points for Mach number.
    /*!
     * Sets the default analysis points for Mach number, depending on which Mach regime has been
     * set.
     */
    void setDefaultMachPoints( );

    //! Set the default analysis points for angle of attack
    /*!
     * Set the default analysis points for angle of attack.
     */
    void setDefaultAngleOfAttackPoints( );

    //! Set the default analysis points for angle of sideslip.
    /*!
     * Sets the default analysis points for angle of sideslip.
     */
    void setDefaultAngleOfSideslipPoints( );

    //! Array of vehicle parts.
    /*!
     * Array of vehicle parts.
     */
    LawgsPartGeometry* vehicleParts_;

    //! Number of entries in vehicleParts_ member variable.
    /*!
     * Number of entries in vehicleParts_ member variable.
     */
    int numberOfVehicleParts_;

    //! Three-dimensional array of panel inclination angles.
    /*!
     * Three-dimensional array of panel inclination angles at current values of
     * independent variables. Indices indicate part-line-point.
     */
    double*** inclination_;

    //! Three-dimensional array of panel pressure coefficients.
    /*!
     * Three-dimensional array of panel pressure coefficients at current values
     * of independent variables. Indices indicate part-line-point.
     */
    double*** pressureCoefficient_;

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
    //! Vehicle name.
    /*!
     * Vehicle name.
     */
    std::string vehicleName_;

    //! Mach regime.
    /*!
     * Mach regime, permissible values are "Full", "High" or "Low", default is "Full".
     */
    std::string machRegime_;

    //! Array of selected methods.
    /*!
     * Array of selected methods, first index represents compression/expansion,
     * second index represents vehicle part.
     */
    int** selectedMethods_;
};

#endif // HYPERSONICLOCALINCLINATIONANALYSIS_H

// End of file.
