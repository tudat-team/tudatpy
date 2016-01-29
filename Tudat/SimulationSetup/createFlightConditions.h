#ifndef CREATEFLIGHTCONDITIONS_H
#define CREATEFLIGHTCONDITIONS_H

#include <boost/multi_array.hpp>

#include <vector>

#include "Tudat/Astrodynamics/Aerodynamics/flightConditions.h"
#include "Tudat/SimulationSetup/body.h"

namespace tudat
{

namespace simulation_setup
{


//! List of aerodynamic coefficient models available in simulations
/*!
 *  List of aerodynamic coefficient models available in simulations. Aerodynamic coefficient models
 *  not defined by this given enum cannot be used for automatic model setup.
 */
enum AerodynamicCoefficientTypes
{
    constant_aerodynamic_coefficients,
    hypersonic_local_inclincation_coefficients,
    tabulated_coefficients
};

//! Class for providing settings for aerodynamic coefficient model.
/*!
 *  Class for providing settings for automatic aerodynamic coefficient model creation. This class is
 *  a functional (base) class for settings of aerodynamic coefficient models that require no
 *  information in addition to their type. Aerodynamic coefficient model classes defining requiring
 *  additional information must be created using an object derived from this class.
 */
class AerodynamicCoefficientSettings
{
public:

    //! Constructor, sets type of aerodynamic coefficient model.
    /*!
     *  Constructor, sets type of aerodynamic coefficient model. Settings for aerodynamic coefficient
     *  models requiring additional information should be defined in a derived class.
     *  \param aerodynamicCoefficientTypes Type of aerodynamic coefficient model that is to be created.
     *  \param referenceLength Reference length with which aerodynamic moments (about x- and z- axes) are non-dimensionalized.
     *  \param referenceArea Reference area with which aerodynamic forces and moments are non-dimensionalized.
     *  \param lateralReferenceLength Reference length with which aerodynamic moments (about y-axis) is non-dimensionalized.
     *  \param momentReferencePoint Point w.r.t. aerodynamic moment is calculated
     *  \param independentVariableNames Vector with identifiers the physical meaning of each
     *  independent variable of the aerodynamic coefficients.
     */
    AerodynamicCoefficientSettings(
            const AerodynamicCoefficientTypes aerodynamicCoefficientTypes,
            const double referenceLength,
            const double referenceArea,
            const double lateralReferenceLength,
            const Eigen::Vector3d& momentReferencePoint,
            const std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables >
            independentVariableNames,
            const bool areCoefficientsInAerodynamicFrame = 1,
            const bool areCoefficientsInNegativeAxisDirection = 1 ):
        aerodynamicCoefficientTypes_( aerodynamicCoefficientTypes ),
        referenceLength_( referenceLength ), referenceArea_( referenceArea ),
        lateralReferenceLength_( lateralReferenceLength ),
        momentReferencePoint_( momentReferencePoint ),
        independentVariableNames_( independentVariableNames ),
        areCoefficientsInAerodynamicFrame_( areCoefficientsInAerodynamicFrame ),
        areCoefficientsInNegativeAxisDirection_( areCoefficientsInNegativeAxisDirection ){ }

    //! Destructor
    virtual ~AerodynamicCoefficientSettings( ){ }

    //! Function to return type of aerodynamic coefficient model that is to be created.
    /*!
     *  Function to return type of aerodynamic coefficient model that is to be created.
     *  \return Type of aerodynamic coefficient model that is to be created.
     */
    AerodynamicCoefficientTypes getAerodynamicCoefficientType( ){ return aerodynamicCoefficientTypes_; }

    //! Get reference area.
    /*!
     * Returns reference area used to non-dimensionalize aerodynamic forces and moments.
     * \return Aerodynamic reference area.
     */
    double getReferenceArea( ) { return referenceArea_; }

    //! Get reference length.
    /*!
     * Returns reference length used to non-dimensionalize aerodynamic moments.
     * \return Aerodynamic reference length.
     */
    double getReferenceLength( ) { return referenceLength_; }

    //! Get lateral reference length.
    /*!
     * Returns lateral reference length used to non-dimensionalize aerodynamic moments.
     * \return Aerodynamic lateral reference length.
     */
    double getLateralReferenceLength( ) { return lateralReferenceLength_; }

    //! Get moment reference point.
    /*!
     * Returns the point w.r.t. which the arm of the aerodynamic moment on a vehicle panel is
     * determined.
     * \return Aerodynamic reference point.
     */
    Eigen::VectorXd getMomentReferencePoint( ) { return momentReferencePoint_; }


    std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables >
    getIndependentVariableNames( )
    {
        return independentVariableNames_;
    }

    bool getAreCoefficientsInAerodynamicFrame( )
    {
        return areCoefficientsInAerodynamicFrame_;
    }

    bool getAreCoefficientsInNegativeAxisDirection( )
    {
        return areCoefficientsInNegativeAxisDirection_;
    }

private:

    //!  Type of atmosphere model that is to be created.
    AerodynamicCoefficientTypes aerodynamicCoefficientTypes_;

    //! Aerodynamic reference length.
    /*!
     * Reference length with which aerodynamic moments are non-dimensionalized.
     */
    double referenceLength_;

    //! Aerodynamic reference area.
    /*!
     * Reference area with which aerodynamic forces and moments are non-dimensionalized.
     */
    double referenceArea_;

    //! Lateral aerodynamic reference length.
    /*!
     * Lateral reference length with which aerodynamic moments are non-dimensionalized.
     */
    double lateralReferenceLength_;

    //! Aerodynamic moment reference point.
    /*!
     * Point w.r.t. which the arm of the moment on a vehicle panel is determined.
     */
    Eigen::Vector3d momentReferencePoint_;

    std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables >
    independentVariableNames_;

    bool areCoefficientsInAerodynamicFrame_;

    bool areCoefficientsInNegativeAxisDirection_;
};

//! AerodynamicCoefficientSettings for defining a constant aerodynamic coefficients
class ConstantAerodynamicCoefficientSettings: public AerodynamicCoefficientSettings
{
public:
    //! Constructor.
    /*!
     *  Constructor.
     *  \param constantForceCoefficient Constant force coefficients.
     *  \param constantMomentCoefficient Constant moment coefficients.
     *  \param referenceLength Reference length with which aerodynamic moments (about x- and z- axes) are non-dimensionalized.
     *  \param referenceArea Reference area with which aerodynamic forces and moments are non-dimensionalized.
     *  \param lateralReferenceLength Reference length with which aerodynamic moments (about y-axis) is non-dimensionalized.
     *  \param momentReferencePoint Point w.r.t. aerodynamic moment is calculated
     *  \param independentVariableNames Vector with identifiers the physical meaning of each
     *  independent variable of the aerodynamic coefficients.
     *  \param areCoefficientsInAerodynamicFrame Boolean to define whether the aerodynamic
     *  coefficients are defined in the aerodynamic frame (lift, drag, side force) or in the body
     *  frame (typically denoted as Cx, Cy, Cz).
     *  \param areCoefficientsInNegativeAxisDirection Boolean to define whether the aerodynamic
     *  coefficients are positiver along tyhe positive axes of the body or aerodynamic frame
     *  (see areCoefficientsInAerodynamicFrame). Note that for (lift, drag, side force), the
     *  coefficients are typically defined in negative direction.
     */
    ConstantAerodynamicCoefficientSettings(
            const double referenceLength,
            const double referenceArea,
            const double lateralReferenceLength,
            const Eigen::Vector3d& momentReferencePoint,
            const Eigen::Vector3d constantForceCoefficient,
            const Eigen::Vector3d constantMomentCoefficient = Eigen::Vector3d::Zero( ),
            const bool areCoefficientsInAerodynamicFrame = 0,
            const bool areCoefficientsInNegativeAxisDirection = 1  ):
        AerodynamicCoefficientSettings(
            constant_aerodynamic_coefficients, referenceLength, referenceArea,
            lateralReferenceLength, momentReferencePoint,
            std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables >( ),
            areCoefficientsInAerodynamicFrame, areCoefficientsInNegativeAxisDirection ),
        constantForceCoefficient_( constantForceCoefficient ),
        constantMomentCoefficient_( constantMomentCoefficient )
    { }


    Eigen::Vector3d getConstantForceCoefficient( )
    {
        return  constantForceCoefficient_;
    }

    Eigen::Vector3d getConstantMomentCoefficient( )
    {
        return constantMomentCoefficient_;
    }

private:

    Eigen::Vector3d constantForceCoefficient_;
    Eigen::Vector3d constantMomentCoefficient_;
};

template< int NumberOfDimensions >
class TabulatedAerodynamicCoefficientSettings: public AerodynamicCoefficientSettings
{
public:
    TabulatedAerodynamicCoefficientSettings(
            const std::vector< std::vector< double > > independentVariables,
            const boost::multi_array< Eigen::Vector3d, NumberOfDimensions > forceCoefficients,
            const boost::multi_array< Eigen::Vector3d, NumberOfDimensions > momentCoefficients,
            const double referenceLength,
            const double referenceArea,
            const double lateralReferenceLength,
            const Eigen::Vector3d& momentReferencePoint,
            const std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables >
            independentVariableNames,
            const bool areCoefficientsInAerodynamicFrame = 1,
            const bool areCoefficientsInNegativeAxisDirection = 1 ):
        AerodynamicCoefficientSettings(
            tabulated_coefficients, referenceLength, referenceArea,
            lateralReferenceLength, momentReferencePoint,
            independentVariableNames, areCoefficientsInAerodynamicFrame,
            areCoefficientsInNegativeAxisDirection ),
        independentVariables_( independentVariables ),
        forceCoefficients_( forceCoefficients ),
        momentCoefficients_( momentCoefficients ){ }

    std::vector< std::vector< double > > getIndependentVariables( )
    {
        return independentVariables_;
    }

    boost::multi_array< Eigen::Vector3d, NumberOfDimensions > getForceCoefficients( )
    {
        return forceCoefficients_;
    }

    boost::multi_array< Eigen::Vector3d, NumberOfDimensions > getMomentCoefficients( )
    {
        return momentCoefficients_;
    }

private:

    std::vector< std::vector< double > > independentVariables_;
    boost::multi_array< Eigen::Vector3d, NumberOfDimensions > forceCoefficients_;
    boost::multi_array< Eigen::Vector3d, NumberOfDimensions > momentCoefficients_;


};


boost::shared_ptr< aerodynamics::AerodynamicCoefficientInterface > createAerodynamicCoefficientInterface(
        const boost::shared_ptr< AerodynamicCoefficientSettings > shapeSettings,
        const std::string& body );

boost::shared_ptr< aerodynamics::FlightConditions > createFlightConditions(
        const boost::shared_ptr< Body > bodyWithFlightConditions,
        const boost::shared_ptr< Body > centralBody,
        const boost::function< double( ) > angleOfAttackFunction =
        boost::lambda::constant ( 0.0 ),
        const boost::function< double( ) > angleOfSideslipFunction =
        boost::lambda::constant ( 0.0 ),
        const boost::function< double( ) > bankAngleFunction =
        boost::lambda::constant ( 0.0 ) );


}

}
#endif // CREATEACCELERATIONMODELS_H
