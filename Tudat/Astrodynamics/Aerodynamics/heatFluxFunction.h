
#include <Tudat/Mathematics/BasicMathematics/basicFunction.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h>

namespace tudat
{

namespace aerodynamics
{

class EquilibriumTemperatureFunction: public tudat::basic_mathematics::Function<double,double>
{
public:
    //! Constructor.
    EquilibriumTemperatureFunction(
            const boost::function< double( const double ) > heatTransferFunction,
            const double wallEmissivity,
            double adiabaticWallTemperature ):
       heatTransferFunction_( heatTransferFunction ), wallEmissivity_( wallEmissivity ),
       adiabaticWallTemperature_( adiabaticWallTemperature ){ }


    //! Destructor.
    ~EquilibriumTemperatureFunction(){}

    double evaluate( const double currentWallTemperature )
    {
        return heatTransferFunction_( currentWallTemperature )
                - wallEmissivity_* tudat::physical_constants::STEFAN_BOLTZMANN_CONSTANT * std::pow( currentWallTemperature, 4.0 );
    }

    // FUNCTION NOT IMPLEMENTED
    double computeDerivative( const unsigned int order, const double independentVariable )
    {
        throw std::runtime_error( "Error, derivative of heat flux not defined" );
        return TUDAT_NAN;
    }

    // FUNCTION NOT IMPLEMENTED
    double computeDefiniteIntegral( const unsigned int order, const double lowerBound, const double upperbound )
    {
        throw std::runtime_error( "Error, integrall of heat flux not defined" );
        return TUDAT_NAN;
    }

    double getLowerBound( ) { return 0.0; }

    double getUpperBound( ) { return adiabaticWallTemperature_*0.25; }

    double getInitialGuess( ) { return adiabaticWallTemperature_*0.01; }

protected:

private:
    // Private variables.
    boost::function< double( const double ) > heatTransferFunction_;

    const double wallEmissivity_;

    double adiabaticWallTemperature_;
};

} //namespace_aerodynamics
} //namespace_tudat
