
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
            double currentAirspeed , double currentDensity , double noseRadius, double wallEmissivity,
            double adiabaticWallTemperature, double currentFreestreamTemperature ):
        currentAirspeed_( currentAirspeed ), currentDensity_( currentDensity ), noseRadius_( noseRadius ),
        wallEmissivity_( wallEmissivity ),adiabaticWallTemperature_( adiabaticWallTemperature ),
        currentFreestreamTemperature_( currentFreestreamTemperature ){ }


    //! Destructor.
    ~EquilibriumTemperatureFunction(){}

    double evaluate( const double currentWallTemperature )
    {
        return heatFluxConstant_ * sqrt( currentDensity_ * std::pow(currentAirspeed_,2.0) / noseRadius_ )
                * ( 0.5 * std::pow(currentAirspeed_,2.0) + 1004.0 * ( currentFreestreamTemperature_ - currentWallTemperature ) )
                - wallEmissivity_* tudat::physical_constants::STEFAN_BOLTZMANN_CONSTANT * std::pow(currentWallTemperature,4);
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

    double getLowerBound( ) { return 0.5; }

    double getUpperBound( ) { return adiabaticWallTemperature_*0.25; }

    double getInitialGuess( ) { return adiabaticWallTemperature_*0.01; }

protected:

private:
    // Private variables.
    double currentAirspeed_;
    double currentDensity_;
    double noseRadius_;
    double wallEmissivity_;
    double adiabaticWallTemperature_;
    double currentFreestreamTemperature_;

    // Private constants.
    const double heatFluxConstant_ = 3.53E-4;
};

} //namespace_aerodynamics
} //namespace_tudat
