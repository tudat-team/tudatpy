//! TUDAT libraries.
#include <Tudat/Mathematics/BasicMathematics/basicFunction.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h>

/*!
* Function that is used in the root finder. This function finds the elevator deflection for trimmed flight.
* Function : CM( angleOfAttack , machNumber ) + dCM( angleOfAttack , machNumber , elevatorDeflection ) = 0
*/

namespace tudat
{

namespace aerodynamics
{

class HeatFluxFunction: public tudat::basic_mathematics::Function<double,double>
{
    public:
        //! Constructor.
        HeatFluxFunction( double currentAirspeed , double currentDensity , double noseRadius )
        {
            currentAirspeed_ = currentAirspeed;
            currentDensity_  = currentDensity;
            noseRadius_      = noseRadius;
        }

        //! Destructor.
        virtual ~HeatFluxFunction(){}

        double evaluate( const double inputValue )
        {
            return heatFluxConstant_ * sqrt( currentDensity_ * std::pow(currentAirspeed_,2.0) / noseRadius_ )
                    * ( 0.5 * std::pow(currentAirspeed_,2.0) + 1004 * ( currentFreestreamTemperature_ - inputValue ) )
                    - wallEmissivity_*stefanBoltzmannConstant_*std::pow(inputValue,4);
        }

        // FUNCTION NOT IMPLEMENTED
        double computeDerivative( const unsigned int order, const double independentVariable )
        {
            return TUDAT_NAN;
        }

        // FUNCTION NOT IMPLEMENTED
        double computeDefiniteIntegral( const unsigned int order, const double lowerBound, const double upperbound )
        {
            return TUDAT_NAN;
        }

        void setWallEmissivity( double wallEmissivity ) { wallEmissivity_ = wallEmissivity; }

        void setAdiabaticWallTemperature( double adiabaticWallTemperature ) { adiabaticWallTemperature_ = adiabaticWallTemperature; }

        void setCurrentFreestreamTemperature( double currentFreestreamTemperature ) { currentFreestreamTemperature_ = currentFreestreamTemperature; }

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
        const double stefanBoltzmannConstant_ = tudat::physical_constants::STEFAN_BOLTZMANN_CONSTANT;
};

} //namespace_aerodynamics
} //namespace_tudat
