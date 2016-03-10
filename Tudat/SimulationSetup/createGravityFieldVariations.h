#ifndef TUDAT_CREATEGRAVITYFIELDVARIATIONS_H
#define TUDAT_CREATEGRAVITYFIELDVARIATIONS_H

#include <string>
#include <vector>

#include <boost/assign/list_of.hpp>


#include "Tudat/SimulationSetup/body.h"
#include "Tudat/Astrodynamics/Gravitation/gravityFieldVariations.h"
#include "Tudat/Mathematics/Interpolators/createInterpolator.h"

namespace tudat
{

namespace simulation_setup
{

class ModelInterpolationSettings
{
public:
    ModelInterpolationSettings(
            const double initialTime = 0.0,
            const double finalTime = 0.0,
            const double timeStep = 0.0,
            const boost::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings =
            boost::make_shared< interpolators::LagrangeInterpolatorSettings >( 6 ) ):
        interpolatorSettings_( interpolatorSettings ), initialTime_( initialTime ), finalTime_( finalTime ), timeStep_( timeStep ){ }

    boost::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings_;
    double initialTime_;
    double finalTime_;
    double timeStep_;
};

class GravityFieldVariationSettings
{
public:
    GravityFieldVariationSettings( const gravitation::BodyDeformationTypes bodyDeformationType,
                                   const boost::shared_ptr< ModelInterpolationSettings > interpolatorSettings = NULL ):
        bodyDeformationType_( bodyDeformationType ),
        interpolatorSettings_( interpolatorSettings ){ }

    virtual ~GravityFieldVariationSettings( ){ }

    gravitation::BodyDeformationTypes getBodyDeformationType( ){ return bodyDeformationType_;  }

    boost::shared_ptr< ModelInterpolationSettings > getInterpolatorSettings( ){ return interpolatorSettings_; }

protected:

    gravitation::BodyDeformationTypes bodyDeformationType_;

    boost::shared_ptr< ModelInterpolationSettings > interpolatorSettings_;

};

class BasicSolidBodyGravityFieldVariationSettings: public GravityFieldVariationSettings
{
public:
    BasicSolidBodyGravityFieldVariationSettings(
            const std::vector< std::string > deformingBodies,
            const std::vector< std::vector< std::complex< double > > > loveNumbers,
            const double bodyReferenceRadius,
            const boost::shared_ptr< ModelInterpolationSettings > interpolatorSettings = NULL ):
        GravityFieldVariationSettings( gravitation::basic_solid_body, interpolatorSettings ),
        deformingBodies_( deformingBodies ), loveNumbers_( loveNumbers ), bodyReferenceRadius_( bodyReferenceRadius ){ }

    virtual ~BasicSolidBodyGravityFieldVariationSettings( ){ }

    std::vector< std::string > getDeformingBodies( ){ return deformingBodies_;}
    std::vector< std::vector< std::complex< double > > > getLoveNumbers( ){ return loveNumbers_; }
    double getBodyReferenceRadius( ){ return bodyReferenceRadius_; }

    void resetDeformingBodies( const std::vector< std::string >& deformingBodies ){ deformingBodies_ = deformingBodies; }

protected:
    std::vector< std::string > deformingBodies_;
    std::vector< std::vector< std::complex< double > > > loveNumbers_;
    double bodyReferenceRadius_;

};


class TabulatedGravityFieldVariationSettings: public GravityFieldVariationSettings
{
public:
    TabulatedGravityFieldVariationSettings(
            const std::map< double, Eigen::MatrixXd > cosineCoefficientCorrections,
            const std::map< double, Eigen::MatrixXd > sineCoefficientCorrections,
            const int minimumDegree, const int minimumOrder,
            const boost::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings ):
        GravityFieldVariationSettings(
            gravitation::tabulated_variation, boost::make_shared< ModelInterpolationSettings >(
                TUDAT_NAN, TUDAT_NAN, TUDAT_NAN, interpolatorSettings ) ),
        cosineCoefficientCorrections_( cosineCoefficientCorrections ),
        sineCoefficientCorrections_( sineCoefficientCorrections ),
        minimumDegree_( minimumDegree ), minimumOrder_( minimumOrder ){ }

    std::map< double, Eigen::MatrixXd > getCosineCoefficientCorrections( )
    {
        return cosineCoefficientCorrections_;
    }

    std::map< double, Eigen::MatrixXd > getSineCoefficientCorrections( )
    {
        return sineCoefficientCorrections_;
    }

    int getMinimumDegree( )
    {
        return minimumDegree_;
    }

    int getMinimumOrder( )
    {
        return minimumOrder_;
    }
private:

    std::map< double, Eigen::MatrixXd > cosineCoefficientCorrections_;
    std::map< double, Eigen::MatrixXd > sineCoefficientCorrections_;

    int minimumDegree_;
    int minimumOrder_;

};

boost::shared_ptr< gravitation::GravityFieldVariationsSet > createGravityFieldModelVariationsSet(
        const std::string& body,
        const NamedBodyMap& bodyMap,
        const std::vector< boost::shared_ptr< GravityFieldVariationSettings > >& gravityFieldVariationSettings );

boost::shared_ptr< gravitation::GravityFieldVariations > createGravityFieldVariationsModel(
        const boost::shared_ptr< GravityFieldVariationSettings > gravityFieldVariationSettings,
        const std::string body,
        const NamedBodyMap bodyMap );

}


}
#endif // TUDAT_CREATEGRAVITYFIELDVARIATIONS_H
