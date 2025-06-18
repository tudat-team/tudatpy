#ifndef TUDAT_TABULATED_IONOSPHERE_MODEL_H
#define TUDAT_TABULATED_IONOSPHERE_MODEL_H

#include <memory>
#include <vector>
#include <utility>
#include <stdexcept>
#include <Eigen/Core>

#include "tudat/math/interpolators/multiDimensionalInterpolator.h"
#include "tudat/math/interpolators/createInterpolator.h"

namespace tudat
{

namespace environment
{

//! Base class for ionospheric model.
class IonosphereModel
{
public:
    virtual ~IonosphereModel( ) { }

    //! Get vertical total electron content in TECU (1e16 e-/m²)
    virtual double getVerticalTotalElectronContent( const double latitude,
                                                    const double longitude,
                                                    const double time ) = 0;

    virtual double getReferenceIonosphereHeight() const = 0;

};

//! Class defining a VTEC model using tabulated IONEX data.
class TabulatedIonosphereModel : public IonosphereModel
{
public:

    //! Constructor from grid + interpolator settings
    TabulatedIonosphereModel(
        const std::vector< double >& times,
        const std::vector< double >& latitudes,
        const std::vector< double >& longitudes,
        const boost::multi_array< double, 3 >& tecGrid,
        const std::shared_ptr< interpolators::InterpolatorSettings >& interpolatorSettings )
    {
        std::vector< interpolators::BoundaryInterpolationType > boundaryHandling =
            interpolatorSettings->getBoundaryHandling( );

        // ✅ Default to warning if unspecified
        if ( boundaryHandling.empty( ) )
        {
            boundaryHandling = std::vector< interpolators::BoundaryInterpolationType >(
                3, interpolators::use_boundary_value_with_warning );
        }

        tecInterpolator_ = interpolators::createMultiDimensionalInterpolator< double, double, 3 >(
            { times, latitudes, longitudes },
            tecGrid,
            interpolatorSettings );
    }

    //! Constructor directly from interpolator
    explicit TabulatedIonosphereModel(
        const std::shared_ptr< interpolators::MultiDimensionalInterpolator< double, double, 3 > >& interpolator,
        const double referenceIonosphereHeight )
        : tecInterpolator_( interpolator ), referenceIonosphereHeight_( referenceIonosphereHeight )
    { }

    //! Get vertical TEC at given lat [deg], lon [deg], and time [s since J2000]
    double getVerticalTotalElectronContent( const double latitude,
                                            const double longitude,
                                            const double time ) override
    {
        std::vector< double > query = { time, latitude, longitude };
        lastQuery_ = query;
        return tecInterpolator_->interpolate( query );
    }

    //! For debugging — return last query
    std::vector< double > getLastQuery( ) const
    {
        return lastQuery_;
    }

    double getReferenceIonosphereHeight( ) const override
    {
        return referenceIonosphereHeight_;
    }

private:

    std::shared_ptr< interpolators::MultiDimensionalInterpolator< double, double, 3 > > tecInterpolator_;
    std::vector< double > lastQuery_;
    double referenceIonosphereHeight_;
};

} // namespace environment
} // namespace tudat

#endif // TUDAT_TABULATED_IONOSPHERE_MODEL_H
