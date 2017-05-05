#include <boost/make_shared.hpp>

#include "Tudat/Astrodynamics/Relativity/metric.h"

namespace tudat
{

namespace relativity
{

boost::shared_ptr< PPNParameterSet > ppnParameterSet = boost::make_shared< PPNParameterSet >( 1.0, 1.0 );

}

}
