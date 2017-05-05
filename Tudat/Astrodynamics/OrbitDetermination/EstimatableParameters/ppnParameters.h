#ifndef TUDAT_PPNPARAMETERS_H
#define TUDAT_PPNPARAMETERS_H

#include "Tudat/Astrodynamics/Relativity/metric.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"

using std::string;

namespace tudat
{

namespace estimatable_parameters
{

class PPNParameterGamma: public EstimatableParameter< double >
{

public:
    PPNParameterGamma( const boost::shared_ptr< relativity::PPNParameterSet > ppnParameterSet ):
        EstimatableParameter< double >( ppn_parameter_gamma, "global_metric" ),
      ppnParameterSet_( ppnParameterSet ){ }


    ~PPNParameterGamma( ) { }

    double getParameterValue( )
    {
        return ppnParameterSet_->getParameterGamma( );
    }

    void setParameterValue( double parameterValue )
    {
        ppnParameterSet_->setParameterGamma( parameterValue );
    }

    int getParameterSize( )
    {
        return 1;
    }

protected:

private:
    boost::shared_ptr< relativity::PPNParameterSet > ppnParameterSet_;

};

class PPNParameterBeta: public EstimatableParameter< double >
{

public:
    PPNParameterBeta( const boost::shared_ptr< relativity::PPNParameterSet > ppnParameterSet ):
        EstimatableParameter< double >( ppn_parameter_beta, "global_metric" ),
      ppnParameterSet_( ppnParameterSet ){ }

    ~PPNParameterBeta( ) { }

    double getParameterValue( )
    {
        return ppnParameterSet_->getParameterBeta( );
    }

    void setParameterValue( double parameterValue )
    {
        ppnParameterSet_->setParameterBeta( parameterValue );
    }

    int getParameterSize( )
    {
        return 1;
    }

protected:

private:
    boost::shared_ptr< relativity::PPNParameterSet > ppnParameterSet_;

};

}

}

#endif // TUDAT_PPNPARAMETERS_H
