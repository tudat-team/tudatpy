#ifndef FASTFOURIERTRANSFORM_H
#define FASTFOURIERTRANSFORM_H

#include <complex>
#include <vector>
#include <map>

#include <fftw3.h>

namespace tudat
{

namespace fftw_interface
{

// Output is vector of length realTimeDomainData.size( ) / 2 + 1, i.e. without the conjugate terms
std::vector< std::complex< double > > performFftOfRealData( const std::vector< double >& realTimeDomainData );

// Output is vector of length numberOfDataPoints (= size of realTimeDomainData), with conjugate terms set to zero.
fftw_complex* performRawFftOfRealData( double* realTimeDomainData, const int numberOfDataPoints );

// Output is vector of length 2 * realTimeDomainData.size( ) / 2 - 2, i.e. input is without conjugate terms
std::vector< double > performInverseFftToRealData( const std::vector< std::complex< double > >& frequencyDomainData );

// Output is vector of length numberOfDataPoints, only first numberOfDataPoints/2 + 1 terms of frequencyDomainData are used
// i.e. assuming other terms are conjugates.
double* performRawInverseFFtToRealData( fftw_complex* frequencyDomainData, const int numberOfDataPoints );


}

}

#endif // FASTFOURIERTRANSFORM_H
