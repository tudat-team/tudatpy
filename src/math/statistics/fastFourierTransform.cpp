#include <iostream>

#include "tudat/math/statistics/fastFourierTransform.h"

namespace tudat
{

namespace fftw_interface
{

std::vector< std::complex< double > > performFftOfRealData(
        const std::vector< double >& realTimeDomainData )
{
    // Convert input data to raw array.
    double* rawArrayTimeDomainData = new double[ realTimeDomainData.size( ) ];
    for( unsigned int i = 0; i < realTimeDomainData.size( ); i++ )
    {
        rawArrayTimeDomainData[ i ] = realTimeDomainData[ i ];
    }

    // Perform fft of real data using raw array.
    fftw_complex* rawArrayFrequencyDomainData = performRawFftOfRealData(
                rawArrayTimeDomainData, realTimeDomainData.size( ) );

    // Copy raw array to vector of complex doubles.
    std::vector< std::complex< double > > frequencyDomainData;
    frequencyDomainData.resize( realTimeDomainData.size( ) / 2 + 1 );
    for( unsigned int i = 0; i < realTimeDomainData.size( ) / 2 + 1; i++ )
    {
        frequencyDomainData[ i ] = std::complex< double >(
                    rawArrayFrequencyDomainData[ i ][ 0 ], rawArrayFrequencyDomainData[ i ][ 1 ] );
    }

    // Deallocate raw arrays.
    fftw_free( rawArrayFrequencyDomainData );
    delete[ ] rawArrayTimeDomainData;

    // Return Fourier transform of real data.
    return frequencyDomainData;
}

fftw_complex* performRawFftOfRealData( double* realTimeDomainData, const int numberOfDataPoints )
{
    // Create and allocate raw array for fourier transform return data.
    fftw_complex *complexFrequencyDomainData;
    complexFrequencyDomainData = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ( std::ceil( numberOfDataPoints ) ) );

    // Create fftw plan for transform of real data.
    fftw_plan p;
    p = fftw_plan_dft_r2c_1d( numberOfDataPoints, realTimeDomainData, complexFrequencyDomainData, FFTW_ESTIMATE );

    // Perform Fourier transform;
    fftw_execute(p);

    // De-allocate plan.
    fftw_destroy_plan(p);

    return complexFrequencyDomainData;
}


std::vector< double > performInverseFftToRealData( const std::vector< std::complex< double > >& frequencyDomainData )
{
    // Allocate raw complex data type.
    fftw_complex* rawArrayFrequencyDomainData;
    rawArrayFrequencyDomainData = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * frequencyDomainData.size( ) );

    // Fill values of raw frequency domain array.
    for( unsigned int i = 0; i < frequencyDomainData.size( ); i++ )
    {
        rawArrayFrequencyDomainData[ i ][ 0 ] = frequencyDomainData[ i ].real( );
        rawArrayFrequencyDomainData[ i ][ 1 ] = frequencyDomainData[ i ].imag( );
    }

    int timeDomainDataSize = 2 * frequencyDomainData.size( ) - 2;

    // Perform fft to obtain real data
    double* rawArrayTimeDomainData = performRawInverseFFtToRealData(
                rawArrayFrequencyDomainData, timeDomainDataSize );

    // Copy raw array data into vector.
    std::vector< double > timeDomainData;
    timeDomainData.resize( timeDomainDataSize );
    for( int i = 0; i < timeDomainDataSize; i++ )
    {
        timeDomainData[ i ] = rawArrayTimeDomainData[ i ] / static_cast< double >( timeDomainDataSize );
    }

    // Deallocate raw arrays.
    fftw_free( rawArrayFrequencyDomainData );
    delete[ ] rawArrayTimeDomainData;

    return timeDomainData;
}


double* performRawInverseFFtToRealData( fftw_complex* frequencyDomainData, const int numberOfDataPoints )
{
    // Create and allocate real time doman output data.
    double *realTimeDomainData;
    realTimeDomainData = new double[ numberOfDataPoints ];

    // Create fourier transform plan.
    fftw_plan p;
    p = fftw_plan_dft_c2r_1d( numberOfDataPoints, frequencyDomainData, realTimeDomainData, FFTW_ESTIMATE );

    // Perform fourier transform
    fftw_execute(p);

    // Delete plan.
    fftw_destroy_plan(p);

    // Return time domain data.
    return realTimeDomainData;
}


}

}
