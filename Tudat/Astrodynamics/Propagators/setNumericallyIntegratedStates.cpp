#include "Tudat/Astrodynamics/Propagators/setNumericallyIntegratedStates.h"
#include "Tudat/Mathematics/Interpolators/lagrangeInterpolator.h"

namespace tudat
{

namespace propagators
{

void checkTranslationalStatesFeasibility(
        const std::vector< std::string >& bodiesToIntegrate,
        const simulation_setup::NamedBodyMap& bodyMap )
{
    for( simulation_setup::NamedBodyMap::const_iterator bodyIterator = bodyMap.begin( ); bodyIterator != bodyMap.end( ); bodyIterator++ )
    {
        if( std::find( bodiesToIntegrate.begin( ), bodiesToIntegrate.end( ), bodyIterator->first ) == bodiesToIntegrate.end( ) )
        {
            std::string ephemerisOrigin = bodyIterator->second->getEphemeris( )->getReferenceFrameOrigin( );
            if( std::find( bodiesToIntegrate.begin( ), bodiesToIntegrate.end( ), ephemerisOrigin ) != bodiesToIntegrate.end( ) )
            {
                std::cerr<<"Warning, found non-integrated body with an integrated body as ephemeris origin"<<" "<<
                           bodyIterator->second->getEphemeris( )->getReferenceFrameOrigin( )<<" "<<bodyIterator->first<<std::endl;
                //boost::throw_exception(
                //            boost::enable_error_info(
                //                std::runtime_error( "Error, found non-integrated body with an integrated body as ephemeris origin" ) ) );
            }
        }
    }
}

template< >
boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Matrix< double, 6, 1 > > >
createStateInterpolator( const std::map< double, Eigen::Matrix< double, 6, 1 > >& stateMap )
{
    return boost::make_shared< interpolators::LagrangeInterpolator< double, Eigen::Matrix< double, 6, 1 > > >( stateMap, 6 );
}

template< >
boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Matrix< long double, 6, 1 > > >
createStateInterpolator( const std::map< double, Eigen::Matrix< long double, 6, 1 > >& stateMap )
{
    return boost::make_shared< interpolators::LagrangeInterpolator< double, Eigen::Matrix< long double, 6, 1 > > >( stateMap, 6 );
}

std::vector< std::string > determineEphemerisUpdateorder( std::vector< std::string > integratedBodies,
                                                          std::vector< std::string > centralBodies,
                                                          std::vector< std::string > ephemerisOrigins )
{

    std::vector< std::string > updateOrder;

    bool isFinished = 0;
    int currentIndex = 0;
    std::vector< std::string >::const_iterator centralBodyIterator;
    std::vector< std::string >::const_iterator ephemerisOriginIterator;
    int counter = 0;

    while( !isFinished )
    {
        centralBodyIterator = std::find( integratedBodies.begin( ), integratedBodies.end( ), centralBodies.at( currentIndex ) );
        ephemerisOriginIterator = std::find( integratedBodies.begin( ), integratedBodies.end( ), ephemerisOrigins.at( currentIndex ) );

        if( centralBodyIterator == integratedBodies.end( ) && ephemerisOriginIterator == integratedBodies.end( ) )
        {
            updateOrder.push_back( integratedBodies.at( currentIndex ) );
            integratedBodies.erase( integratedBodies.begin( ) + currentIndex );
            centralBodies.erase( centralBodies.begin( ) + currentIndex );
            ephemerisOrigins.erase( ephemerisOrigins.begin( ) + currentIndex );

            currentIndex = 0;


            if( integratedBodies.size( ) == 0 )
            {
                isFinished = 1;
            }
        }
        else
        {
            if( centralBodyIterator != integratedBodies.end( ) && ephemerisOriginIterator != integratedBodies.end( ) )
            {

                currentIndex = std::min(
                            std::distance( integratedBodies.begin( ), std::find( integratedBodies.begin( ), integratedBodies.end( ),
                                                                                 centralBodies.at( currentIndex ) ) ),
                            std::distance( integratedBodies.begin( ), std::find( integratedBodies.begin( ), integratedBodies.end( ),
                                                                                 ephemerisOrigins.at( currentIndex ) ) ));

            }
            else if( centralBodyIterator != integratedBodies.end( ) )
            {
                currentIndex = std::distance( integratedBodies.begin( ), std::find( integratedBodies.begin( ), integratedBodies.end( ),
                                                                                    centralBodies.at( currentIndex ) ) );
            }
            else if( ephemerisOriginIterator != integratedBodies.end( ) )
            {
                currentIndex = std::distance( integratedBodies.begin( ), std::find( integratedBodies.begin( ), integratedBodies.end( ),
                                                                                    ephemerisOrigins.at( currentIndex ) ) );
            }
        }
        counter++;
        if( counter > 10000 )
        {
            std::cerr<<"Warning, ephemeris update order determination now at iteration "<<counter<<std::endl;
        }
    }

    return updateOrder;
}

}

}
