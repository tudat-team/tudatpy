namespace tudat
{
namespace observation_models
{
template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::vector< double > ObservationDataset< ObservationScalarType, TimeType, Dummy >::getNumericalObservationMetadata(
        const std::string& key,
        const std::vector< unsigned int >& observationIds ) const
{
    const auto metadataIterator = numericalObservationMetadata_.find( key );
    if( metadataIterator == numericalObservationMetadata_.end( ) )
    {
        throw std::out_of_range( "Error when retrieving numerical observation metadata, there is no metadata with key '" + key + "'." );
    }
    std::vector< double > values;
    values.reserve( observationIds.size( ) );
    for( const auto id : observationIds )
    {
        const auto value = metadataIterator->second.find( id );
        if( value == metadataIterator->second.end( ) )
        {
            throw std::out_of_range( "Error when retrieving numerical observation metadata, observation " + std::to_string( id ) +
                                     " has no value for key '" + key + "'." );
        }
        values.push_back( value->second );
    }
    return values;
}
