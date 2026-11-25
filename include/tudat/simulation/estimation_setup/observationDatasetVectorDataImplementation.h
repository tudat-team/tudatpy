template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::setLinkEndReferencePoint(
        const std::string& bodyName,
        const std::string& referencePointName,
        const LinkEndType linkEndType,
        const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition )
{
    const std::vector< unsigned int > matchingSetIds = getSetIds( condition );
    for( const unsigned int setId : std::set< unsigned int >( matchingSetIds.begin( ), matchingSetIds.end( ) ) )
    {
        std::map< LinkEndType, LinkEndId > linkEnds = getLinkDefinition( setMetadata_.at( setId ).linkDefinitionId_ ).linkEnds_;
        typename std::map< LinkEndType, LinkEndId >::iterator linkEndIterator = linkEnds.find( linkEndType );
        if( linkEndIterator == linkEnds.end( ) || linkEndIterator->second.bodyName_ != bodyName ||
            linkEndIterator->second.getReferencePointName( ) == referencePointName )
        {
            continue;
        }

        linkEndIterator->second = LinkEndId( linkEndIterator->second.bodyName_, referencePointName );
        resetLinkDefinitionForSet( setId, LinkDefinition( linkEnds ) );
    }
}
template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::resetDependentVariableBookkeepingForSet(
        const unsigned int setId,
        const std::shared_ptr< simulation_setup::ObservationDependentVariableBookkeeping >& dependentVariableBookkeeping )
{
    const auto& metadata = getObservationSetMetadata( setId );
    if( dependentVariableBookkeeping )
    {
        if( dependentVariableBookkeeping->getObservableType( ) != metadata.observableType_ ||
            !( dependentVariableBookkeeping->getLinkEnds( ) == getLinkDefinition( metadata.linkDefinitionId_ ) ) )
        {
            throw std::runtime_error( "Dependent-variable bookkeeping is incompatible with the observation set." );
        }
        for( const auto id : getObservationIdsForSet( setId ) )
        {
            const auto size = getDependentVariables( id ).size( );
            if( size != 0 && size != dependentVariableBookkeeping->getTotalDependentVariableSize( ) )
            {
                throw std::runtime_error( "Dependent-variable layout does not match the stored values." );
            }
        }
    }
    setMetadata_.at( setId ).dependentVariableLayoutId_ = registerEntry( dependentVariableLayoutRegistry_, dependentVariableBookkeeping );
    ++structuralVersion_;
}
template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::vector< unsigned int > ObservationDataset< ObservationScalarType, TimeType, Dummy >::getAllObservationIds( ) const
{
    std::vector< unsigned int > result;
    result.reserve( observationRows_.size( ) );
    for( const auto& row : observationRows_ )
    {
        result.push_back( row.observationId_ );
    }
    return result;
}
