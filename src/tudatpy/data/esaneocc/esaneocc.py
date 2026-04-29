from ESANEOCC import neocc
import pandas as pd
class NEOCCQuery:
    """Handles raw data retrieval from ESA NEOCC."""

    @staticmethod
    def get_observational_ephemerides(
            designator: str,
            start: str,
            stop: str,
            step: int,
            unit: str = 'days',
            center: str = "Earth"
    ) -> pd.DataFrame:
        return neocc.query_object(
            designator,
            tab='ephemerides',
            observatory=center,
            start=start,
            stop=stop,
            step=step,
            step_unit=unit
        ).ephemerides


