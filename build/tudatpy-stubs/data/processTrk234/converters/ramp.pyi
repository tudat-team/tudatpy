import typing
from . import Converter as Converter

class RampConverter(Converter):

    def extract(self, sfdu_list):
        ...

    def process(self, ramp_df):
        """
        Concatenate ramp records and merge them per station following these rules:

        1. Each ramp DataFrame has an \'epoch\' column and a \'type\' column.
           A row with type==1 indicates the start of a new ramp;
           rows with type==4 or type==5 mark ramp end events.
        2. If a ramp has no end time (i.e. a start event is followed by a new start event),
           update the previous ramp\'s end time to the new ramp\'s start time.
        3. If a new start event is a continuation of the previous ramp, meaning that its
           frequency equals the frequency extrapolated from the previous ramp\'s (freq, rate)
           and its rate matches the previous ramp\'s rate (within a tolerance),
           then update the previous ramp\'s end time and do not add a new record.
        4. If two ramp intervals (with fully specified end times) overlap,
           the new ramp takes priority and the previous ramp\'s end time is updated
           to the new ramp\'s start time.

        Ramp merging is done separately per station.

        Parameters
        ----------
        rampDf : pd.DataFrame
            Ramp DataFrame to process. Each DataFrame must have at least the columns:
            "station", "epoch" (datetime-like), "type" (int), "freq" (float), "rate" (float).
        tolerance : float, optional
            Tolerance used when comparing frequency and rate values. Default is 1e-6 Hz that is the
            precision from TNF data.

        Returns
        -------
        pd.DataFrame
            A merged ramp DataFrame with ramp intervals merged per station.
        """