"""
Ramp converter
"""

from . import Converter
from trk234 import bands

from pandas import DataFrame, concat


class RampConverter(Converter):
    def extract(self, sfdu_list):
        # Filter SFDU objects that represent ramp data.
        # - Ramp format_code == 9
        # - Only keep decoded ones
        # - Only keep ramp types 1, 4, 5
        ramp_sfdu_list = [sfdu for sfdu in sfdu_list if sfdu.pri_chdo.format_code == 9]
        for sfdu in ramp_sfdu_list:
            sfdu.decode(sfdu.binarydata, label=False, agg_chdo=False, pri_chdo=False)
        ramp_sfdu_list = [sfdu for sfdu in ramp_sfdu_list if sfdu.is_decoded]
        ramp_sfdu = [
            sfdu
            for sfdu in ramp_sfdu_list
            if sfdu.is_decoded and (sfdu.trk_chdo.ramp_type in [1, 4, 5])
        ]

        data = {
            "epoch": [sfdu.timestamp() for sfdu in ramp_sfdu],
            "station": ["DSS-" + str(sfdu.sec_chdo.ul_dss_id) for sfdu in ramp_sfdu],
            "band": [bands[sfdu.sec_chdo.ul_band] for sfdu in ramp_sfdu],
            "type": [sfdu.trk_chdo.ramp_type for sfdu in ramp_sfdu],
            "freq": [sfdu.trk_chdo.ramp_freq for sfdu in ramp_sfdu],
            "rate": [sfdu.trk_chdo.ramp_rate for sfdu in ramp_sfdu],
        }

        return DataFrame(data)

    def process(self, ramp_df):
        """
        Concatenate ramp records and merge them per station following these rules:

        1. Each ramp DataFrame has an 'epoch' column and a 'type' column.
           A row with type==1 indicates the start of a new ramp;
           rows with type==4 or type==5 mark ramp end events.
        2. If a ramp has no end time (i.e. a start event is followed by a new start event),
           update the previous ramp's end time to the new ramp's start time.
        3. If a new start event is a continuation of the previous ramp, meaning that its
           frequency equals the frequency extrapolated from the previous ramp's (freq, rate)
           and its rate matches the previous ramp's rate (within a tolerance),
           then update the previous ramp's end time and do not add a new record.
        4. If two ramp intervals (with fully specified end times) overlap,
           the new ramp takes priority and the previous ramp's end time is updated
           to the new ramp's start time.

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
        tolerance = 1e-6  # Precision from TNF data
        merged_dfs = []
        for station, group in ramp_df.groupby("station"):
            group = group.sort_values("epoch").reset_index(drop=True)
            merged_intervals = []
            current_interval = None
            for _, row in group.iterrows():
                event_type = row["type"]
                event_time = row["epoch"]
                event_freq = row["freq"]
                event_rate = row["rate"]
                if event_type == 1:
                    if current_interval is None:
                        current_interval = row.to_dict()
                        current_interval["end_time"] = None
                    else:
                        delta_t = (
                            event_time - current_interval["epoch"]
                        ).total_seconds()
                        expected_freq = (
                            current_interval["freq"]
                            + current_interval["rate"] * delta_t
                        )
                        if (
                            abs(event_freq - expected_freq) <= tolerance
                            and abs(event_rate - current_interval["rate"]) <= tolerance
                        ):
                            current_interval["end_time"] = event_time
                            continue
                        else:
                            current_interval["end_time"] = event_time
                            merged_intervals.append(current_interval)
                            current_interval = row.to_dict()
                            current_interval["end_time"] = None
                elif event_type in [4, 5]:
                    if current_interval is not None:
                        current_interval["end_time"] = event_time
                        merged_intervals.append(current_interval)
                        current_interval = None
            if current_interval is not None:
                merged_intervals.append(current_interval)
            merged_intervals = sorted(merged_intervals, key=lambda x: x["epoch"])
            final_intervals = []
            for interval in merged_intervals:
                if final_intervals:
                    last = final_intervals[-1]
                    if (
                        last.get("end_time") is not None
                        and interval["epoch"] < last["end_time"]
                    ):
                        last["end_time"] = interval["epoch"]
                final_intervals.append(interval)
            merged_dfs.append(DataFrame(final_intervals))
        merged_df = concat(merged_dfs, ignore_index=True)
        cols = list(merged_df.columns)
        if "epoch" in cols and "end_time" in cols:
            remaining = [c for c in cols if c not in ["epoch", "end_time"]]
            merged_df = merged_df[["epoch", "end_time"] + remaining]
        merged_df = merged_df.rename(columns={"epoch": "start_time"})

        return merged_df
