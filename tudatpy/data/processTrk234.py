# %%
from tudatpy.numerical_simulation import (
    environment_setup,
    estimation_setup,
)  # type:ignore
from tudatpy.numerical_simulation import estimation, environment  # type:ignore
from tudatpy.numerical_simulation.estimation_setup import observation  # type:ignore
from tudatpy.astro import time_conversion
from tudatpy.interface import spice
from tudatpy.data.mission_data_downloader import *

import trk234
import requests
import os
import pandas as pd
import numpy as np

# %%
trkModeDict = {
    0: "Unknown",
    1: "1W",
    2: "2W",
    3: "3W",
}

time_scale_converter = time_conversion.default_time_scale_converter()

frequencyBandsDict = {
    "S": observation.FrequencyBands.s_band,
    "X": observation.FrequencyBands.x_band,
    # "K": observation.FrequencyBands.ku_band,
    "Ka": observation.FrequencyBands.ka_band,
}

stationDict = environment.get_approximate_dsn_ground_station_positions()


def getLinkEnds(sfdu):
    """
    Returns the uplink, spacecraft, and downlink IDs for a given SFDU record.
    The secondary CHDO has to be decoded before calling this function.

    Parameters
    ----------
    sfdu : trk234.SFDU
        The SFDU record to extract the link ends from.

    Returns
    -------
    tuple(int, int, int)
        A tuple containing the uplink, spacecraft, and downlink IDs.
        If the uplink or downlink are unknown or not valid, the function returns NaN.
    """
    upLink = (
        sfdu.sec_chdo.vld_ul_stn
        if sfdu.sec_chdo.vld_ul_stn != 0
        else sfdu.sec_chdo.ul_prdx_stn if sfdu.sec_chdo.ul_prdx_stn != 0 else np.nan
    )
    upLink = "DSS-" + str(upLink) if upLink != np.nan else str(upLink)

    # Add minus sign to comply with NAIF convention
    scId = str(-sfdu.sec_chdo.scft_id)
    dlLink = "DSS-" + str(sfdu.sec_chdo.dl_dss_id)
    return (upLink, scId, dlLink)


def getRadioBand(sfdu):
    """
    Returns the uplink and downlink radio bands for a given SFDU record.
    The secondary CHDO has to be decoded before calling this function.

    Parameters
    ----------
    sfdu : trk234.SFDU
        The SFDU record to extract the radio bands from.

    Returns
    -------
    tuple(str, str)
        A tuple containing the uplink and downlink radio bands.
    """
    return (
        trk234.bands[sfdu.sec_chdo.ul_band_dl],
        trk234.bands[sfdu.sec_chdo.vld_dl_band],
    )


def getTrackingMode(sfdu):
    """
    Returns the tracking mode for a given SFDU record.
    The secondary CHDO has to be decoded before calling this function.

    Parameters
    ----------
    sfdu : trk234.SFDU
        The SFDU record to extract the tracking mode from.

    Returns
    -------
    str
        The tracking mode of the SFDU record.
    """
    trkMode = (
        trkModeDict[sfdu.sec_chdo.vld_dop_mode]
        if sfdu.sec_chdo.vld_dop_mode != 0
        else sfdu.tracking_mode()
    )
    return trkMode


def getTurnAroundRatio(sfdu):
    """
    Calculate the turn-around ratio from the given SFDU record.
    The secondary CHDO has to be decoded before calling this function.

    Parameters
    ----------
    sfdu : trk234.SFDU
        The SFDU record to extract the turn-around ratio from.

    Returns
    -------
    float
        The turn-around ratio of the SFDU record. If the turn-around ratio is not valid, the function returns NaN.
    """
    return (
        sfdu.sec_chdo.scft_transpd_turn_num / sfdu.sec_chdo.scft_transpd_turn_den
        if sfdu.sec_chdo.scft_transpd_turn_den != 0.0
        or sfdu.sec_chdo.scft_transpd_turn_num != 0.0
        else np.nan
    )


def getLinkDelays(sfdu):
    """
    Returns the transmit time tag delay, spacecraft transmit delay, and receive time tag delay for a given SFDU record.
    The secondary CHDO has to be decoded before calling this function.

    Parameters
    ----------
    sfdu : trk234.SFDU
        The SFDU record to extract the time tag delays from.

    Returns
    -------
    tuple(float, float, float)
        A tuple containing the transmit time tag delay, spacecraft transmit delay, and receive time tag delay.
        If the values are not valid or not provided, the function returns NaN.
    """
    uplinkDelay = 0.0
    uplinkDelay += (
        sfdu.sec_chdo.transmit_time_tag_delay
        if sfdu.sec_chdo.transmit_time_tag_delay != -1.0
        else 0.0
    )
    uplinkDelay += (
        sfdu.sec_chdo.ul_zheight_corr if sfdu.sec_chdo.ul_zheight_corr != -99.0 else 0.0
    )

    downlinkDelay = 0.0
    downlinkDelay += (
        sfdu.sec_chdo.rcv_time_tag_delay
        if sfdu.sec_chdo.rcv_time_tag_delay != -1.0
        else 0.0
    )
    downlinkDelay += (
        sfdu.sec_chdo.dl_zheight_corr if sfdu.sec_chdo.dl_zheight_corr != -99.0 else 0.0
    )

    scft_transpd_delay = (
        sfdu.sec_chdo.scft_transpd_delay
        if sfdu.sec_chdo.scft_transpd_delay != -1.0
        else 0.0
    )

    return (uplinkDelay, scft_transpd_delay, downlinkDelay)
    # return (
    #     (
    #         sfdu.sec_chdo.transmit_time_tag_delay
    #         if sfdu.sec_chdo.transmit_time_tag_delay != -1.0
    #         else np.nan
    #     ),
    #     (
    #         sfdu.sec_chdo.scft_transpd_delay
    #         if sfdu.sec_chdo.scft_transpd_delay != -1.0
    #         else np.nan
    #     ),
    #     (
    #         sfdu.sec_chdo.rcv_time_tag_delay
    #         if sfdu.sec_chdo.rcv_time_tag_delay != -1.0
    #         else np.nan
    #     ),
    # )


def getArrayDelay(sfdu):
    """
    Returns the time delay added to path by arraying equipment for a given SFDU record.
    The secondary CHDO has to be decoded before calling this function.

    Parameters
    ----------
    sfdu : trk234.SFDU
        The SFDU record to extract the array delay from.

    Returns
    -------
    float
        The array delay of the SFDU record. If the array flag is false, the function returns NaN.

    """
    return sfdu.sec_chdo.array_delay if sfdu.sec_chdo.array_flag != 0.0 else np.nan


# Group by station and apply transformations
def set_start_stop_time_per_station(group):
    station = group.name
    # Sort the group by the "epoch" column
    group = group.sort_values(by="epoch").reset_index(drop=True)

    # Set the end time for each row to the start time of the next ramp
    group["end_time"] = group["epoch"].shift(-1)

    # Drop rows with type 4 and 5 which are "end of ramp" type records
    group = group[~group["type"].isin([4, 5])]

    group["station"] = station

    return group


def extractRamps(sfdu_list):
    """
    Extract ramp data from the given list of SFDU decoded records.
    The label, agg_chdo, pri_chdo CHDOs have to be decoded.

    Parameters
    ----------
    sfdu_list : list(trk234.SFDU)
        The list of SFDU records to extract the ramp data from.

    Returns
    -------
    pd.DataFrame
        A DataFrame containing the ramp data extracted from the SFDU records.
    """
    # Initialize a list to collect the data
    ramp_data = []

    # Loop over all SFDU ramps records
    ramp_sfdu_list = [sfdu for sfdu in sfdu_list if sfdu.pri_chdo.format_code == 9]
    for sfdu in ramp_sfdu_list:
        # Now decode the secondary chdo and tracking chdo - (label, agg_chdo, pri_chdo already decoded)
        sfdu.decode(sfdu.binarydata, label=False, agg_chdo=False, pri_chdo=False)

        # Skip invalid SFDUs
        if not sfdu.is_decoded:
            continue

        # Extract time
        epoch = sfdu.timestamp()

        # Get DSS info
        dss = "DSS-" + str(sfdu.sec_chdo.ul_dss_id)
        band = trk234.bands[sfdu.sec_chdo.ul_band]

        # Get link delays
        uplinkDelay = sfdu.sec_chdo.transmit_time_tag_delay
        uplinkZheight = sfdu.sec_chdo.ul_zheight_corr

        # Get ramp info
        ramp_freq = sfdu.trk_chdo.ramp_freq
        ramp_rate = sfdu.trk_chdo.ramp_rate
        ramp_type = sfdu.trk_chdo.ramp_type

        # Collect the data in a dictionary and add it to the list
        ramp_data.append(
            {
                "epoch": epoch,
                "station": dss,
                "band": band,
                "freq": ramp_freq,
                "rate": ramp_rate,
                "type": ramp_type,
                "uplink_delay": uplinkDelay,
                "uplink_zheight": uplinkZheight,
            }
        )

    # Convert the list of dictionaries to a DataFrame
    rampDf = pd.DataFrame(ramp_data)

    # Drop rows with types 0, 2, 3 and 6
    rampDf = rampDf[~rampDf["type"].isin([0, 2, 3, 6])]

    return rampDf


def extractDoppler(sfdu_list):
    """
    Extract Doppler data from the given list of SFDU decoded records.
    The label, agg_chdo, pri_chdo CHDOs have to be decoded.

    Parameters
    ----------
    sfdu_list : list(trk234.SFDU)
        The list of SFDU records to extract the Doppler data from.

    Returns
    -------
    pd.DataFrame
        A DataFrame containing the Doppler data extracted from the SFDU records.
    """
    # Initialize a list to collect the data
    doppler_data = []

    # Loop over all SFDU doppler records
    doppler_sfdu_list = [sfdu for sfdu in sfdu_list if sfdu.pri_chdo.format_code == 16]
    for sfdu in doppler_sfdu_list:
        # Now decode the secondary chdo and tracking chdo - (label, agg_chdo, pri_chdo already decoded)
        sfdu.decode(sfdu.binarydata, label=False, agg_chdo=False, pri_chdo=False)

        # Skip invalid SFDUs
        if not sfdu.is_decoded:
            continue

        # Collect the data in a dictionary and add it to the list
        doppler_data.append(
            {
                "epoch": sfdu.timestamp(),
                "link_ends": getLinkEnds(sfdu),
                "band": getRadioBand(sfdu),
                "tracking_mode": getTrackingMode(sfdu),
                "turn_around_ratio": getTurnAroundRatio(sfdu),
                "link_delays": getLinkDelays(sfdu),
                "array_delay": getArrayDelay(sfdu),
                "count_time": sfdu.trk_chdo.obs_cnt_time,
                "obs": sfdu.trk_chdo.rcv_carr_obs[0],
            }
        )

    # Convert the list of dictionaries to a DataFrame
    doppler_df = pd.DataFrame(doppler_data)

    return doppler_df


def fromDateTimeToTBD(epoch, station):
    """
    Convert a datetime object in UTC into seconds since J2000 in TDB.

    Parameters
    ----------
    epoch : datetime
        The datetime object to convert.

    Returns
    -------
    float
        The time in seconds since J2000 in TDB.
    """
    if station not in stationDict:
        raise KeyError(
            "Error when processing TNF file, converting time from UTC to TDB: \n"
            + "the position of the ground station {} was not specified.".format(station)
        )

    epoch_utc = time_conversion.datetime_to_tudat(epoch).epoch()
    epoch_tdb = time_scale_converter.convert_time(
        input_scale=time_conversion.utc_scale,
        output_scale=time_conversion.tdb_scale,
        input_value=epoch_utc,
        earth_fixed_position=stationDict[station],
    )
    return epoch_tdb


def merge_ramp_dfs(rampDfs, tolerance=1e-6):
    """
    Concatenate a list of ramp DataFrames and merge them per station following these rules:

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
    rampDfs : list[pd.DataFrame]
        List of ramp DataFrames to concatenate. Each DataFrame must have at least the columns:
        "station", "epoch" (datetime-like), "type" (int), "freq" (float), "rate" (float).
    tolerance : float, optional
        Tolerance used when comparing frequency and rate values.

    Returns
    -------
    pd.DataFrame
        A merged ramp DataFrame with ramp intervals merged per station.
    """
    import numpy as np
    import pandas as pd

    # Concatenate and sort by epoch globally (grouping is per station later)
    all_ramps = pd.concat(rampDfs, ignore_index=True)
    all_ramps.sort_values("epoch", inplace=True)
    all_ramps.reset_index(drop=True, inplace=True)

    merged_dfs = []  # stores merged DataFrames for each station

    # Process each station group separately
    for station, group in all_ramps.groupby("station"):
        group = group.sort_values("epoch").reset_index(drop=True)
        merged_intervals = []  # list of dicts for merged ramp intervals here
        current_interval = None

        for _, row in group.iterrows():
            event_type = row["type"]
            event_time = row["epoch"]
            event_freq = row["freq"]
            event_rate = row["rate"]

            if event_type == 1:  # Start event
                if current_interval is None:
                    # Start a new ramp interval.
                    current_interval = row.to_dict()
                    current_interval["end_time"] = None
                else:
                    # There is an open ramp.
                    delta_t = (event_time - current_interval["epoch"]).total_seconds()
                    # Extrapolated frequency based on previous ramp:
                    expected_freq = (
                        current_interval["freq"] + current_interval["rate"] * delta_t
                    )

                    if np.isclose(
                        event_freq, expected_freq, atol=tolerance
                    ) and np.isclose(
                        event_rate, current_interval["rate"], atol=tolerance
                    ):
                        # This event is a continuation: update the current ramp's end_time.
                        current_interval["end_time"] = event_time
                        # Do not add a new record.
                        continue
                    else:
                        # Not a continuation: finalize current ramp by setting its end time.
                        current_interval["end_time"] = event_time
                        merged_intervals.append(current_interval)
                        # Start a new ramp with the current event.
                        current_interval = row.to_dict()
                        current_interval["end_time"] = None

            elif event_type in [4, 5]:  # End event
                if current_interval is not None:
                    current_interval["end_time"] = event_time
                    merged_intervals.append(current_interval)
                    current_interval = None
                else:
                    # No open ramp; ignore or handle as needed.
                    pass

        # If an open interval remains at the end, add it.
        if current_interval is not None:
            merged_intervals.append(current_interval)

        # Post-process: handle overlapping intervals.
        # If an interval's start is before the previous interval's end, update the previous end.
        merged_intervals = sorted(merged_intervals, key=lambda x: x["epoch"])
        final_intervals = []
        for interval in merged_intervals:
            if final_intervals:
                last = final_intervals[-1]
                if (
                    last["end_time"] is not None
                    and interval["epoch"] < last["end_time"]
                ):
                    last["end_time"] = interval["epoch"]
            final_intervals.append(interval)

        merged_dfs.append(pd.DataFrame(final_intervals))

    merged_df = pd.concat(merged_dfs, ignore_index=True)
    # Optionally, reorder columns to bring 'epoch' and 'end_time' to the front.
    cols = list(merged_df.columns)
    if "epoch" in cols and "end_time" in cols:
        remaining = [c for c in cols if c not in ["epoch", "end_time"]]
        merged_df = merged_df[["epoch", "end_time"] + remaining]

    # Rename 'epoch' column to 'start_time'
    merged_df = merged_df.rename(columns={"epoch": "start_time"})

    return merged_df


def create_ObservationCollection(df, spacecraftName=None):

    observation_set_list = []

    # Extract the unique link ends

    unique_link_ends = df["link_ends"].unique()
    for link_end_tuple in unique_link_ends:

        # Set custom spacecraft name if provided
        if spacecraftName is not None:
            spacecraft = observation.body_origin_link_end_id(spacecraftName)

        else:
            spacecraft = observation.body_origin_link_end_id(link_end_tuple[1])

        link_ends_dict = dict()
        if link_end_tuple[0] == "nan":
            link_ends_dict[observation.transmitter] = spacecraft
            link_ends_dict[observation.receiver] = (
                observation.body_reference_point_link_end_id("Earth", link_end_tuple[2])
            )
        else:
            link_ends_dict[observation.transmitter] = (
                observation.body_reference_point_link_end_id("Earth", link_end_tuple[0])
            )
            link_ends_dict[observation.reflector1] = spacecraft
            link_ends_dict[observation.receiver] = (
                observation.body_reference_point_link_end_id("Earth", link_end_tuple[2])
            )

        link_definition = observation.link_definition(link_ends_dict)

        # Filter out DataFrame rows corresponding to current link ends
        link_df = df[df["link_ends"] == link_end_tuple]

        # Extract unique bands for current link ends
        unique_bands = link_df["band"].unique()
        for band in unique_bands:
            band_df = link_df[link_df["band"] == band]

            # Extract unique time tag delays for this band
            unique_ttd_values = band_df["link_delays"].unique()
            for ttd in unique_ttd_values:
                ttd_df = band_df[band_df["link_delays"] == ttd]

                # Extract unique count times for this band
                unique_count_times = ttd_df["count_time"].unique()
                for count_time in unique_count_times:

                    ct_df = ttd_df[ttd_df["count_time"] == count_time]

                    ancillary_settings = (
                        observation.dsn_n_way_doppler_ancilliary_settings(
                            [frequencyBandsDict[band[0]], frequencyBandsDict[band[1]]],
                            frequencyBandsDict[band[1]],
                            0.0,
                            count_time,
                            ttd,
                        )
                    )
                    # Convert "obs" to a list of one (m,1)-shaped float64 array
                    obs_values = ttd_df["obs"].to_numpy(dtype=float).reshape((-1, 1))

                    # Convert "epoch" column to seconds since reference_epoch
                    epoch_seconds = (
                        ct_df["epoch"]
                        .apply(lambda t: fromDateTimeToTBD(t, link_end_tuple[2]))
                        .tolist()
                    )

                    # Create one single_observation_set per unique link ends, band, and time_tag_delay
                    observation_set = estimation.single_observation_set(
                        observation.dsn_n_way_averaged_doppler,
                        link_definition,
                        obs_values,
                        epoch_seconds,
                        observation.receiver,
                        ancillary_settings,
                    )

                    observation_set_list.append(observation_set)

    observation_collection = estimation.ObservationCollection(observation_set_list)

    return observation_collection


def setTnfInformationInBodies(
    bodies,
    rampDf,
):
    """
    Set the information extracted from the TNF file in the bodies of the simulation.

    Parameters
    ----------
    tnfDf : pd.DataFrame
        The DataFrame containing the TNF information extracted from the file.
    """

    # Convert epochs to seconds since J2000 (UTC)
    rampDf["start_time_seconds"] = rampDf["start_time"].apply(
        lambda x: time_conversion.datetime_to_tudat(x).epoch()
    )
    rampDf["end_time_seconds"] = rampDf["end_time"].apply(
        lambda x: time_conversion.datetime_to_tudat(x).epoch()
    )

    # Set the frequency interpolator for each ground station
    bodyWithGroundStations = bodies.get("Earth")
    for station in rampDf["station"].unique():
        df = rampDf[rampDf["station"] == station]

        frequencyInterpolator = environment.PiecewiseLinearFrequencyInterpolator(
            df["start_time_seconds"].tolist(),
            df["end_time_seconds"].tolist(),
            df["rate"].tolist(),
            df["freq"].tolist(),
        )

        groundStation = bodyWithGroundStations.get_ground_station(station)
        groundStation.set_transmitting_frequency_calculator(frequencyInterpolator)

    spacecraft = bodies.get("MRO")
    spacecraft.system_models.set_default_transponder_turnaround_ratio_function()


def create_observation_collection_from_tnf(tnf_file_paths, bodies, spacecraftName=None):

    rampDfs = []
    dopplerDfs = []
    for tnf_file_path in tnf_file_paths:
        # Read and decode the file
        tnfReader = trk234.Reader(tnf_file_path)
        tnfReader.decode(sec_chdo=False, trk_chdo=False)
        sfduList = tnfReader.sfdu_list

        rampDfs.append(extractRamps(sfduList))
        dopplerDfs.append(extractDoppler(sfduList))

    rampDf = merge_ramp_dfs(rampDfs)
    dopplerDf = pd.concat(dopplerDfs, ignore_index=True)

    dopplerDf = dopplerDf[dopplerDf["tracking_mode"] == "2W"]

    obsColl = create_ObservationCollection(dopplerDf, spacecraftName)

    setTnfInformationInBodies(bodies, rampDf)

    return obsColl, rampDf, dopplerDf
