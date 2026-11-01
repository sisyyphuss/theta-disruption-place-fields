"""Session discovery and NWB loading for Dandiset 000059 (Petersen & Buzsaki 2020).

Mirrors the loading pattern from notebooks/01-loading-nwb-data.py and
notebooks/02a-spike-stimulus-analysis.py in the ncbs-neural-circuits-navigation
course repo, generalized to iterate over *all* sessions for a subject rather
than a single hard-coded session.
"""
from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime

import h5py
import numpy as np
from dandi.dandiapi import DandiAPIClient
from pynwb import NWBHDF5IO
from pynwb.file import NWBFile
from remfile import File as RemoteFile

DANDISET_ID = "000059"
DANDISET_VERSION = "0.230907.2101"


@dataclass
class SessionHandle:
    subject: str
    path: str
    nwbfile: NWBFile
    io: NWBHDF5IO

    def close(self) -> None:
        self.io.close()


def list_processed_sessions(subject: str) -> list[str]:
    """Return sorted asset paths for a subject's processed behavior+ecephys files."""
    with DandiAPIClient() as client:
        dandiset = client.get_dandiset(DANDISET_ID, DANDISET_VERSION)
        assets = dandiset.get_assets_with_path_prefix(f"sub-{subject}")
        paths = [a.path for a in assets if "processed_behavior" in a.path]
    return sorted(paths)


def open_session(path: str) -> SessionHandle:
    """Stream-open one NWB session by its DANDI asset path."""
    subject = path.split("/")[0].removeprefix("sub-")
    with DandiAPIClient() as client:
        dandiset = client.get_dandiset(DANDISET_ID, DANDISET_VERSION)
        asset = dandiset.get_asset_by_path(path)
        url = asset.get_content_url(follow_redirects=1, strip_query=True)

    remote_file = RemoteFile(url)
    h5_file = h5py.File(remote_file, "r")
    io = NWBHDF5IO(file=h5_file, load_namespaces=True)
    nwbfile = io.read()
    return SessionHandle(subject=subject, path=path, nwbfile=nwbfile, io=io)


def session_datetime(session: SessionHandle) -> datetime:
    return session.nwbfile.session_start_time


def get_good_units(session: SessionHandle):
    """Return (units_df, good_unit_indices) filtered to quality in {good, good2}.

    Returns (None, empty array) if this session has no units table at all —
    observed to happen for a subset of 000059 sessions.
    """
    if session.nwbfile.units is None:
        return None, np.array([], dtype=int)
    units_df = session.nwbfile.units.to_dataframe()
    good_mask = units_df["quality"].isin(["good", "good2"])
    good_indices = np.where(good_mask)[0]
    return units_df, good_indices


def get_position(session: SessionHandle):
    """Return (x, y, timestamps) from the behavior module's position interface."""
    behavior_module = session.nwbfile.processing["behavior"]
    position_interface = next(
        interface
        for name, interface in behavior_module.data_interfaces.items()
        if "position" in name.lower()
    )
    spatial_series = next(iter(position_interface.spatial_series.values()))
    position_data = spatial_series.data[:]
    timestamps = spatial_series.timestamps[:]
    return position_data[:, 0], position_data[:, 1], timestamps
