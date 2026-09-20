#!/usr/bin/env python3
"""Reusable COSI BTO spectral-response tools.

Main public API
---------------
``get_bto_response``
    Build one BTO1 or BTO2 response for a fixed spacecraft direction or an
    orientation/time-averaged celestial source.
``get_bto_background``
    Put a constant background spectrum on the response channel grid.
``simulate_bto_spectrum``
    Fold a photon model, add background, draw optional Poisson counts, and
    optionally write a Type-I PHA.
``write_ogip_response``
    Write ARF, RMF, combined RSP, and a metadata sidecar.
``write_ogip_background``
    Write a channel-compatible OGIP background BAK.
``load_bto_response``
    Read a module-written RSP or RMF+ARF pair back into memory.

Configuration objects are :class:`BTOCalibration` and
:class:`BTOResolution`.  Returned data containers are
:class:`BTOResponseProduct`, :class:`BTOBackground`, and
:class:`BTOSpectrumSimulation`.  Everything whose name starts with ``_`` is a
supporting helper, not the main user API.  A normal analysis should call the
main functions above one step at a time; no all-in-one dataset builder is
needed.

This module turns the direction-dependent, unsmeared Geant4/Cosima deposit
response in ``bto_response_order2.h5`` into science-analysis products for one
of the two COSI Background and Transient Observer detectors.  The HDF5 master
response stores

``response/area_matrix_cm2[detector, direction, E_true, E_deposit]``.

The detector axis is read from ``detector/name``; BTO1 and BTO2 are never
identified by an assumed numeric order alone.  They are separate NaI(Tl)
detectors on opposite sides of the payload, so surrounding-mass shadowing
makes their angular responses different.

Direction convention
--------------------
Directions are source line-of-sight vectors in the spacecraft frame.  ``theta``
is colatitude from +Z and ``phi`` is azimuth from +X toward +Y::

    n_source = (sin(theta) cos(phi), sin(theta) sin(phi), cos(theta))

The incident photon momentum is ``-n_source``.  A fixed response can be
requested with ``theta_deg``/``phi_deg`` or a unit ``spacecraft_vector``.
Alternatively, a celestial source and COSI orientation FITS file can be given.
The latter uses the installed COSIpy ``SpacecraftHistory`` transformation and
averages responses with interval livetime (and optionally Earth visibility).
Both modes use the four nearest response directions with inverse angular
distance squared weights.  An exact HEALPix-center query returns that center
with weight one.

Calibration and resolution
--------------------------
Gain supports a monotonic quadratic channel calibration::

    channel_coordinate = a * E_measured_keV**2
                         + b * E_measured_keV + c

where ``a=gain_quadratic_ch_per_keV2``, ``b=gain_slope_ch_per_keV``, and
``c=gain_offset_ch``.  The default ``a=0`` preserves the original linear
calibration.  The inverse conversion uses the physical monotonic branch of
the quadratic equation; users should not invert the polynomial themselves.

The bundled default has 4096 channels, 0.8125 keV/channel, a 30 keV lower
threshold, and a 3000 keV upper saturation.  These values are a configurable
science-analysis PLACEHOLDER, not flight calibration.  BTO1 and BTO2 can be
given separate :class:`BTOCalibration` objects.

The resolution model is

``sigma(E) = sqrt(a0 + a1*E + a2*E**2)`` and
``FWHM(E) = 2*sqrt(2*ln(2))*sigma(E)``.

Nagasawa et al. (Proc. SPIE 13625, DOI 10.1117/12.3063967;
arXiv:2509.00725) report this empirical form and about 65 keV FWHM at 662 keV
for a cylindrical prototype whose geometry differs from the flight detector.
The paper does not publish numerical ``a0/a1/a2`` coefficients.  Consequently,
``Prototype_Resolution`` is explicitly an approximation with
``a0=a1=0`` and constant fractional FWHM 65/662; it is not the paper's full fit
and is not flight calibration.  Arbitrary coefficients can be supplied.

Each deposit-energy bin is redistributed by integrating a Gaussian between
measured-channel edges with the Gaussian CDF.  Probability below threshold,
above saturation, or outside the channel range is lost from the observed
response.  With ``A[Etrue,Edep]`` denoting the interpolated Geant4 area matrix,
the returned quantities are

``rsp_cm2[n_photon,n_channel]``
    Area redistribution after electronics effects.
``arf_cm2[n_photon]``
    ``rsp_cm2.sum(axis=1)``; threshold and saturation losses are included.
``rmf_probability[n_photon,n_channel]``
    ``rsp_cm2 / arf_cm2[:,None]`` for nonzero rows.  It is the conditional
    redistribution of recorded events and each nonzero row sums to one.

Thus ``rsp_cm2 = arf_cm2[:,None] * rmf_probability``.  A standalone RMF stores
the conditional probability and its companion ARF stores the observed area.
A combined RSP stores ``rsp_cm2`` directly.

The HDF5 response does not retain Geant4 interaction-process labels.  The
``photopeak_arf_cm2`` quantity therefore means *full-energy/photopeak effective
area*: deposit energy lies in the corresponding true-energy bin.  It is a
useful proxy for photoelectric/full-energy absorption, but it must not be
interpreted as “the first interaction was photoelectric”.

Examples
--------
Make a fixed spacecraft-direction response without writing files::

    from bto_response import get_bto_response

    response = get_bto_response(
        "response/bto_response_order2.h5", "BTO1",
        theta_deg=35.0, phi_deg=120.0, write_ogip=False,
    )

Make an orientation-averaged response and OGIP files::

    from astropy.coordinates import SkyCoord
    import astropy.units as u

    response = get_bto_response(
        "response/bto_response_order2.h5",
        "BTO1",
        orientation_file="response/20280301_3_month_with_orbital_info.fits",
        source_coord=SkyCoord(l=123.4*u.deg, b=-25*u.deg, frame="galactic"),
        tstart="2028-03-01T12:00:00",
        tstop="2028-03-01T12:10:00",
        earth_occultation=True,
        output_prefix="output/bto1_source",
    )

Fold a photon spectrum and write a fake Type-I PHA::

    from bto_response import (
        get_bto_background, simulate_bto_spectrum, write_ogip_background,
    )

    def powerlaw(E_keV):
        return 1.0e-2 * (E_keV / 100.0)**-2.0

    background = get_bto_background(response)  # time-independent reference BAK
    write_ogip_background(
        background, response, "output/bto1_background.bak", exposure_s=40.0,
    )
    fake = simulate_bto_spectrum(
        response, powerlaw, exposure_s=10.0, background=background,
        poisson=True, random_seed=12345, output_pha="output/bto1_fake.pha",
    )

With the locally installed threeML 2.5 API, separate RMF/ARF products load as::

    from threeML.plugins.OGIPLike import OGIPLike
    plugin = OGIPLike(
        "bto1", observation="output/bto1_fake.pha",
        background="output/bto1_background.bak",
        response="output/bto1_source.rmf",
        arf_file="output/bto1_source.arf", verbose=False,
    )
    plugin.set_active_measurements("30-2000")

Core HDF5 response calculation does not depend on threeML.  Astropy is needed
only for orientation coordinates and FITS I/O; COSIpy is needed only for an
orientation file.
"""

from __future__ import annotations

import hashlib
import json
import math
import warnings
from dataclasses import asdict, dataclass, field, replace
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version as package_version
from pathlib import Path
from typing import Any, Callable, Mapping, Sequence

import h5py
import numpy as np

try:
    from scipy.special import ndtr
except ImportError as exc:  # pragma: no cover - dependency is present in production
    raise ImportError("bto_response requires scipy (scipy.special.ndtr)") from exc


__version__ = "2.2.0"

FWHM_TO_SIGMA = 1.0 / (2.0 * math.sqrt(2.0 * math.log(2.0)))
DEFAULT_REFERENCE_BACKGROUND = (
    Path(__file__).resolve().parents[1] / "response" / "BTO_background_template_10-3000keV.fits"
)


@dataclass(frozen=True)
class BTOCalibration:
    """Quadratic channel calibration and acquisition limits.

    The convention is ``channel = a*E_keV**2 + b*E_keV + c`` with
    ``a=gain_quadratic_ch_per_keV2``, ``b=gain_slope_ch_per_keV``, and
    ``c=gain_offset_ch``.  Set ``a=0`` for a linear calibration.
    ``lower_threshold_keV`` and ``upper_saturation_keV`` act on measured energy.
    The default object is a placeholder and is not flight calibration.
    """

    n_channels: int = 4096
    gain_quadratic_ch_per_keV2: float = 0.0
    gain_slope_ch_per_keV: float = 1.0 / 0.8125
    gain_offset_ch: float = 0.0
    lower_threshold_keV: float = 30.0
    upper_saturation_keV: float | None = 3000.0
    name: str = "BTO science placeholder 4096ch 0.8125keV/ch"
    provenance: str = (
        "Configurable analysis placeholder; not measured flight BTO calibration"
    )

    def __post_init__(self) -> None:
        if self.n_channels <= 0:
            raise ValueError("n_channels must be positive")
        if not np.isfinite(self.gain_quadratic_ch_per_keV2):
            raise ValueError("gain_quadratic_ch_per_keV2 must be finite")
        if not np.isfinite(self.gain_slope_ch_per_keV) or self.gain_slope_ch_per_keV <= 0:
            raise ValueError("gain_slope_ch_per_keV must be finite and positive")
        if not np.isfinite(self.gain_offset_ch):
            raise ValueError("gain_offset_ch must be finite")
        if not np.isfinite(self.lower_threshold_keV):
            raise ValueError("lower_threshold_keV must be finite")
        if self.upper_saturation_keV is not None:
            if not np.isfinite(self.upper_saturation_keV):
                raise ValueError("upper_saturation_keV must be finite or None")
            if self.upper_saturation_keV <= self.lower_threshold_keV:
                raise ValueError("upper_saturation_keV must exceed lower_threshold_keV")
        # The physical inverse branch must exist and remain increasing over all
        # digitizer channel edges.  This rejects a quadratic vertex inside the
        # acquisition range instead of silently creating reversed EBOUNDS.
        channels = np.asarray([0.0, float(self.n_channels)])
        a = float(self.gain_quadratic_ch_per_keV2)
        b = float(self.gain_slope_ch_per_keV)
        c = float(self.gain_offset_ch)
        if np.isclose(a, 0.0, rtol=0.0, atol=1.0e-18):
            energies = (channels - c) / b
        else:
            discriminant = b * b + 4.0 * a * (channels - c)
            if np.any(discriminant < 0.0):
                raise ValueError("quadratic gain has no real inverse over all channel edges")
            roots = np.sqrt(np.clip(discriminant, 0.0, None))
            denominator = b + roots
            if np.any(np.isclose(denominator, 0.0)):
                raise ValueError("quadratic gain inverse is singular over the channel range")
            energies = 2.0 * (channels - c) / denominator
        if not np.all(np.isfinite(energies)) or energies[1] <= energies[0]:
            raise ValueError("quadratic gain must map channel edges to increasing energies")
        if np.any(b + 2.0 * a * energies <= 0.0):
            raise ValueError("quadratic gain derivative must stay positive over the channel range")

    @property
    def keV_per_channel(self) -> float:
        """Local keV/channel at E=0; nonlinear calibrations have no global value."""

        return 1.0 / self.gain_slope_ch_per_keV

    @property
    def channel_edges(self) -> np.ndarray:
        return np.arange(self.n_channels + 1, dtype=np.float64)

    @property
    def measured_energy_edges_keV(self) -> np.ndarray:
        return _channel_to_energy(self.channel_edges, self)


@dataclass(frozen=True)
class BTOResolution:
    """Empirical BTO energy-resolution coefficients.

    ``sigma_keV(E) = sqrt(a0_keV2 + a1_keV*E + a2*E**2)``.
    """

    a0_keV2: float
    a1_keV: float
    a2: float
    name: str
    provenance: str

    def __post_init__(self) -> None:
        values = (self.a0_keV2, self.a1_keV, self.a2)
        if not np.all(np.isfinite(values)):
            raise ValueError("resolution coefficients must be finite")


Prototype_Resolution = BTOResolution(
    a0_keV2=0.0,
    a1_keV=0.0,
    a2=((65.0 / 662.0) * FWHM_TO_SIGMA) ** 2,
    name="Prototype resolution (approximate constant fraction)",
    provenance=(
        "Approximation normalized to the reported 65 keV FWHM at 662 keV; "
        "Nagasawa et al. 2025 prototype; not full-fit coefficients or flight calibration"
    ),
)

# Backward-compatible alias; new examples use Prototype_Resolution.
NAGASAWA2025_APPROXIMATE_RESOLUTION = Prototype_Resolution

NO_SMEARING_RESOLUTION = BTOResolution(
    0.0,
    0.0,
    0.0,
    "No detector smearing",
    "Delta-function redistribution for validation or custom downstream smearing",
)

DEFAULT_CALIBRATION = BTOCalibration()


@dataclass
class BTOResponseProduct:
    """One detector response.

    Shapes are ``arf_cm2[n_photon]``,
    ``rmf_probability[n_photon,n_channel]``, and
    ``rsp_cm2[n_photon,n_channel]``.  ``unsmeared_area_matrix_cm2`` has shape
    ``[n_photon,n_deposit]`` and its row sum is ``unsmeared_arf_cm2``.
    For an OGIP-loaded product, deposit diagnostics are not stored in FITS;
    their legacy array slots are zero placeholders, not physical zero areas.
    Check ``metadata['unsmeared_components_available']`` before plotting them.
    """

    detector_name: str
    photon_energy_edges_keV: np.ndarray
    channel_edges: np.ndarray
    measured_energy_edges_keV: np.ndarray
    arf_cm2: np.ndarray
    rmf_probability: np.ndarray
    rsp_cm2: np.ndarray
    unsmeared_area_matrix_cm2: np.ndarray
    unsmeared_arf_cm2: np.ndarray
    deposit_energy_edges_keV: np.ndarray
    photopeak_arf_cm2: np.ndarray
    calibration: BTOCalibration
    resolution: BTOResolution
    metadata: dict[str, Any] = field(default_factory=dict)
    arf_path: Path | None = None
    rmf_path: Path | None = None
    rsp_path: Path | None = None
    metadata_path: Path | None = None

    def __post_init__(self) -> None:
        n_photon = len(self.photon_energy_edges_keV) - 1
        n_channel = len(self.channel_edges) - 1
        n_deposit = len(self.deposit_energy_edges_keV) - 1
        expected = {
            "arf_cm2": (n_photon,),
            "rmf_probability": (n_photon, n_channel),
            "rsp_cm2": (n_photon, n_channel),
            "unsmeared_area_matrix_cm2": (n_photon, n_deposit),
            "unsmeared_arf_cm2": (n_photon,),
            "photopeak_arf_cm2": (n_photon,),
        }
        for name, shape in expected.items():
            if np.shape(getattr(self, name)) != shape:
                raise ValueError(f"{name} has shape {np.shape(getattr(self, name))}; expected {shape}")


@dataclass
class BTOBackground:
    """Time-independent background rate in counts/s/channel."""

    rate_per_channel: np.ndarray
    measured_energy_edges_keV: np.ndarray
    metadata: dict[str, Any] = field(default_factory=dict)
    path: Path | None = None


@dataclass
class BTOSpectrumSimulation:
    """Folded expectation and optional Poisson realization."""

    expected_source_counts: np.ndarray
    expected_background_counts: np.ndarray
    total_expected_counts: np.ndarray
    simulated_counts: np.ndarray
    integrated_photon_flux: np.ndarray
    exposure_s: float
    channel_edges: np.ndarray
    measured_energy_edges_keV: np.ndarray
    response_metadata: dict[str, Any]
    pha_path: Path | None = None


__all__ = [
    "BTOBackground",
    "BTOCalibration",
    "BTOResolution",
    "BTOResponseProduct",
    "BTOSpectrumSimulation",
    "DEFAULT_CALIBRATION",
    "Prototype_Resolution",
    "NO_SMEARING_RESOLUTION",
    "get_bto_background",
    "get_bto_response",
    "load_bto_response",
    "simulate_bto_spectrum",
    "write_ogip_background",
    "write_ogip_response",
]


# ---------------------------------------------------------------------------
# Main public API -- read this section first
# ---------------------------------------------------------------------------


def get_bto_response(
    response_h5,
    bto_id,
    *,
    theta_deg=None,
    phi_deg=None,
    spacecraft_vector=None,
    orientation_file=None,
    source_coord=None,
    tstart=None,
    tstop=None,
    calibration=None,
    resolution=None,
    threshold_keV=None,
    output_prefix=None,
    write_ogip=True,
    earth_occultation=False,
    overwrite=False,
) -> BTOResponseProduct:
    """Required inputs: master HDF5 ``response_h5`` and ``bto_id``.

    Optional inputs
        Supply exactly one direction mode: (1) ``theta_deg`` + ``phi_deg``,
        (2) ``spacecraft_vector``, or (3) ``orientation_file`` +
        ``source_coord`` + ``tstart`` + ``tstop``.  Calibration, resolution,
        threshold, Earth occultation, and OGIP writing options are optional.
    Output
        :class:`BTOResponseProduct` containing the interpolated unsmeared area
        matrix and detector-folded ARF, RMF, and RSP.  When ``write_ogip=True``,
        ``output_prefix`` is required and the returned path fields are filled.
    """

    return _get_bto_response_impl(
        response_h5,
        bto_id,
        theta_deg=theta_deg,
        phi_deg=phi_deg,
        spacecraft_vector=spacecraft_vector,
        orientation_file=orientation_file,
        source_coord=source_coord,
        tstart=tstart,
        tstop=tstop,
        calibration=calibration,
        resolution=resolution,
        threshold_keV=threshold_keV,
        output_prefix=output_prefix,
        write_ogip=write_ogip,
        earth_occultation=earth_occultation,
        overwrite=overwrite,
    )


def get_bto_background(
    response: BTOResponseProduct,
    *,
    rate_per_channel: float | Sequence[float] | None = None,
    total_rate_hz: float | None = None,
    background_template: str | Path | None = None,
    background_pha: str | Path | None = None,
    background_rmf: str | Path | None = None,
    background_energy_edges_keV: Sequence[float] | None = None,
) -> BTOBackground:
    """Required input: one :class:`BTOResponseProduct`.

    Optional inputs
        Choose at most one background source: ``rate_per_channel``,
        ``total_rate_hz``, or ``background_template``. A template contains its
        own SPECTRUM (RATE or COUNTS+EXPOSURE) and EBOUNDS tables. No native
        energy-grid constants or response file are needed. ``background_pha``
        is retained for PHA compatibility. For a file without EBOUNDS, provide
        ``background_energy_edges_keV`` (n_channel+1 boundaries in keV), or
        ``background_rmf`` for compatibility, but not both. This background
        grid is independent of the new detector's gain/channel calibration.
        With no source argument, the module's 10--3000 keV reference template
        is used. Rates outside the template energy range are zero, not extrapolated.
    Output
        :class:`BTOBackground` in counts/s/channel on exactly the response
        measured-energy grid.  No Poisson realization is made here.
    """

    return _get_bto_background_impl(
        response,
        rate_per_channel=rate_per_channel,
        total_rate_hz=total_rate_hz,
        background_template=background_template,
        background_pha=background_pha,
        background_rmf=background_rmf,
        background_energy_edges_keV=background_energy_edges_keV,
    )


def simulate_bto_spectrum(
    response: BTOResponseProduct,
    photon_model: Any,
    exposure_s: float,
    *,
    background=None,
    poisson=True,
    random_seed=None,
    output_pha=None,
    overwrite=False,
) -> BTOSpectrumSimulation:
    """Required inputs: ``response``, ``photon_model``, and ``exposure_s``.

    Optional inputs
        ``background`` may be a :class:`BTOBackground` or an explicitly
        documented array/mapping.  ``poisson``, ``random_seed``, and optional
        Type-I ``output_pha`` writing control the fake realization.
        With ``poisson=False``, PHA stores the unrounded expectation as RATE,
        POISSERR=False, and explicit Poisson-equivalent STAT_ERR=sqrt(expected
        counts)/exposure (an assumed counting error, not a sampled fluctuation).
    Output
        :class:`BTOSpectrumSimulation` with source, background, total expected,
        and simulated counts per channel.  A callable photon model is treated
        as differential ph cm^-2 s^-1 keV^-1 and integrated in true-energy bins.
    """

    return _simulate_bto_spectrum_impl(
        response,
        photon_model,
        exposure_s,
        background=background,
        poisson=poisson,
        random_seed=random_seed,
        output_pha=output_pha,
        overwrite=overwrite,
    )


def write_ogip_response(
    product: BTOResponseProduct,
    output_prefix,
    *,
    overwrite: bool = False,
) -> BTOResponseProduct:
    """Required inputs: response ``product`` and ``output_prefix``.

    Optional input
        ``overwrite=False`` protects existing products.
    Output
        The same response object with ``arf_path``, ``rmf_path``, ``rsp_path``,
        and ``metadata_path`` filled after writing ARF/RMF/RSP/JSON files.
    """

    return _write_ogip_response_impl(product, output_prefix, overwrite=overwrite)


def write_ogip_background(
    background: BTOBackground,
    response: BTOResponseProduct,
    output_bak,
    *,
    exposure_s: float = 40.0,
    overwrite: bool = False,
) -> BTOBackground:
    """Required inputs: ``background``, matching ``response``, and ``output_bak``.

    Optional inputs
        ``exposure_s`` defines the Poisson-equivalent ``STAT_ERR``; ``overwrite``
        protects existing files.
    Output
        The same :class:`BTOBackground` with ``path`` set to the OGIP BAK.
    """

    return _write_ogip_background_impl(
        background,
        response,
        output_bak,
        exposure_s=exposure_s,
        overwrite=overwrite,
    )


def load_bto_response(response_file, *, arf_file=None) -> BTOResponseProduct:
    """Required input: a module-written ``response_file`` (RSP or RMF).

    Optional input
        ``arf_file`` is required only when ``response_file`` is a standalone RMF.
        Its ENERG_LO/ENERG_HI bins must match the RMF (rtol=1e-7, atol=0,
        allowing float32 ARF versus float64 RMF serialization); otherwise ValueError.
    Output
        :class:`BTOResponseProduct` reconstructed from OGIP.  Geant4 deposit
        matrices are unavailable in OGIP and are represented by zero placeholders.
        A same-prefix ``.response.json`` sidecar, when present, restores source/
        interval metadata; its detector and optional FITS checksum are validated.
        Resolution is already applied: no additional smearing is performed.
    """

    return _load_bto_response_impl(response_file, arf_file=arf_file)


# ---------------------------------------------------------------------------
# Private implementation helpers -- not the main API
# ---------------------------------------------------------------------------


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _jsonable(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Mapping):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    return value


def _dependency_versions() -> dict[str, str | None]:
    versions: dict[str, str | None] = {}
    for distribution in ("numpy", "scipy", "h5py", "astropy", "cosipy", "threeML", "astromodels"):
        try:
            versions[distribution] = package_version(distribution)
        except PackageNotFoundError:
            versions[distribution] = None
    return versions


def _detector_name(bto_id: int | str) -> str:
    if isinstance(bto_id, bool):
        raise ValueError("bto_id must be 1, 2, 'BTO1', or 'BTO2'; bool is invalid")
    if isinstance(bto_id, (int, np.integer)) and int(bto_id) in (1, 2):
        return f"BTO{int(bto_id)}"
    if isinstance(bto_id, str) and bto_id.strip().upper() in ("BTO1", "BTO2"):
        return bto_id.strip().upper()
    raise ValueError(f"invalid bto_id {bto_id!r}; expected 1, 2, 'BTO1', or 'BTO2'")


def _select_for_detector(value: Any, detector_name: str, expected_type: type) -> Any:
    if isinstance(value, Mapping):
        candidates = (detector_name, detector_name.lower(), int(detector_name[-1]))
        for key in candidates:
            if key in value:
                value = value[key]
                break
        else:
            raise KeyError(f"no {detector_name} entry in per-detector mapping")
    if not isinstance(value, expected_type):
        raise TypeError(f"expected {expected_type.__name__} for {detector_name}, got {type(value).__name__}")
    return value


def _energy_to_channel(E_measured_keV: Any, calibration: BTOCalibration = DEFAULT_CALIBRATION) -> np.ndarray:
    """Helper: map energy to channel using ``a*E**2 + b*E + c``."""

    energy = np.asarray(E_measured_keV, dtype=np.float64)
    return (
        calibration.gain_quadratic_ch_per_keV2 * energy**2
        + calibration.gain_slope_ch_per_keV * energy
        + calibration.gain_offset_ch
    )


def _channel_to_energy(channel_coordinate: Any, calibration: BTOCalibration = DEFAULT_CALIBRATION) -> np.ndarray:
    """Helper: invert quadratic gain on its physical monotonic branch."""

    channel = np.asarray(channel_coordinate, dtype=np.float64)
    a = float(calibration.gain_quadratic_ch_per_keV2)
    b = float(calibration.gain_slope_ch_per_keV)
    delta = channel - float(calibration.gain_offset_ch)
    if np.isclose(a, 0.0, rtol=0.0, atol=1.0e-18):
        return delta / b
    discriminant = b * b + 4.0 * a * delta
    if np.any(discriminant < -1.0e-12):
        raise ValueError("channel coordinate is outside the real quadratic-gain branch")
    root = np.sqrt(np.clip(discriminant, 0.0, None))
    denominator = b + root
    if np.any(np.isclose(denominator, 0.0)):
        raise ValueError("quadratic gain inverse is singular at the requested channel")
    return 2.0 * delta / denominator


def _sigma_keV(E_keV: Any, resolution: BTOResolution = Prototype_Resolution) -> np.ndarray:
    """Return Gaussian sigma in keV, rejecting a negative variance model."""

    energy = np.asarray(E_keV, dtype=np.float64)
    variance = resolution.a0_keV2 + resolution.a1_keV * energy + resolution.a2 * energy**2
    tolerance = 1.0e-12 * np.maximum(1.0, np.abs(variance))
    if np.any(variance < -tolerance):
        raise ValueError(f"resolution model {resolution.name!r} produces negative variance")
    return np.sqrt(np.clip(variance, 0.0, None))


def _fwhm_keV(E_keV: Any, resolution: BTOResolution = Prototype_Resolution) -> np.ndarray:
    """Return FWHM in keV."""

    return _sigma_keV(E_keV, resolution) / FWHM_TO_SIGMA


def _vector_from_angles(theta_deg: float, phi_deg: float) -> np.ndarray:
    theta = math.radians(float(theta_deg))
    phi = math.radians(float(phi_deg) % 360.0)
    return np.asarray(
        [math.sin(theta) * math.cos(phi), math.sin(theta) * math.sin(phi), math.cos(theta)],
        dtype=np.float64,
    )


def _normalize_vector(vector: Sequence[float]) -> np.ndarray:
    values = np.asarray(vector, dtype=np.float64)
    if values.shape != (3,) or not np.all(np.isfinite(values)):
        raise ValueError("spacecraft_vector must contain three finite values")
    norm = float(np.linalg.norm(values))
    if norm <= 0:
        raise ValueError("spacecraft_vector must be nonzero")
    return values / norm


def _angles_from_vector(vector: np.ndarray) -> tuple[float, float]:
    vector = _normalize_vector(vector)
    return (
        math.degrees(math.acos(float(np.clip(vector[2], -1.0, 1.0)))),
        math.degrees(math.atan2(float(vector[1]), float(vector[0]))) % 360.0,
    )


def _nearest_direction_weights(
    vector: np.ndarray, response_vectors: np.ndarray, neighbors: int = 4
) -> tuple[np.ndarray, np.ndarray]:
    """Return nearest spherical IDW2 weights, with exact-center preservation."""

    vector = _normalize_vector(vector)
    grid = np.asarray(response_vectors, dtype=np.float64)
    grid /= np.linalg.norm(grid, axis=1)[:, None]
    neighbors = max(1, min(int(neighbors), len(grid)))
    dots = np.clip(grid @ vector, -1.0, 1.0)
    chosen = np.argpartition(-dots, neighbors - 1)[:neighbors]
    angles = np.arccos(dots[chosen])
    exact = np.flatnonzero((1.0 - dots[chosen]) <= 2.0e-14)
    if exact.size:
        weights = np.zeros(neighbors, dtype=np.float64)
        weights[int(exact[0])] = 1.0
    else:
        weights = 1.0 / np.maximum(angles, 1.0e-15) ** 2
        weights /= weights.sum()
    order = np.argsort(angles)
    return chosen[order], weights[order]


def _interpolate_area_matrix(
    detector_matrix: Any,
    response_vectors: np.ndarray,
    vector: np.ndarray,
    *,
    complete: np.ndarray | None = None,
    neighbors: int = 4,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    selected, weights = _nearest_direction_weights(vector, response_vectors, neighbors)
    if complete is not None and not np.all(complete[selected]):
        missing = selected[~complete[selected]].tolist()
        raise RuntimeError(f"nearest response directions are incomplete: {missing}")
    # h5py fancy indices must be monotonically increasing.  Reading the four
    # neighbors individually preserves the angular-weight order for both an
    # h5py Dataset and an in-memory ndarray.
    values = np.stack(
        [np.asarray(detector_matrix[int(index), :, :], dtype=np.float64) for index in selected],
        axis=0,
    )
    averaged = np.einsum("p,pij->ij", weights, values, optimize=True)
    return averaged, selected, weights


def _average_orientation_samples(
    detector_matrix: Any,
    response_vectors: np.ndarray,
    sample_vectors: np.ndarray,
    sample_exposure_s: np.ndarray,
    *,
    complete: np.ndarray | None = None,
    neighbors: int = 4,
) -> tuple[np.ndarray, np.ndarray]:
    """Average interpolated matrices for orientation samples (testable core)."""

    vectors = np.asarray(sample_vectors, dtype=np.float64)
    exposures = np.asarray(sample_exposure_s, dtype=np.float64)
    if vectors.ndim != 2 or vectors.shape[1] != 3 or exposures.shape != (len(vectors),):
        raise ValueError("orientation vectors/exposures have inconsistent shapes")
    if np.any(~np.isfinite(exposures)) or np.any(exposures < 0):
        raise ValueError("orientation exposure weights must be finite and nonnegative")
    total = float(exposures.sum())
    if total <= 0:
        raise RuntimeError("source has zero visible livetime in the requested interval")
    dwell = np.zeros(len(response_vectors), dtype=np.float64)
    for vector, exposure in zip(vectors, exposures):
        if exposure <= 0:
            continue
        selected, weights = _nearest_direction_weights(vector, response_vectors, neighbors)
        if complete is not None and not np.all(complete[selected]):
            missing = selected[~complete[selected]].tolist()
            raise RuntimeError(f"nearest response directions are incomplete: {missing}")
        dwell[selected] += exposure * weights
    fractions = dwell / dwell.sum()
    nonzero = np.flatnonzero(fractions > 0)
    values = np.asarray(detector_matrix[nonzero, :, :], dtype=np.float64)
    averaged = np.einsum("p,pij->ij", fractions[nonzero], values, optimize=True)
    return averaged, dwell


def _as_skycoord(source_coord: Any):
    try:
        import astropy.units as u
        from astropy.coordinates import SkyCoord
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise ImportError("source_coord requires astropy") from exc

    if isinstance(source_coord, SkyCoord):
        return source_coord
    if isinstance(source_coord, (tuple, list)) and len(source_coord) == 2:
        return SkyCoord(l=float(source_coord[0]) * u.deg, b=float(source_coord[1]) * u.deg, frame="galactic")
    raise TypeError("source_coord must be SkyCoord or a (Galactic l_deg, b_deg) tuple")


def _orientation_vectors_and_weights(
    orientation_file: Path,
    source_coord: Any,
    tstart: Any,
    tstop: Any,
    earth_occultation: bool,
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    try:
        import astropy.units as u
        from astropy.time import Time
        from cosipy.spacecraftfile import SpacecraftHistory
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise ImportError(
            "orientation mode requires astropy and the installed COSIpy environment"
        ) from exc

    source = _as_skycoord(source_coord)
    start = tstart if isinstance(tstart, Time) else Time(tstart, scale="utc")
    stop = tstop if isinstance(tstop, Time) else Time(tstop, scale="utc")
    if stop <= start:
        raise ValueError("tstop must be later than tstart")
    history = SpacecraftHistory.open(orientation_file, tstart=start, tstop=stop).select_interval(start, stop)
    local = history.get_target_in_sc_frame(source)
    lon = np.asarray(local.lon.to_value(u.rad), dtype=np.float64)
    lat = np.asarray(local.lat.to_value(u.rad), dtype=np.float64)
    vectors = np.column_stack(
        (np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat))
    )
    times = np.asarray(history.obstime.unix, dtype=np.float64)
    dt = np.diff(times)
    if dt.size == 0 or not np.any(dt > 0):
        raise RuntimeError("orientation interval has no positive time steps")
    if np.any(dt < 0):
        raise RuntimeError("orientation timestamps are not monotonic")
    if not hasattr(history, "get_source_visibility"):
        if earth_occultation:
            raise NotImplementedError("installed COSIpy lacks Earth-occultation visibility")
        interval_livetime = dt
    else:
        interval_livetime = np.asarray(
            history.get_source_visibility(source if earth_occultation else None).to_value(u.s),
            dtype=np.float64,
        )
    if interval_livetime.shape != dt.shape:
        raise RuntimeError(
            f"unexpected COSIpy livetime shape {interval_livetime.shape}; expected {dt.shape}"
        )
    # COSIpy 0.4.2 can append a right-edge sample whose floating-point Unix
    # timestamp is identical to the preceding sample when select_interval()
    # lands exactly on a one-second cadence.  Its corresponding livetime is
    # only roundoff.  Retain all physical intervals and give such duplicate
    # edge samples zero weight instead of rejecting an otherwise valid file.
    interval_livetime = np.where(dt > 0, interval_livetime, 0.0)
    point_exposure = np.zeros(len(times), dtype=np.float64)
    point_exposure[:-1] += 0.5 * interval_livetime
    point_exposure[1:] += 0.5 * interval_livetime
    gal = source.galactic
    icrs = source.icrs
    metadata = {
        "tstart": start.isot,
        "tstop": stop.isot,
        "source_galactic_l_deg": float(gal.l.deg),
        "source_galactic_b_deg": float(gal.b.deg),
        "source_icrs_ra_deg": float(icrs.ra.deg),
        "source_icrs_dec_deg": float(icrs.dec.deg),
        "orientation_samples": int(len(vectors)),
        "intervals": int(len(dt)),
        "visible_livetime_s": float(point_exposure.sum()),
    }
    return vectors, point_exposure, metadata


def _photopeak_effective_area(
    area_matrix_cm2: np.ndarray,
    true_energy_edges_keV: np.ndarray,
    deposit_energy_edges_keV: np.ndarray,
) -> np.ndarray:
    """Full-energy/photopeak Aeff from a deposit matrix.

    For true-energy row ``i``, deposit bins whose centers lie in that true
    energy bin are summed.  This is not a Geant4 interaction-process tag.
    """

    matrix = np.asarray(area_matrix_cm2, dtype=np.float64)
    true_edges = np.asarray(true_energy_edges_keV, dtype=np.float64)
    deposit_edges = np.asarray(deposit_energy_edges_keV, dtype=np.float64)
    if matrix.shape != (len(true_edges) - 1, len(deposit_edges) - 1):
        raise ValueError("area_matrix shape is inconsistent with energy edges")
    centers = 0.5 * (deposit_edges[:-1] + deposit_edges[1:])
    result = np.zeros(len(true_edges) - 1, dtype=np.float64)
    for i, (low, high) in enumerate(zip(true_edges[:-1], true_edges[1:])):
        mask = (centers >= low) & (centers < high)
        if i == len(result) - 1:
            mask |= np.isclose(centers, high, rtol=0.0, atol=1.0e-12)
        result[i] = matrix[i, mask].sum()
    return result


def _deposit_to_channel_probability(
    deposit_edges_keV: np.ndarray,
    measured_edges_keV: np.ndarray,
    calibration: BTOCalibration,
    resolution: BTOResolution,
) -> np.ndarray:
    centers = 0.5 * (deposit_edges_keV[:-1] + deposit_edges_keV[1:])
    sigma = _sigma_keV(centers, resolution)
    low = np.maximum(measured_edges_keV[:-1], calibration.lower_threshold_keV)
    high = measured_edges_keV[1:].copy()
    if calibration.upper_saturation_keV is not None:
        high = np.minimum(high, calibration.upper_saturation_keV)
    valid_channel = high > low
    probability = np.zeros((len(centers), len(low)), dtype=np.float64)
    positive_sigma = sigma > 0
    if np.any(positive_sigma):
        mu = centers[positive_sigma, None]
        width = sigma[positive_sigma, None]
        values = ndtr((high[None, :] - mu) / width) - ndtr((low[None, :] - mu) / width)
        values[:, ~valid_channel] = 0.0
        probability[positive_sigma] = np.clip(values, 0.0, 1.0)
    for index in np.flatnonzero(~positive_sigma):
        energy = centers[index]
        if energy < calibration.lower_threshold_keV:
            continue
        if calibration.upper_saturation_keV is not None and energy >= calibration.upper_saturation_keV:
            continue
        channel = int(np.searchsorted(measured_edges_keV, energy, side="right") - 1)
        if 0 <= channel < len(low) and valid_channel[channel]:
            probability[index, channel] = 1.0
    probability[~np.isfinite(probability)] = 0.0
    return np.clip(probability, 0.0, 1.0)


def _apply_detector_model(
    area_matrix_cm2: np.ndarray,
    deposit_edges_keV: np.ndarray,
    calibration: BTOCalibration,
    resolution: BTOResolution,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    measured_edges = calibration.measured_energy_edges_keV
    probability = _deposit_to_channel_probability(
        deposit_edges_keV, measured_edges, calibration, resolution
    )
    clean = np.nan_to_num(np.asarray(area_matrix_cm2, dtype=np.float64), nan=0.0, posinf=0.0, neginf=0.0)
    if np.any(clean < -1.0e-12):
        raise ValueError("area matrix contains negative values")
    clean = np.clip(clean, 0.0, None)
    rsp = clean @ probability
    rsp = np.clip(np.nan_to_num(rsp, nan=0.0, posinf=0.0, neginf=0.0), 0.0, None)
    arf = rsp.sum(axis=1)
    rmf = np.divide(rsp, arf[:, None], out=np.zeros_like(rsp), where=arf[:, None] > 0)
    return measured_edges, arf, rmf, rsp


def _validate_direction_mode(
    theta_deg: float | None,
    phi_deg: float | None,
    spacecraft_vector: Sequence[float] | None,
    orientation_file: str | Path | None,
    source_coord: Any,
    tstart: Any,
    tstop: Any,
) -> str:
    has_angles = theta_deg is not None or phi_deg is not None
    if has_angles and (theta_deg is None or phi_deg is None):
        raise ValueError("theta_deg and phi_deg must be provided together")
    has_vector = spacecraft_vector is not None
    orientation_parts = (orientation_file, source_coord, tstart, tstop)
    has_orientation = any(item is not None for item in orientation_parts)
    if has_orientation and not all(item is not None for item in orientation_parts):
        missing = [
            name for name, value in zip(
                ("orientation_file", "source_coord", "tstart", "tstop"), orientation_parts
            ) if value is None
        ]
        raise ValueError(f"orientation mode is missing: {', '.join(missing)}")
    modes = int(has_angles) + int(has_vector) + int(has_orientation)
    if modes != 1:
        raise ValueError(
            "specify exactly one direction mode: theta_deg+phi_deg, "
            "spacecraft_vector, or orientation_file+source_coord+tstart+tstop"
        )
    return "orientation" if has_orientation else "fixed"


def _get_bto_response_impl(
    response_h5,
    bto_id,
    *,
    theta_deg=None,
    phi_deg=None,
    spacecraft_vector=None,
    orientation_file=None,
    source_coord=None,
    tstart=None,
    tstop=None,
    calibration=None,
    resolution=None,
    threshold_keV=None,
    output_prefix=None,
    write_ogip=True,
    earth_occultation=False,
    overwrite=False,
):
    """Implement fixed-direction or orientation-averaged response creation.

    Direction modes are mutually exclusive and validated without guessing.
    ``calibration`` and ``resolution`` may be objects or mappings keyed by
    ``BTO1``/``BTO2``.  ``None`` selects explicit placeholder defaults.
    If ``write_ogip`` is true, ``output_prefix`` is required.  Existing OGIP
    products are protected unless ``overwrite=True`` is requested explicitly.
    """

    detector_name = _detector_name(bto_id)
    mode = _validate_direction_mode(
        theta_deg, phi_deg, spacecraft_vector, orientation_file, source_coord, tstart, tstop
    )
    calibration = DEFAULT_CALIBRATION if calibration is None else calibration
    resolution = Prototype_Resolution if resolution is None else resolution
    calibration = _select_for_detector(calibration, detector_name, BTOCalibration)
    resolution = _select_for_detector(resolution, detector_name, BTOResolution)
    if threshold_keV is not None:
        calibration = replace(calibration, lower_threshold_keV=float(threshold_keV))
    if write_ogip and output_prefix is None:
        raise ValueError("output_prefix is required when write_ogip=True")

    response_path = Path(response_h5).expanduser().resolve()
    if not response_path.exists():
        raise FileNotFoundError(response_path)
    metadata: dict[str, Any] = {
        "module_version": __version__,
        "created_utc": _utc_now(),
        "detector_name": detector_name,
        "direction_mode": mode,
        "unsmeared_components_available": True,
        "earth_occultation": bool(earth_occultation),
        "input_response_h5": str(response_path),
        "input_response_h5_sha256": _sha256(response_path),
        "angular_interpolation": "4-nearest inverse-angular-distance-squared",
        "gain_convention": (
            "channel = gain_quadratic_ch_per_keV2 * E_keV**2 + "
            "gain_slope_ch_per_keV * E_keV + gain_offset_ch"
        ),
        "calibration": asdict(calibration),
        "resolution": asdict(resolution),
        "calibration_status": "PLACEHOLDER; not flight calibration",
        "software_versions": _dependency_versions(),
    }

    with h5py.File(response_path, "r") as handle:
        required = (
            "detector/name", "direction/vector", "energy/true_edges_keV",
            "energy/deposit_edges_keV", "response/area_matrix_cm2",
        )
        missing = [name for name in required if name not in handle]
        if missing:
            raise KeyError(f"response HDF5 is missing datasets: {missing}")
        names = [item.decode() if isinstance(item, bytes) else str(item) for item in handle["detector/name"][:]]
        if detector_name not in names:
            raise KeyError(f"{detector_name} not present in detector/name={names}")
        detector_index = names.index(detector_name)
        response_vectors = np.asarray(handle["direction/vector"][:], dtype=np.float64)
        complete = np.asarray(
            handle["simulation/complete"][:] if "simulation/complete" in handle else np.ones(len(response_vectors)),
            dtype=bool,
        )
        detector_matrix = handle["response/area_matrix_cm2"][detector_index]
        true_edges = np.asarray(handle["energy/true_edges_keV"][:], dtype=np.float64)
        deposit_edges = np.asarray(handle["energy/deposit_edges_keV"][:], dtype=np.float64)
        for key in ("format", "schema_version", "config_sha256", "geometry_sha256", "healpix_nside", "direction_scheme", "normalization"):
            if key in handle.attrs:
                metadata[f"hdf5_{key}"] = _jsonable(handle.attrs[key])

        if mode == "fixed":
            if spacecraft_vector is not None:
                vector = _normalize_vector(spacecraft_vector)
                input_form = "spacecraft_vector"
            else:
                if not (0.0 <= float(theta_deg) <= 180.0):
                    raise ValueError("theta_deg must be in [0, 180]")
                vector = _vector_from_angles(float(theta_deg), float(phi_deg))
                input_form = "theta_phi"
            area_matrix, selected, weights = _interpolate_area_matrix(
                detector_matrix, response_vectors, vector, complete=complete, neighbors=4
            )
            query_theta, query_phi = _angles_from_vector(vector)
            metadata.update(
                {
                    "fixed_direction_input": input_form,
                    "spacecraft_vector": vector.tolist(),
                    "theta_deg": query_theta,
                    "phi_deg": query_phi,
                    "neighbor_direction_indices": selected.tolist(),
                    "neighbor_weights": weights.tolist(),
                    "orientation_samples": 0,
                }
            )
        else:
            orientation_path = Path(orientation_file).expanduser().resolve()
            if not orientation_path.exists():
                raise FileNotFoundError(orientation_path)
            vectors, exposures, orientation_metadata = _orientation_vectors_and_weights(
                orientation_path, source_coord, tstart, tstop, bool(earth_occultation)
            )
            area_matrix, dwell = _average_orientation_samples(
                detector_matrix, response_vectors, vectors, exposures, complete=complete, neighbors=4
            )
            metadata.update(orientation_metadata)
            metadata.update(
                {
                    "orientation_file": str(orientation_path),
                    "orientation_file_sha256": _sha256(orientation_path),
                    "nonzero_response_directions": int(np.count_nonzero(dwell)),
                    "direction_dwell_seconds": dwell.tolist(),
                }
            )

    area_matrix = np.asarray(area_matrix, dtype=np.float64)
    unsmeared_arf = area_matrix.sum(axis=1)
    photopeak_arf = _photopeak_effective_area(area_matrix, true_edges, deposit_edges)
    measured_edges, arf, rmf, rsp = _apply_detector_model(
        area_matrix, deposit_edges, calibration, resolution
    )
    product = BTOResponseProduct(
        detector_name=detector_name,
        photon_energy_edges_keV=true_edges,
        channel_edges=calibration.channel_edges,
        measured_energy_edges_keV=measured_edges,
        arf_cm2=arf,
        rmf_probability=rmf,
        rsp_cm2=rsp,
        unsmeared_area_matrix_cm2=area_matrix,
        unsmeared_arf_cm2=unsmeared_arf,
        deposit_energy_edges_keV=deposit_edges,
        photopeak_arf_cm2=photopeak_arf,
        calibration=calibration,
        resolution=resolution,
        metadata=metadata,
    )
    if write_ogip:
        write_ogip_response(product, output_prefix, overwrite=bool(overwrite))
    return product


def _fits_module():
    try:
        from astropy.io import fits
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise ImportError("OGIP FITS I/O requires astropy") from exc
    return fits


def _metadata_history(header, metadata: Mapping[str, Any]) -> None:
    keys = (
        "module_version", "created_utc", "detector_name", "direction_mode",
        "theta_deg", "phi_deg", "orientation_file", "tstart", "tstop",
        "earth_occultation", "input_response_h5", "input_response_h5_sha256",
        "hdf5_geometry_sha256", "hdf5_config_sha256", "angular_interpolation",
        "gain_convention", "orientation_samples",
    )
    for key in keys:
        if key in metadata:
            header.add_history(f"{key}={metadata[key]}")
    if "software_versions" in metadata:
        header.add_history(
            "software_versions="
            + json.dumps(metadata["software_versions"], sort_keys=True, separators=(",", ":"))
        )


def _base_response_header(header, product: BTOResponseProduct) -> None:
    header["HDUCLASS"] = "OGIP"
    header["HDUCLAS1"] = "RESPONSE"
    header["TELESCOP"] = "COSI"
    header["INSTRUME"] = "BTO"
    header["DETNAM"] = product.detector_name
    header["FILTER"] = "NONE"
    header["CHANTYPE"] = "PI"
    header["DETCHANS"] = len(product.channel_edges) - 1
    header["GAINA"] = product.calibration.gain_quadratic_ch_per_keV2
    header["GAINB"] = product.calibration.gain_slope_ch_per_keV
    header["GAINC"] = product.calibration.gain_offset_ch
    # Retain the legacy linear keywords so older readers still understand
    # products whose quadratic coefficient is zero.
    header["GAINSLP"] = product.calibration.gain_slope_ch_per_keV
    header["GAINOFF"] = product.calibration.gain_offset_ch
    header["LOTHRESH"] = product.calibration.lower_threshold_keV
    if product.calibration.upper_saturation_keV is not None:
        header["HISAT"] = product.calibration.upper_saturation_keV
    header["RESOLMOD"] = product.resolution.name[:68]
    header["RESA0"] = product.resolution.a0_keV2
    header["RESA1"] = product.resolution.a1_keV
    header["RESA2"] = product.resolution.a2
    header.add_history("Calibration preset is a configurable placeholder, not flight calibration.")
    header.add_history(f"Calibration provenance: {product.calibration.provenance}")
    header.add_history(f"Resolution provenance: {product.resolution.provenance}")
    _metadata_history(header, product.metadata)


def _ebounds_hdu(product: BTOResponseProduct):
    fits = _fits_module()
    n_channel = len(product.channel_edges) - 1
    columns = [
        fits.Column(name="CHANNEL", format="J", array=np.arange(n_channel, dtype=np.int32)),
        fits.Column(name="E_MIN", format="D", unit="keV", array=product.measured_energy_edges_keV[:-1]),
        fits.Column(name="E_MAX", format="D", unit="keV", array=product.measured_energy_edges_keV[1:]),
    ]
    hdu = fits.BinTableHDU.from_columns(columns, name="EBOUNDS")
    _base_response_header(hdu.header, product)
    hdu.header["HDUCLAS2"] = "EBOUNDS"
    hdu.header["HDUVERS"] = "1.1.0"
    hdu.header["TLMIN1"] = 0
    hdu.header["TLMAX1"] = n_channel - 1
    return hdu


def _matrix_hdu(product: BTOResponseProduct, matrix: np.ndarray, *, combined: bool):
    fits = _fits_module()
    n_true, n_channel = matrix.shape
    columns = [
        fits.Column(name="ENERG_LO", format="D", unit="keV", array=product.photon_energy_edges_keV[:-1]),
        fits.Column(name="ENERG_HI", format="D", unit="keV", array=product.photon_energy_edges_keV[1:]),
        fits.Column(name="N_GRP", format="I", array=np.ones(n_true, dtype=np.int16)),
        fits.Column(name="F_CHAN", format="J", array=np.zeros(n_true, dtype=np.int32)),
        fits.Column(name="N_CHAN", format="J", array=np.full(n_true, n_channel, dtype=np.int32)),
        fits.Column(
            name="MATRIX", format=f"{n_channel}E", unit="cm**2" if combined else None,
            array=np.asarray(matrix, dtype=np.float32),
        ),
    ]
    hdu = fits.BinTableHDU.from_columns(columns, name="MATRIX")
    _base_response_header(hdu.header, product)
    hdu.header["HDUCLAS2"] = "RSP_MATRIX"
    hdu.header["HDUCLAS3"] = "FULL" if combined else "REDIST"
    hdu.header["HDUVERS"] = "1.3.0"
    hdu.header["TLMIN4"] = 0
    hdu.header["TLMAX4"] = n_channel - 1
    hdu.header["LO_THRES"] = 0.0
    hdu.header.add_history(
        "MATRIX is RSP [cm2]" if combined else "MATRIX is conditional RMF probability"
    )
    return hdu


def _arf_hdu(product: BTOResponseProduct):
    fits = _fits_module()
    columns = [
        fits.Column(name="ENERG_LO", format="E", unit="keV", array=product.photon_energy_edges_keV[:-1].astype(np.float32)),
        fits.Column(name="ENERG_HI", format="E", unit="keV", array=product.photon_energy_edges_keV[1:].astype(np.float32)),
        fits.Column(name="SPECRESP", format="E", unit="cm**2", array=product.arf_cm2.astype(np.float32)),
    ]
    hdu = fits.BinTableHDU.from_columns(columns, name="SPECRESP")
    _base_response_header(hdu.header, product)
    hdu.header["HDUCLAS2"] = "SPECRESP"
    hdu.header["HDUVERS"] = "1.1.0"
    hdu.header.add_history("SPECRESP includes lower-threshold, saturation, and channel-range losses.")
    return hdu


def _write_ogip_response_impl(
    product: BTOResponseProduct,
    output_prefix,
    *,
    overwrite: bool = False,
) -> BTOResponseProduct:
    """Write standalone ARF, conditional RMF, combined RSP, and JSON metadata."""

    fits = _fits_module()
    prefix = Path(output_prefix).expanduser().resolve()
    prefix.parent.mkdir(parents=True, exist_ok=True)
    arf_path = Path(f"{prefix}.arf")
    rmf_path = Path(f"{prefix}.rmf")
    rsp_path = Path(f"{prefix}.rsp")
    metadata_path = Path(f"{prefix}.response.json")
    for path in (arf_path, rmf_path, rsp_path, metadata_path):
        if path.exists() and not overwrite:
            raise FileExistsError(f"output exists; pass overwrite=True: {path}")

    primary = fits.PrimaryHDU()
    primary.header["TELESCOP"] = "COSI"
    primary.header["INSTRUME"] = "BTO"
    primary.header["DETNAM"] = product.detector_name
    fits.HDUList([primary.copy(), _arf_hdu(product)]).writeto(arf_path, overwrite=overwrite, checksum=True)
    rmf_matrix = _matrix_hdu(product, product.rmf_probability, combined=False)
    rmf_matrix.header["ANCRFILE"] = arf_path.name
    fits.HDUList([primary.copy(), rmf_matrix, _ebounds_hdu(product)]).writeto(
        rmf_path, overwrite=overwrite, checksum=True
    )
    fits.HDUList([primary.copy(), _matrix_hdu(product, product.rsp_cm2, combined=True), _ebounds_hdu(product)]).writeto(
        rsp_path, overwrite=overwrite, checksum=True
    )
    sidecar = {
        "format": "COSI-BTO-OGIP-response-metadata",
        "module_version": __version__,
        "detector_name": product.detector_name,
        "arf_path": str(arf_path),
        "rmf_path": str(rmf_path),
        "rsp_path": str(rsp_path),
        "rsp_sha256": _sha256(rsp_path),
        "rmf_sha256": _sha256(rmf_path),
        "calibration": asdict(product.calibration),
        "resolution": asdict(product.resolution),
        "metadata": product.metadata,
    }
    metadata_path.write_text(json.dumps(_jsonable(sidecar), indent=2, sort_keys=True) + "\n")
    product.arf_path = arf_path
    product.rmf_path = rmf_path
    product.rsp_path = rsp_path
    product.metadata_path = metadata_path
    return product


def _read_dense_ogip_matrix(hdu) -> np.ndarray:
    data = hdu.data
    n_channel = int(hdu.header["DETCHANS"])
    first = int(hdu.header.get(f"TLMIN{data.columns.names.index('F_CHAN') + 1}", 1))
    output = np.zeros((len(data), n_channel), dtype=np.float64)
    for row_index in range(len(data)):
        n_group = int(np.squeeze(data["N_GRP"][row_index]))
        fchan = np.atleast_1d(data["F_CHAN"][row_index]).astype(int) - first
        nchan = np.atleast_1d(data["N_CHAN"][row_index]).astype(int)
        values = np.atleast_1d(data["MATRIX"][row_index]).astype(float)
        offset = 0
        for group in range(n_group):
            count = int(nchan[group])
            start = int(fchan[group])
            output[row_index, start : start + count] = values[offset : offset + count]
            offset += count
    return output


def _load_bto_response_impl(response_file, *, arf_file=None) -> BTOResponseProduct:
    """Implement loading of a module-written RSP or RMF+ARF pair."""

    fits = _fits_module()
    response_path = Path(response_file).expanduser().resolve()
    with fits.open(response_path, memmap=False, checksum=True) as hdul:
        matrix_hdu = hdul["MATRIX"]
        matrix = _read_dense_ogip_matrix(matrix_hdu)
        true_lo = matrix_hdu.data["ENERG_LO"].astype(float)
        true_hi = matrix_hdu.data["ENERG_HI"].astype(float)
        true_edges = np.append(true_lo, true_hi[-1])
        ebounds = hdul["EBOUNDS"].data
        measured_edges = np.append(ebounds["E_MIN"].astype(float), float(ebounds["E_MAX"][-1]))
        detector_name = matrix_hdu.header.get("DETNAM", "BTO1")
        combined = str(matrix_hdu.header.get("HDUCLAS3", "")).strip().upper() == "FULL"
        cal = BTOCalibration(
            n_channels=len(measured_edges) - 1,
            gain_quadratic_ch_per_keV2=float(matrix_hdu.header.get("GAINA", 0.0)),
            gain_slope_ch_per_keV=float(
                matrix_hdu.header.get(
                    "GAINB",
                    matrix_hdu.header.get("GAINSLP", 1.0 / np.median(np.diff(measured_edges))),
                )
            ),
            gain_offset_ch=float(
                matrix_hdu.header.get("GAINC", matrix_hdu.header.get("GAINOFF", 0.0))
            ),
            lower_threshold_keV=float(matrix_hdu.header.get("LOTHRESH", measured_edges[0])),
            upper_saturation_keV=matrix_hdu.header.get("HISAT", None),
            name="Loaded from OGIP header",
            provenance=str(response_path),
        )
        res = BTOResolution(
            float(matrix_hdu.header.get("RESA0", 0.0)),
            float(matrix_hdu.header.get("RESA1", 0.0)),
            float(matrix_hdu.header.get("RESA2", 0.0)),
            str(matrix_hdu.header.get("RESOLMOD", "Loaded from OGIP header")),
            str(response_path),
        )
    if combined:
        rsp = matrix
        arf = rsp.sum(axis=1)
        rmf = np.divide(rsp, arf[:, None], out=np.zeros_like(rsp), where=arf[:, None] > 0)
    else:
        rmf = matrix
        if arf_file is None:
            raise ValueError("arf_file is required when loading a standalone RMF")
        with fits.open(Path(arf_file).expanduser().resolve(), memmap=False, checksum=True) as hdul:
            arf_data = hdul["SPECRESP"].data
            # Check every lower AND upper edge before multiplying rows by area.
            if len(arf_data) != len(true_lo) or any(
                not np.allclose(arf_data[column], edges, rtol=1e-7, atol=0.0)
                for column, edges in (("ENERG_LO", true_lo), ("ENERG_HI", true_hi))
            ):
                raise ValueError("ARF true-energy bins (ENERG_LO/ENERG_HI) do not match RMF")
            arf = np.asarray(arf_data["SPECRESP"], dtype=float)
        rsp = arf[:, None] * rmf
    empty_deposit_edges = true_edges.copy()
    empty_area = np.zeros((len(true_edges) - 1, len(true_edges) - 1), dtype=float)
    # A matching module-written sidecar preserves the source/interval provenance.
    # The actual numerical response and detector model always come from FITS.
    metadata = {"direction_mode": "loaded_ogip"}
    metadata_path = response_path.with_suffix('.response.json')
    if metadata_path.exists():
        sidecar = json.loads(metadata_path.read_text())
        if sidecar.get('detector_name') != detector_name:
            raise ValueError("response metadata detector does not match FITS DETNAM")
        digest_key = 'rsp_sha256' if combined else 'rmf_sha256'
        if digest_key in sidecar and sidecar[digest_key] != _sha256(response_path):
            raise ValueError("response metadata checksum does not match FITS response")
        metadata.update(sidecar.get('metadata', {}))
    metadata.update(loaded_from=str(response_path), unsmeared_components_available=False)
    return BTOResponseProduct(
        detector_name=detector_name,
        photon_energy_edges_keV=true_edges,
        channel_edges=np.arange(len(measured_edges), dtype=float),
        measured_energy_edges_keV=measured_edges,
        arf_cm2=arf,
        rmf_probability=rmf,
        rsp_cm2=rsp,
        unsmeared_area_matrix_cm2=empty_area,
        unsmeared_arf_cm2=np.zeros(len(true_edges) - 1),
        deposit_energy_edges_keV=empty_deposit_edges,
        photopeak_arf_cm2=np.zeros(len(true_edges) - 1),
        calibration=cal,
        resolution=res,
        metadata=metadata,
        arf_path=Path(arf_file).resolve() if arf_file is not None else None,
        rmf_path=response_path if not combined else None,
        rsp_path=response_path if combined else None,
        metadata_path=metadata_path if metadata_path.exists() else None,
    )


def _integrate_callable(model: Callable, edges: np.ndarray, order: int = 16) -> np.ndarray:
    nodes, weights = np.polynomial.legendre.leggauss(order)
    midpoint = 0.5 * (edges[:-1] + edges[1:])
    halfwidth = 0.5 * (edges[1:] - edges[:-1])
    energies = midpoint[:, None] + halfwidth[:, None] * nodes[None, :]
    # Evaluate on a flat vector.  In particular, astromodels 2.5's numba Band
    # implementation accepts a one-dimensional energy array but fails on the
    # natural [true-bin, quadrature-node] grid used by this integrator.
    evaluation_energies = energies.reshape(-1)
    try:
        values = model(evaluation_energies)
    except Exception as first_error:
        try:
            import astropy.units as u
            values = model(evaluation_energies * u.keV)
        except Exception:
            raise TypeError(
                "photon_model callable failed for numeric keV and Astropy keV inputs"
            ) from first_error
    if hasattr(values, "to_value"):
        try:
            import astropy.units as u
            values = values.to_value(1 / (u.cm**2 * u.s * u.keV))
        except Exception:
            values = values.value
    values = np.asarray(values, dtype=np.float64)
    try:
        values = np.broadcast_to(values, evaluation_energies.shape).reshape(energies.shape)
    except ValueError as exc:
        raise ValueError(
            f"photon_model returned shape {values.shape}; expected a scalar or "
            f"{evaluation_energies.shape} values"
        ) from exc
    if np.any(~np.isfinite(values)) or np.any(values < 0):
        raise ValueError("photon_model must return finite nonnegative dN/dE")
    return halfwidth * np.sum(weights[None, :] * values, axis=1)


def _integrated_flux(photon_model: Any, edges: np.ndarray) -> np.ndarray:
    if callable(photon_model):
        return _integrate_callable(photon_model, edges)
    flux = np.asarray(photon_model, dtype=np.float64)
    if flux.shape != (len(edges) - 1,):
        raise ValueError("integrated photon flux array must have one value per true-energy bin")
    if np.any(~np.isfinite(flux)) or np.any(flux < 0):
        raise ValueError("integrated photon flux must be finite and nonnegative")
    return flux


def _fold_photon_spectrum(response: BTOResponseProduct, photon_model: Any, exposure_s: float = 1.0) -> tuple[np.ndarray, np.ndarray]:
    """Return integrated true-bin flux and expected source counts/channel."""

    if not np.isfinite(exposure_s) or exposure_s < 0:
        raise ValueError("exposure_s must be finite and nonnegative")
    flux = _integrated_flux(photon_model, response.photon_energy_edges_keV)
    counts = float(exposure_s) * np.einsum("e,ec->c", flux, response.rsp_cm2, optimize=True)
    return flux, np.clip(counts, 0.0, None)


def _integrate_piecewise_histogram(counts: np.ndarray, source_edges: np.ndarray, target_edges: np.ndarray) -> np.ndarray:
    counts = np.asarray(counts, dtype=np.float64)
    source_edges = np.asarray(source_edges, dtype=np.float64)
    target_edges = np.asarray(target_edges, dtype=np.float64)
    density = np.divide(counts, np.diff(source_edges), out=np.zeros_like(counts), where=np.diff(source_edges) > 0)
    cumulative = np.concatenate(([0.0], np.cumsum(counts)))

    def cumulative_at(energy):
        energy = np.asarray(energy, dtype=float)
        clipped = np.clip(energy, source_edges[0], source_edges[-1])
        index = np.searchsorted(source_edges, clipped, side="right") - 1
        index = np.clip(index, 0, len(counts) - 1)
        value = cumulative[index] + density[index] * (clipped - source_edges[index])
        value = np.where(energy <= source_edges[0], 0.0, value)
        return np.where(energy >= source_edges[-1], cumulative[-1], value)

    return np.clip(cumulative_at(target_edges[1:]) - cumulative_at(target_edges[:-1]), 0.0, None)


def _get_bto_background_impl(
    response: BTOResponseProduct,
    *,
    rate_per_channel: float | Sequence[float] | None = None,
    total_rate_hz: float | None = None,
    background_template: str | Path | None = None,
    background_pha: str | Path | None = None,
    background_rmf: str | Path | None = None,
    background_energy_edges_keV: Sequence[float] | None = None,
) -> BTOBackground:
    """Return a time-independent BTO background in counts/s/channel.

    Exactly one source can be supplied.  With no explicit source, the existing
    Takashima single-BTO mean template and its embedded EBOUNDS are used. This is constant
    in time but retains the simulated energy spectrum.  ``rate_per_channel`` is
    counts/s/channel; a scalar is repeated.  ``total_rate_hz`` is spread
    uniformly over channels whose energy overlaps the response acquisition
    threshold/saturation range.
    """

    explicit = sum(x is not None for x in (rate_per_channel, total_rate_hz, background_pha, background_template))
    if explicit > 1:
        raise ValueError("choose only one of rate_per_channel, total_rate_hz, background_template or background_pha")
    if background_rmf is not None and background_energy_edges_keV is not None:
        raise ValueError("choose background_energy_edges_keV or background_rmf, not both")
    if (rate_per_channel is not None or total_rate_hz is not None) and (
        background_rmf is not None or background_energy_edges_keV is not None
    ):
        raise ValueError("background energy-grid options apply only to a PHA background")
    target_edges = response.measured_energy_edges_keV
    n_channel = len(target_edges) - 1
    if rate_per_channel is not None:
        rate = np.asarray(rate_per_channel, dtype=float)
        if rate.ndim == 0:
            rate = np.full(n_channel, float(rate))
        if rate.shape != (n_channel,):
            raise ValueError("rate_per_channel must be scalar or have n_channel values")
        provenance = "user counts/s/channel"
        path = None
    elif total_rate_hz is not None:
        total = float(total_rate_hz)
        if not np.isfinite(total) or total < 0:
            raise ValueError("total_rate_hz must be finite and nonnegative")
        active = target_edges[1:] > response.calibration.lower_threshold_keV
        if response.calibration.upper_saturation_keV is not None:
            active &= target_edges[:-1] < response.calibration.upper_saturation_keV
        if not np.any(active):
            raise ValueError("response has no channels inside threshold/saturation")
        rate = np.zeros(n_channel)
        rate[active] = total / np.count_nonzero(active)
        provenance = "user total rate uniformly distributed across active channels"
        path = None
    else:
        fits = _fits_module()
        pha_path = Path(background_template or background_pha or DEFAULT_REFERENCE_BACKGROUND).expanduser().resolve()
        if not pha_path.exists():
            raise FileNotFoundError(
                f"reference background not found: {pha_path}; pass rate_per_channel or total_rate_hz"
            )
        embedded_edges = None
        template_origin = "user background spectrum"
        with fits.open(pha_path, memmap=False) as hdul:
            template_origin = hdul[0].header.get('ORIGIN', template_origin)
            data = hdul["SPECTRUM"].data
            names = set(data.names)
            exposure = float(hdul["SPECTRUM"].header.get("EXPOSURE", 1.0))
            channels = np.asarray(data["CHANNEL"], dtype=int)
            if "RATE" in names:
                native_rate = np.asarray(data["RATE"], dtype=float)
            elif "COUNTS" in names:
                if not np.isfinite(exposure) or exposure <= 0:
                    raise ValueError("COUNTS background requires finite positive EXPOSURE")
                native_rate = np.asarray(data["COUNTS"], dtype=float) / exposure
            else:
                raise ValueError("background PHA contains neither RATE nor COUNTS")
            if 'EBOUNDS' in hdul:
                ebounds = hdul['EBOUNDS'].data
                if not np.array_equal(ebounds['CHANNEL'], channels):
                    raise ValueError('template EBOUNDS and SPECTRUM CHANNEL values must agree')
                import astropy.units as u
                lo = np.asarray(ebounds['E_MIN'], float) * u.Unit(hdul['EBOUNDS'].columns['E_MIN'].unit or 'keV').to(u.keV)
                hi = np.asarray(ebounds['E_MAX'], float) * u.Unit(hdul['EBOUNDS'].columns['E_MAX'].unit or 'keV').to(u.keV)
                if np.any(hi <= lo) or not np.allclose(hi[:-1], lo[1:], rtol=0, atol=1e-8):
                    raise ValueError('template EBOUNDS must have positive, contiguous bins')
                embedded_edges = np.append(lo, hi[-1])
        if np.any(~np.isfinite(native_rate)) or np.any(native_rate < 0):
            raise ValueError("background input rate must be finite and nonnegative")
        if background_energy_edges_keV is not None:
            source_edges = np.asarray(background_energy_edges_keV, dtype=float)
            grid_provenance = "user-supplied background energy boundaries [keV]"
        elif background_rmf is not None:
            rmf_path = Path(background_rmf).expanduser().resolve()
            with fits.open(rmf_path, memmap=False) as hdul:
                ebounds = hdul["EBOUNDS"].data
                if not np.array_equal(ebounds["CHANNEL"], channels):
                    raise ValueError("PHA and background RMF CHANNEL values must agree")
                if not np.allclose(ebounds["E_MAX"][:-1], ebounds["E_MIN"][1:]):
                    raise ValueError("background RMF energy bins must be contiguous")
                source_edges = np.append(np.asarray(ebounds["E_MIN"], float), float(ebounds["E_MAX"][-1]))
            grid_provenance = "explicit background RMF EBOUNDS"
        elif embedded_edges is not None:
            source_edges = embedded_edges
            grid_provenance = "embedded EBOUNDS in background template"
        else:
            raise ValueError("background spectrum requires embedded EBOUNDS, background_energy_edges_keV or background_rmf")
        if (source_edges.shape != (len(native_rate) + 1,)
                or np.any(~np.isfinite(source_edges)) or np.any(np.diff(source_edges) <= 0)):
            raise ValueError("background energy edges must be finite, increasing and have n_channel+1 values")
        rate = _integrate_piecewise_histogram(native_rate, source_edges, target_edges)
        provenance = f"{template_origin}; constant in time"
        # The source BAK has the old response's channelization.  Keep it in
        # metadata, but do not advertise it as BACKFILE until the rebinned rate
        # has been written with write_ogip_background().
        path = None
    if np.any(~np.isfinite(rate)) or np.any(rate < 0):
        raise ValueError("background rate contains non-finite or negative values")
    return BTOBackground(
        rate_per_channel=rate,
        measured_energy_edges_keV=target_edges.copy(),
        metadata={
            "unit": "count/s/channel",
            "provenance": provenance,
            "total_rate_hz": float(rate.sum()),
            "time_dependence": "constant",
            "source_background_pha": str(pha_path) if 'pha_path' in locals() else None,
            "source_background_template": str(pha_path) if 'pha_path' in locals() else None,
            "source_background_sha256": _sha256(pha_path) if 'pha_path' in locals() else None,
            "source_background_rmf": str(rmf_path) if 'rmf_path' in locals() else None,
            "background_energy_grid": grid_provenance if 'grid_provenance' in locals() else "response channel grid",
            "source_energy_edges_keV": source_edges.tolist() if 'source_edges' in locals() else None,
        },
        path=path,
    )


def _write_ogip_background_impl(
    background: BTOBackground,
    response: BTOResponseProduct,
    output_bak,
    *,
    exposure_s: float = 40.0,
    overwrite: bool = False,
) -> BTOBackground:
    """Write a response-channelized, time-independent OGIP background BAK.

    ``RATE`` and ``STAT_ERR`` are counts/s/channel.  ``STAT_ERR`` is the
    Poisson-equivalent uncertainty ``sqrt(rate*exposure)/exposure``; the input
    mean itself is not randomized.
    """

    fits = _fits_module()
    exposure_s = float(exposure_s)
    if not np.isfinite(exposure_s) or exposure_s <= 0:
        raise ValueError("exposure_s must be finite and positive")
    rate = np.asarray(background.rate_per_channel, dtype=np.float64)
    n_channel = response.rsp_cm2.shape[1]
    if rate.shape != (n_channel,):
        raise ValueError("background channel count does not match response")
    if not np.allclose(
        background.measured_energy_edges_keV,
        response.measured_energy_edges_keV,
        rtol=0.0,
        atol=1.0e-10,
    ):
        raise ValueError("background measured-energy edges do not match response")
    rate_error = np.sqrt(np.clip(rate * exposure_s, 0.0, None)) / exposure_s
    columns = [
        fits.Column(name="CHANNEL", format="J", array=np.arange(n_channel, dtype=np.int32)),
        fits.Column(name="RATE", format="D", unit="count/s", array=rate),
        fits.Column(name="STAT_ERR", format="D", unit="count/s", array=rate_error),
        fits.Column(name="SYS_ERR", format="D", array=np.zeros(n_channel)),
        fits.Column(name="QUALITY", format="I", array=np.zeros(n_channel, dtype=np.int16)),
        fits.Column(name="GROUPING", format="I", array=np.ones(n_channel, dtype=np.int16)),
        fits.Column(name="AREASCAL", format="E", array=np.ones(n_channel, dtype=np.float32)),
        fits.Column(name="BACKSCAL", format="E", array=np.ones(n_channel, dtype=np.float32)),
    ]
    spectrum = fits.BinTableHDU.from_columns(columns, name="SPECTRUM")
    header = spectrum.header
    header["HDUCLASS"] = "OGIP"
    header["HDUCLAS1"] = "SPECTRUM"
    header["HDUCLAS2"] = "BKG"
    header["HDUCLAS3"] = "RATE"
    header["HDUVERS"] = "1.2.1"
    header["TELESCOP"] = "COSI"
    header["INSTRUME"] = "BTO"
    header["DETNAM"] = response.detector_name
    header["FILTER"] = "NONE"
    header["CHANTYPE"] = "PI"
    header["DETCHANS"] = n_channel
    header["TLMIN1"] = 0
    header["TLMAX1"] = n_channel - 1
    header["EXPOSURE"] = exposure_s
    header["POISSERR"] = False
    header["AREASCAL"] = 1.0
    header["BACKSCAL"] = 1.0
    header["RESPFILE"] = response.rmf_path.name if response.rmf_path else (response.rsp_path.name if response.rsp_path else "NONE")
    header["ANCRFILE"] = response.arf_path.name if response.arf_path else "NONE"
    header["BACKFILE"] = "NONE"
    header["CORRFILE"] = "NONE"
    header["CORRSCAL"] = 1.0
    header.add_history("Time-independent mean BTO background; RATE was not Poisson randomized.")
    header.add_history("STAT_ERR=sqrt(RATE*EXPOSURE)/EXPOSURE.")
    header.add_history(f"Background provenance: {background.metadata.get('provenance', 'unspecified')}")
    _metadata_history(header, response.metadata)
    primary = fits.PrimaryHDU()
    primary.header["TELESCOP"] = "COSI"
    primary.header["INSTRUME"] = "BTO"
    primary.header["DETNAM"] = response.detector_name
    path = Path(output_bak).expanduser().resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    fits.HDUList([primary, spectrum]).writeto(path, overwrite=overwrite, checksum=True)
    background.path = path
    background.metadata["written_background_bak"] = str(path)
    background.metadata["statistical_exposure_s"] = exposure_s
    return background


def _background_counts(background: Any, response: BTOResponseProduct, exposure_s: float) -> tuple[np.ndarray, str]:
    n_channel = response.rsp_cm2.shape[1]
    if background is None:
        return np.zeros(n_channel), "NONE"
    if isinstance(background, BTOBackground):
        if background.rate_per_channel.shape != (n_channel,):
            raise ValueError("BTOBackground channel count does not match response")
        return background.rate_per_channel * exposure_s, background.path.name if background.path else "NONE"
    if isinstance(background, Mapping):
        allowed = set(background) & {"rate_per_channel", "counts_per_channel"}
        if len(allowed) != 1:
            raise ValueError("background mapping needs exactly one of rate_per_channel/counts_per_channel")
        key = next(iter(allowed))
        values = np.asarray(background[key], dtype=float)
        counts = values * exposure_s if key == "rate_per_channel" else values
        source = "NONE"
    elif isinstance(background, (str, Path)):
        loaded = get_bto_background(response, background_pha=background)
        return loaded.rate_per_channel * exposure_s, Path(background).name
    else:
        values = np.asarray(background, dtype=float)
        counts = values * exposure_s  # bare arrays are explicitly rates
        source = "NONE"
    if counts.shape != (n_channel,):
        raise ValueError("background array must have n_channel values")
    if np.any(~np.isfinite(counts)) or np.any(counts < 0):
        raise ValueError("background counts must be finite and nonnegative")
    return counts, source


def _write_pha(
    output_path: Path,
    counts: np.ndarray,
    exposure_s: float,
    response: BTOResponseProduct,
    background_file: str,
    *,
    poisson: bool,
    overwrite: bool,
) -> None:
    fits = _fits_module()
    n_channel = len(counts)
    columns = [
        fits.Column(name="CHANNEL", format="J", array=np.arange(n_channel, dtype=np.int32)),
        (fits.Column(name="COUNTS", format="J", unit="count", array=np.asarray(counts, dtype=np.int64))
         if poisson else fits.Column(name="RATE", format="D", unit="count/s", array=counts / exposure_s)),
        fits.Column(name="QUALITY", format="I", array=np.zeros(n_channel, dtype=np.int16)),
        fits.Column(name="GROUPING", format="I", array=np.ones(n_channel, dtype=np.int16)),
        fits.Column(name="AREASCAL", format="E", array=np.ones(n_channel, dtype=np.float32)),
        fits.Column(name="BACKSCAL", format="E", array=np.ones(n_channel, dtype=np.float32)),
    ]
    if not poisson:
        # Preserve fractional expectations; specify counting errors explicitly.
        columns.append(fits.Column(name="STAT_ERR", format="D", unit="count/s",
                                   array=np.sqrt(counts) / exposure_s))
    spectrum = fits.BinTableHDU.from_columns(columns, name="SPECTRUM")
    header = spectrum.header
    header["HDUCLASS"] = "OGIP"
    header["HDUCLAS1"] = "SPECTRUM"
    header["HDUCLAS2"] = "TOTAL"
    header["HDUCLAS3"] = "COUNT" if poisson else "RATE"
    header["HDUVERS"] = "1.2.1"
    header["TELESCOP"] = "COSI"
    header["INSTRUME"] = "BTO"
    header["DETNAM"] = response.detector_name
    header["FILTER"] = "NONE"
    header["CHANTYPE"] = "PI"
    header["DETCHANS"] = n_channel
    header["TLMIN1"] = 0
    header["TLMAX1"] = n_channel - 1
    header["EXPOSURE"] = float(exposure_s)
    header["POISSERR"] = bool(poisson)
    if not poisson:
        header.add_history("Unsampled expectation; STAT_ERR=sqrt(expected counts)/EXPOSURE.")
        header.add_history("STAT_ERR is a Poisson-equivalent assumption, not measured scatter.")
    header["AREASCAL"] = 1.0
    header["BACKSCAL"] = 1.0
    header["RESPFILE"] = response.rmf_path.name if response.rmf_path else (response.rsp_path.name if response.rsp_path else "NONE")
    header["ANCRFILE"] = response.arf_path.name if response.arf_path else "NONE"
    header["BACKFILE"] = background_file
    header["CORRFILE"] = "NONE"
    header["CORRSCAL"] = 1.0
    _metadata_history(header, response.metadata)
    primary = fits.PrimaryHDU()
    primary.header["TELESCOP"] = "COSI"
    primary.header["INSTRUME"] = "BTO"
    primary.header["DETNAM"] = response.detector_name
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fits.HDUList([primary, spectrum]).writeto(output_path, overwrite=overwrite, checksum=True)


def _simulate_bto_spectrum_impl(
    response,
    photon_model,
    exposure_s,
    *,
    background=None,
    poisson=True,
    random_seed=None,
    output_pha=None,
    overwrite=False,
):
    """Fold a source model, add background, optionally draw Poisson counts.

    A bare background array is interpreted as counts/s/channel.  A mapping can
    instead specify ``counts_per_channel`` explicitly.  Dead time is not applied.
    Existing ``output_pha`` files are protected unless ``overwrite=True``.
    """

    exposure_s = float(exposure_s)
    if not np.isfinite(exposure_s) or exposure_s <= 0:
        raise ValueError("exposure_s must be finite and positive")
    flux, source_counts = _fold_photon_spectrum(response, photon_model, exposure_s)
    background_counts, background_file = _background_counts(background, response, exposure_s)
    total = np.clip(source_counts + background_counts, 0.0, None)
    if poisson:
        simulated = np.random.default_rng(random_seed).poisson(total).astype(np.int64)
    else:
        simulated = total.copy()
    pha_path = None
    if output_pha is not None:
        pha_path = Path(output_pha).expanduser().resolve()
        _write_pha(
            pha_path,
            simulated,
            exposure_s,
            response,
            background_file,
            poisson=bool(poisson),
            overwrite=bool(overwrite),
        )
    return BTOSpectrumSimulation(
        expected_source_counts=source_counts,
        expected_background_counts=background_counts,
        total_expected_counts=total,
        simulated_counts=simulated,
        integrated_photon_flux=flux,
        exposure_s=exposure_s,
        channel_edges=response.channel_edges.copy(),
        measured_energy_edges_keV=response.measured_energy_edges_keV.copy(),
        response_metadata=dict(response.metadata),
        pha_path=pha_path,
    )


def _apply_nonparalyzable_deadtime(rate_hz: Any, deadtime_s: float = 20.0e-6) -> np.ndarray:
    """Apply ``r_observed = r_true/(1+r_true*tau)``; disabled unless called.

    The default 20 microseconds/event is the prototype-system value reported by
    Nagasawa et al.; it is a placeholder, not silently part of the response.
    """

    rate = np.asarray(rate_hz, dtype=np.float64)
    if np.any(~np.isfinite(rate)) or np.any(rate < 0):
        raise ValueError("rate_hz must be finite and nonnegative")
    if not np.isfinite(deadtime_s) or deadtime_s < 0:
        raise ValueError("deadtime_s must be finite and nonnegative")
    return rate / (1.0 + rate * float(deadtime_s))
