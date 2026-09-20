# BTO APIs and Tutorials for COSIpy spectral fitting

This is a project for how to use the BTO in COSIpy spectral fitting.   
This includes a set of Jupyter notebooks that demonstrate how to construct BTO responses, generate fake observations, and compare COSI Compton-camera-only (CC) spectral fits with joint CC+BTO1+BTO2 fits.  
The notebooks use standard COSIpy likelihood classes, standard threeML `OGIPLike` plugins and the public `bto_response` module. 

## Notebook Examples

| Notebook                                                                        | What you learn                                                                                                                                                                |
| ------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| [01 — BTO API and response quick looks](notebooks/01_bto_response_api.ipynb)    | Inspect HDF5, plot directional area maps, apply detector effects, write/read OGIP, regrid background, simulate a fake spectrum and load it in threeML. No CC input is needed. |
| [02 — Fixed cross-normalization](notebooks/02_grb_cc_bto_fixed_crossnorm.ipynb) | Generate CC and both BTO fake datasets, then fit CC only and CC+BTO1+BTO2 with both BTO effective-area correction factors fixed to one.                                       |
| [03 — Free cross-normalization](notebooks/03_grb_cc_bto_free_crossnorm.ipynb)   | Run the same experiment while independently fitting BTO1-to-CC and BTO2-to-CC factors in 0.5–1.5. CC remains the reference, not a frozen spectrum.                            |

Use the installed `cosipy` kernel with COSIpy, HistPy, Astropy, NumPy/SciPy,
h5py, Matplotlib, pandas, astromodels and threeML. In Section 1, set only `project_path`. All input paths use its `response/`
subdirectory, and the BTO module is bundled in `api/`. 

## Input files

All paths below are relative to `project_path`, the directory containing this README.
The five scientific inputs are ordinary local copies (about 1.48 GB total), not symlinks.
Move the whole project and change only `project_path`; the installed Python environment is still required.

```text
project_path/
├── api/bto_response.py
├── response/
│   ├── SMEXv12.Continuum.HEALPixO3_10bins_log_flat.binnedimaging.imagingresponse.h5
│   ├── 20280301_3_month_with_orbital_info.fits
│   ├── bkg_binned_data_1s_local.hdf5
│   ├── bto_response_order2.h5
│   ├── BTO_background_template_10-3000keV.fits
├── notebooks/                # Jupyter notebooks for the tutorial
└── notebook_output/          # Generated products, not inputs
```


## Scientific example and COSIpy reference

This tutorial extends the COSIpy [GRB spectral-fitting example](https://cositools-cosipy.readthedocs.io/en/latest/tutorials/spectral_fits/continuum_fit/grb/SpectralFit_GRB.html)
([source notebook](https://github.com/cositools/cosipy/blob/main/docs/tutorials/spectral_fits/continuum_fit/grb/SpectralFit_GRB.ipynb)).
The reference example identifies its simulated burst as **GRB090206620** and
uses Galactic `(l,b)=(93,-53)` degrees. We adopt that simulated position and
injected Band spectrum, generate new response-consistent Poisson data, and add
the two BTO detectors. The 40 s interval starts at **2028-05-22 08:36:50 UTC**
in the simulated spacecraft history; this is not the historical GRB date or
a flight observation.

The supplied **CC background contains the albedo-photon simulation** described
by the COSIpy example. It is not a total satellite-background model. Cosmic
diffuse photons, charged particles, neutrons and activation/SAA components must
not be assumed present. The background file matches the official download
checksum. Its `PsiChi` axis is explicitly `spacecraftframe`, matching the local
coordinate convention used by the installed response implementation.

The BTO background is a different, 11-component simulation. Therefore the
CC-versus-joint comparison teaches joint fitting and the effect of BTO in this
specified example; it is **not a full-background mission sensitivity forecast**.



```python
from pathlib import Path
import sys

project_path = Path('/path/to/bto')
response_path = project_path / 'response'
bto_api_path = project_path / 'api'
sys.path.insert(0, str(bto_api_path))

response_h5_path = response_path / 'bto_response_order2.h5'
orientation_path = response_path / '20280301_3_month_with_orbital_info.fits'
background_template_path = response_path / 'BTO_background_template_10-3000keV.fits'
```

All inputs are read-only; notebook outputs go under this project's
`notebook_output/{fixed_crossnorm,free_crossnorm,bto_api}`.

| Variable                                                                    | Default file                                                                            | Meaning                                                                                                                                                                                                |
| --------------------------------------------------------------------------- | --------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `cc_response_path`                                                          | `response/SMEXv12.Continuum.HEALPixO3_10bins_log_flat.binnedimaging.imagingresponse.h5` | COSIpy's full simulated Compton response from incident direction/true energy to measured energy, scatter angle and scattered-photon direction. Neither source counts nor a 1D RMF. GRB notebooks only. |
| `orientation_path`                                                          | `response/20280301_3_month_with_orbital_info.fits`                                      | Simulated timestamps, X/Z pointing, Earth zenith, altitude and livetime. Transforms the celestial source into spacecraft directions and supplies interval exposure. Shared by CC and BTO.              |
| `cc_time_resolved_background_path`                                          | `response/bkg_binned_data_1s_local.hdf5`                                                | One-second albedo-photon CC background histogram in `Time, Em, Phi, PsiChi`. Supplies both measurement axes and OFF-background shape. GRB notebooks only.                                              |
| `bto_response_h5_path` (`response_h5_path` in API notebook)                 | `response/bto_response_order2.h5`                                                       | This project's all-direction, **unsmeared** Geant4 BTO1/BTO2 deposit response. The input to detector-model application and new OGIP generation.                                                        |
| `bto_background_template_path` (`background_template_path` in API notebook) | `response/BTO_background_template_10-3000keV.fits`                                      | Simulated mean single-BTO background spectrum in counts/s/channel. Constant in time for this example; used separately for both BTOs.                                                                   |

The CC inputs are distributed by COSIpy; the BTO response and BAK are project
products. COSIpy's [example source](https://github.com/cositools/cosipy/blob/main/docs/api/interfaces/examples/grb/example_grb_fit_threeml_plugin_interfaces.py)
lists their download locations. The local CC input checksums were verified:

| Input         | Distribution key below `COSI-SMEX/`                                                                            | MD5                                |
| ------------- | -------------------------------------------------------------------------------------------------------------- | ---------------------------------- |
| CC response   | `cosipy_tutorials/Data/Responses/SMEXv12.Continuum.HEALPixO3_10bins_log_flat.binnedimaging.imagingresponse.h5` | `eb72400a1279325e9404110f909c7785` |
| History       | `DC3/Data/Orientation/20280301_3_month_with_orbital_info.fits`                                                 | `5e69bc1d55fab9390f90635690f62896` |
| CC background | `cosipy_tutorials/grb_spectral_fit_local_frame/bkg_binned_data_1s_local.hdf5`                                  | `b842a7444e6fc1a5dd567b395c36ae7f` |

The measured CC axes come directly from the background: `Em` has 10 bins
from 100–10000 keV, `Phi` has 36 bins from 0–180 degrees, and `PsiChi` uses
NSIDE=8/RING (768 pixels) in spacecraft coordinates. The likelihood operates
on all `(10,36,768)` CDS bins. Count spectra are energy projections for plotting.


## Calibration and resolution

The GRB notebooks explicitly define a configurable placeholder:

```python
calibration = BTOCalibration(
    n_channels=4096,
    gain_quadratic_ch_per_keV2=0.0,
    gain_slope_ch_per_keV=1.0/0.8125,
    gain_offset_ch=0.0,
    lower_threshold_keV=30.0,
    upper_saturation_keV=3000.0,
)
resolution = Prototype_Resolution
```

Gain is `channel = a*E_measured² + b*E_measured + c`. The API notebook also
demonstrates a nonzero quadratic coefficient. Resolution is
`sigma(E) = sqrt(a0 + a1*E + a2*E²)` and `FWHM=2*sqrt(2*ln(2))*sigma`.
The approximate preset uses FWHM/E = 65/662, anchored to the prototype's
approximately 65 keV FWHM at 662 keV. It is not the full fitted coefficient set
or flight calibration. Replace these dataclasses when validated detector-specific
calibration becomes available; the likelihood/OGIP workflow remains the same.

## Background model and energy grid

`BTO_background_template_10-3000keV.fits` was produced by the project's `make_bto_background_ogip.py`
from the reference single-BTO simulation used in `opt_bkd_bto.ipynb`. It combines:

- cosmic photons and SAA protons;
- primary protons, electrons, positrons and alphas;
- albedo photons and neutrons;
- secondary protons, electrons and positrons.

The preparation applies the component scale factors and five-bin smoothing.
The result is a mean spectrum, not a Poisson realization or flight measurement.
It does not encode the full systematic covariance from smoothing, environmental
variation or component normalization.

The self-contained template stores:

| Extension  | Contents                                                                                          |
| ---------- | ------------------------------------------------------------------------------------------------- |
| EBOUNDS    | CHANNEL, E_MIN, E_MAX [keV]; 310 native logarithmic bins covering 10–3000 keV                     |
| SPECTRUM   | Mean RATE [counts/s/bin], COUNTS=RATE×EXPOSURE, and the 7,979,955 s simulation reference exposure |
| COMPONENTS | Individual mean component rates on the same grid                                                  |
| PROVENANCE | Component names, scale factors, input paths and SHA-256 checksums                                 |

COUNTS is a weighted/smoothed expectation, not raw integer Poisson events.
It must not be used to infer naive simulation statistical errors. No XSPEC
response association or fixed ADC grid is needed in this input template.

    background = get_bto_background(
        response, background_template=background_template_path,
    )

The API reads the energy edges from the file and integrates fractional bin overlap
onto the response grid. It contains no template-specific channel-count/energy-step
constants. 

The prepared background template is included in `response/`; no raw component
simulations or background-building script are needed for the tutorials. Its
PROVENANCE extension records the simulation inputs and scale factors.
The notebooks read this template and write a response-channelized OGIP BAK.

## API recipe and fitting

The main component calls are:

| Function                | Required inputs                                                        | Output                                                                                  |
| ----------------------- | ---------------------------------------------------------------------- | --------------------------------------------------------------------------------------- |
| `get_bto_response`      | HDF5, BTO ID, one direction mode                                       | Response arrays after applying the selected detector effects model; optional OGIP files |
| `write_ogip_response`   | Response, output prefix                                                | ARF, RMF, RSP and metadata JSON                                                         |
| `load_bto_response`     | RSP or RMF+ARF                                                         | Recorded response ready for reuse, without re-smearing                                  |
| `get_bto_background`    | Response; optional self-contained background template or explicit rate | Background counts/s/channel on the response grid                                        |
| `write_ogip_background` | Background, response, output path                                      | BAK on that grid                                                                        |
| `simulate_bto_spectrum` | Response, photon model, exposure                                       | Expectations, sampled counts, optional PHA                                              |

```python
response = get_bto_response(
    response_h5_path, 'BTO1',
    orientation_file=orientation_path, source_coord=source_coord,
    tstart=tstart, tstop=tstop,
    calibration=calibration, resolution=resolution, write_ogip=False,
)
write_ogip_response(response, output_path / 'bto1')
background = get_bto_background(response, background_template=background_template_path)
write_ogip_background(
    background, response, output_path / 'bto1_background.bak', exposure_s=40.,
)
fake = simulate_bto_spectrum(
    response, band, exposure_s=40., background=background,
    poisson=True, random_seed=36512, output_pha=output_path / 'bto1_fake.pha',
)
```

The GRB notebooks perform this separately for BTO1/BTO2, sharing the astrophysical
Band model with CC. Their BTO fits use measured energies 30–2500 keV and background
grouping through `OGIPLike`; CC uses COSIpy `BinnedThreeMLModelFolding`,
`FreeNormBinnedBackground`, `PoissonLikelihood` and `ThreeMLPluginInterface`.

CC source/background expectations, independent Poisson realizations, total ON
data and OFF template remain distinct in-memory quantities. The OFF shape is
the sum of 100 s before and 100 s after ON; fake means scale by ON/OFF duration,
while the fit normalizes that shape with a free rate [Hz] and ON livetime.
Finite-template errors and time-variation systematics are not a separate OFF
likelihood in this demonstration. The Band formula, parameters and units are
explained in the model cell.
