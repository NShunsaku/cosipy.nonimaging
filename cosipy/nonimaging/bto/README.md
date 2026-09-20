# BTO APIs and Tutorials for cosipy spectral fitting

This project provides APIs and tutorials for incorporating BTO data into joint spectral analyses with Ge Compton telescope using the cosipy and threeML frameworks.  
This includes a set of Jupyter notebooks that demonstrate how to construct BTO responses, generate fake observations, and compare COSI Compton-camera-only (CC) spectral fits with joint CC+BTO1+BTO2 fits.  
The notebooks use standard cosipy likelihood classes, standard threeML `OGIPLike` plugins and the public `bto_response` module. 

------- 

## Notebook Examples

| Notebook                                                                       | What you learn                                                                                                                                                                |
| ------------------------------------------------------------------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| [01: BTO API and response quick looks](notebooks/01_bto_response_api.ipynb)    | Inspect HDF5, plot directional area maps, apply detector effects, write/read OGIP, regrid background, simulate a fake spectrum and load it in threeML. No CC input is needed. |
| [02: Fixed cross-normalization](notebooks/02_grb_cc_bto_fixed_crossnorm.ipynb) | Generate CC and both BTO fake datasets, then fit CC only and CC+BTO1+BTO2 with both BTO effective-area correction factors fixed to one.                                       |
| [03: Free cross-normalization](notebooks/03_grb_cc_bto_free_crossnorm.ipynb)   | Run the same experiment while independently fitting BTO1-to-CC and BTO2-to-CC factors in 0.5–1.5. CC remains the reference, not a frozen spectrum.                            |

Use the installed python kernel with cosipy, HistPy, Astropy, NumPy/SciPy, h5py, Matplotlib, pandas, astromodels and threeML. In Section 1, set only `project_path`.   
All input paths use its `response/` subdirectory, and the BTO module is bundled in `api/`. 

## Input files

All paths below are relative to `project_path`, the directory containing this README.

```text
project_path/
├── api/bto_response.py
├── response/
│   ├── SMEXv12.Continuum.HEALPixO3_10bins_log_flat.binnedimaging.imagingresponse.h5
│   ├── 20280301_3_month_with_orbital_info.fits
│   ├── bkg_binned_data_1s_local.hdf5
│   ├── bto_response_order2.h5
│   └── BTO_background_template_10-3000keV.fits
├── notebooks/                # Jupyter notebooks for the tutorial
└── notebook_output/          # Generated products, not inputs
```


## Scientific example and cosipy reference

This tutorial extends the cosipy [GRB spectral-fitting example](https://cositools-cosipy.readthedocs.io/en/latest/tutorials/spectral_fits/continuum_fit/grb/SpectralFit_GRB.html)
([source notebook](https://github.com/cositools/cosipy/blob/main/docs/tutorials/spectral_fits/continuum_fit/grb/SpectralFit_GRB.ipynb)).  
The reference example identifies its simulated burst as **GRB090206620** and uses Galactic `(l,b)=(93,-53)` degrees. We adopt that simulated position and injected Band spectrum, generate new response-consistent Poisson data, and add the two BTO detectors. The 40 s interval starts at **2028-05-22 08:36:50 UTC** in the simulated spacecraft history; this is not the historical GRB date or a flight observation.

The supplied **CC background contains the albedo-photon simulation** described by the cosipy example. It is not a total satellite-background model. Cosmic diffuse photons, charged particles, neutrons and activation/SAA components must not be assumed present. 

The BTO background is a different, 11-component simulation. Therefore the CC-versus-joint comparison teaches joint fitting and the effect of BTO in this specified example; it is **not a full-background mission sensitivity forecast**.


```python

project_path = Path('/path/to/bto')
response_path = project_path / 'response'
bto_api_path = project_path / 'api'
sys.path.insert(0, str(bto_api_path))

response_h5_path = response_path / 'bto_response_order2.h5'
orientation_path = response_path / '20280301_3_month_with_orbital_info.fits'
background_template_path = response_path / 'BTO_background_template_10-3000keV.fits'
```

All inputs are read-only; 
Notebook outputs go under this project's `notebook_output/{fixed_crossnorm free_crossnorm,bto_api}`.  

| Variable                                                                    | Default file                                                                            | Meaning                                                                                                                                                                                   |
| --------------------------------------------------------------------------- | --------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `cc_response_path`                                                          | `response/SMEXv12.Continuum.HEALPixO3_10bins_log_flat.binnedimaging.imagingresponse.h5` | COSI Compton camera continum response used for tutorial.                                                                                                                                  |
| `orientation_path`                                                          | `response/20280301_3_month_with_orbital_info.fits`                                      | Satellite information including timestamps, X/Z pointing, Earth zenith, altitude and livetime. Transforms the celestial source into spacecraft directions and supplies interval exposure. |
| `cc_time_resolved_background_path`                                          | `response/bkg_binned_data_1s_local.hdf5`                                                | One-second albedo-photon CC background histogram in `Time, Em, Phi, PsiChi`. Supplies both measurement axes and OFF-background shape for tutorial.                                        |
| `bto_response_h5_path` (`response_h5_path` in API notebook)                 | `response/bto_response_order2.h5`                                                       | All-direction, **unsmeared** Meaglib/Cosima BTO1/BTO2 deposit response. The input to detector-model application and new OGIP generation.                                                  |
| `bto_background_template_path` (`background_template_path` in API notebook) | `response/BTO_background_template_10-3000keV.fits`                                      | Simulated mean single-BTO background spectrum in counts/s/channel. Constant in time for this example                                                                                      |

The CC inputs are distributed by cosipy Wasabi server.  
See cosipy lists for their download locations:  
https://github.com/cositools/cosipy/blob/main/docs/api/interfaces/examples/grb/example_grb_fit_threeml_plugin_interfaces.py
The BTO response and background file are not distributed; Please contact the BTO team for access.  


| Input             | Distribution key below `COSI-SMEX/`                                                                            | MD5                                |
| ----------------- | -------------------------------------------------------------------------------------------------------------- | ---------------------------------- |
| CC response       | `cosipy_tutorials/Data/Responses/SMEXv12.Continuum.HEALPixO3_10bins_log_flat.binnedimaging.imagingresponse.h5` | `eb72400a1279325e9404110f909c7785` |
| Satellite History | `DC3/Data/Orientation/20280301_3_month_with_orbital_info.fits`                                                 | `5e69bc1d55fab9390f90635690f62896` |
| CC background     | `cosipy_tutorials/grb_spectral_fit_local_frame/bkg_binned_data_1s_local.hdf5`                                  | `b842a7444e6fc1a5dd567b395c36ae7f` |

The measured CC axes come directly from the background: `Em` has 10 bins from 100–10000 keV, `Phi` has 36 bins from 0–180 degrees, and `PsiChi` uses NSIDE=8/RING (768 pixels) in spacecraft coordinates.  
The likelihood operates on all `(10,36,768)` CDS bins. Count spectra are energy projections for plotting.  


## BTO Background model and energy grid

`BTO_background_template_10-3000keV.fits` was produced from the reference single-BTO simulation (conducted by S.Takashima, S.Nagasawa and H.Yoneda).   
It combines:  
- cosmic photons and SAA protons;
- primary protons, electrons, positrons and alphas;
- albedo photons and neutrons;
- secondary protons, electrons and positrons.

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

----

## BTO API Quick Start Recipe

### Purpose of this workflow

The BTO API separates the analysis into three conceptually different stages:

1. **Construct a detector response for a specified BTO and source direction.**
   The master HDF5 contains the unsmeared Geant4 deposited-energy response for
   BTO1 and BTO2 over the simulated spacecraft directions. The API interpolates
   this library to the requested direction or averages it over a source track,
   then applies the selected gain, energy resolution, threshold and saturation.

2. **Place the response, background and simulated observation on one common
   measured-energy/channel grid.**
   The response defines the detector channels. The background template is
   conservatively regridded onto those channels, and an optional fake source
   spectrum is folded through the area-valued response.

3. **Write standard OGIP products for reuse and fitting.**
   The resulting ARF, RMF, RSP, BAK and PHA files can be loaded by threeML
   through `OGIPLike`. The BTO module prepares these detector products; the
   astrophysical parameter fit itself is performed by threeML.

The main data flow is:

```text
unsmeared all-sky BTO HDF5
    + BTO1 or BTO2
    + fixed spacecraft direction
      or celestial source + orientation history + time interval
    + calibration, resolution, threshold and saturation
                         |
                         v
                 get_bto_response
                         |
                         v
                 BTOResponseProduct
            ┌────────────┼─────────────┐
            |            |             |
            v            v             v
      ARF/RMF/RSP   background grid   source folding
            |            |             |
            |            v             |
            |   get_bto_background     |
            |            |             |
            |            v             |
            |        BTOBackground     |
            |            |             |
            v            v             v
      write_ogip_response          simulate_bto_spectrum
            |            |             |
            |   write_ogip_background  |
            |            |             |
            v            v             v
       .arf/.rmf/.rsp    .bak          .pha
             \____________|____________/
                          |
                          v
                   threeML OGIPLike
                          |
                          v
             spectral fit and error evaluation
```

### Main public functions

| Function                | Purpose                                                                                                                  | Required inputs                                                      | Important optional inputs                                                                          | In-memory output                                                                           | Optional file output                                                                  |
| ----------------------- | ------------------------------------------------------------------------------------------------------------------------ | -------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------ | ------------------------------------------------------------------------------------- |
| `get_bto_response`      | Select BTO1 or BTO2, interpolate or time-average the unsmeared HDF5 response, and apply the detector effects model       | Master response HDF5, BTO ID and exactly one direction specification | Calibration, resolution, threshold override, Earth-occultation option and OGIP-writing options     | `BTOResponseProduct` containing unsmeared deposit diagnostics and the recorded ARF/RMF/RSP | Can directly write `.arf`, `.rmf`, `.rsp` and `.response.json` when `write_ogip=True` |
| `write_ogip_response`   | Serialize an already constructed `BTOResponseProduct` into OGIP response files                                           | `BTOResponseProduct` and output prefix                               | `overwrite`                                                                                        | The same response object, with output path attributes populated                            | `.arf`, `.rmf`, `.rsp` and `.response.json`                                           |
| `load_bto_response`     | Reload a previously generated recorded-channel response without rebuilding it from the master HDF5                       | Combined RSP, or standalone RMF                                      | ARF path when loading a standalone RMF                                                             | `BTOResponseProduct` reconstructed from OGIP                                               | None                                                                                  |
| `get_bto_background`    | Put a time-independent background spectrum on exactly the measured-energy grid of one response                           | `BTOResponseProduct`                                                 | Self-contained background template, explicit per-channel rate, total rate, or compatible PHA input | `BTOBackground` in counts s$^{-1}$ channel$^{-1}$                                          | None                                                                                  |
| `write_ogip_background` | Write the mean background on the response channel grid as an OGIP BAK                                                    | `BTOBackground`, matching response and output filename               | Statistical exposure and `overwrite`                                                               | The same background object, with its output path populated                                 | `.bak`                                                                                |
| `simulate_bto_spectrum` | Integrate a photon model in true-energy bins, fold it through the RSP, add background and optionally draw Poisson counts | Response, photon model and exposure                                  | Background, Poisson switch, random seed and PHA output filename                                    | `BTOSpectrumSimulation` containing all source/background expectations and sampled counts   | Type-I `.pha`                                                                         |

### 1. Construct a source- and detector-specific response

```python
response = get_bto_response(
    response_h5_path,
    "BTO1",
    orientation_file=orientation_path,
    source_coord=source_coord,
    tstart=tstart,
    tstop=tstop,
    calibration=calibration,
    resolution=resolution,
    earth_occultation=False,
    write_ogip=False,
)
```

#### Purpose

`get_bto_response` is the central response-construction function. It starts
from the unsmeared all-direction HDF5 library and returns the response for one
named detector and one observing geometry.

BTO1 and BTO2 are stored as separate detector entries in the HDF5. The
function selects only the requested detector; the two responses are never
summed or exchanged internally.

#### Required inputs

- `response_h5`  
  Path to the master `bto_response_order2.h5` file. Its main array is

  ```text
  response/area_matrix_cm2[detector, direction, Etrue, Edeposit]
  ```

  It contains effective area as a function of true photon energy and deposited
  energy. It does **not** yet contain measured-channel energy resolution, gain,
  threshold or saturation.

- `bto_id`  
  Detector identifier. Accepted values are `1`, `2`, `"BTO1"` and `"BTO2"`.

- Exactly one of the following mutually exclusive direction modes:

  1. Fixed spacecraft angles:

     ```python
     theta_deg=...,
     phi_deg=...,
     ```

     Here, `theta_deg` is the colatitude from spacecraft $+Z$, and `phi_deg`
     is the azimuth measured from $+X$ toward $+Y$.

  2. Fixed spacecraft unit vector:

     ```python
     spacecraft_vector=(x, y, z),
     ```

  3. Celestial source tracked through the spacecraft history:

     ```python
     orientation_file=orientation_path,
     source_coord=source_coord,
     tstart=tstart,
     tstop=tstop,
     ```

     `source_coord` is normally an Astropy `SkyCoord`. For every spacecraft
     history sample in the selected interval, the source is transformed into
     the spacecraft frame. The direction-dependent responses are then averaged
     using the corresponding time/livetime weights.

Providing an incomplete direction mode or combining multiple modes raises an
error instead of silently choosing one.

#### Detector-effects inputs

- `calibration`  
  A `BTOCalibration` object defining the channel count, quadratic gain,
  threshold and saturation. The gain convention is

  $$
  {\rm channel}
  = aE_{\rm measured}^{2}+bE_{\rm measured}+c .
  $$

- `resolution`  
  A `BTOResolution` object defining

  $$
  \sigma(E)=\sqrt{a_0+a_1E+a_2E^2},
  \qquad
  {\rm FWHM}(E)=2\sqrt{2\ln2}\,\sigma(E).
  $$

- `threshold_keV`  
  Optional one-call override of the lower threshold stored in `calibration`.

- `earth_occultation`  
  If enabled, requests Earth-occultation filtering in orientation mode.
  This option must not be interpreted as a general Earth-response correction.

The current example calibration and `Prototype_Resolution` are configurable
placeholders, not validated BTO1/BTO2 flight calibration.

#### Direction interpolation

For each requested spacecraft direction, the response uses the nearest four
available simulated directions and inverse-angular-distance-squared weights.
A query exactly at a simulated direction center reproduces that HDF5 response,
within numerical precision.

In orientation mode, this angular interpolation is evaluated for every source
direction sample before the time/livetime-weighted average is calculated.

#### Returned object

The returned `BTOResponseProduct` contains, among other fields:

- `photon_energy_edges_keV[n_photon + 1]`  
  True photon-energy bin boundaries.

- `deposit_energy_edges_keV[n_deposit + 1]`  
  Geant4 deposited-energy bin boundaries.

- `measured_energy_edges_keV[n_channel + 1]`  
  Measured-energy boundaries corresponding to the detector channels.

- `unsmeared_area_matrix_cm2[n_photon, n_deposit]`  
  Direction-interpolated Geant4 area matrix before detector effects.

- `unsmeared_arf_cm2[n_photon]`  
  Sum of the unsmeared area matrix over all deposited-energy bins:

  $$
  A_{\rm unsmeared}(E_i)
  =\sum_k A_{\rm deposit}(E_i,E_{{\rm dep},k}).
  $$

- `photopeak_arf_cm2[n_photon]`  
  A full-energy-deposit proxy based on the deposit-energy bins. This is a
  diagnostic quantity, not a selection using Geant4 interaction-process labels.

- `arf_cm2[n_photon]`  
  Recorded effective area after redistribution and losses below threshold,
  above saturation, or outside the available channel range.

- `rmf_probability[n_photon, n_channel]`  
  Conditional measured-channel distribution for events that remain recorded.
  Rows with nonzero ARF sum to one.

- `rsp_cm2[n_photon, n_channel]`  
  Full area-valued response used to fold a photon spectrum:

  $$
  {\rm RSP}_{ij}
  ={\rm ARF}_i\,{\rm RMF}_{ij}.
  $$

- `metadata`  
  Detector name, direction mode, source/time information, interpolation
  details, input checksums, calibration, resolution and software provenance.

A practical pattern is to construct the response with `write_ogip=False`,
inspect it, and explicitly call `write_ogip_response` in the next step. This
keeps the computational and file-writing stages visible in the notebook.

### 2. Write the ARF, RMF and RSP

```python
write_ogip_response(
    response,
    output_path / "bto1",
    overwrite=True,
)
```

This writes four related products:

```text
bto1.arf
bto1.rmf
bto1.rsp
bto1.response.json
```

Their meanings are:

- **ARF**  
  Stores `arf_cm2`, the recorded effective area as a function of true photon
  energy.

- **RMF**  
  Stores `rmf_probability`, the conditional redistribution from true photon
  energy to measured channel. For every nonzero-ARF row,

  $$
  \sum_j {\rm RMF}_{ij}=1.
  $$

- **RSP**  
  Stores `rsp_cm2`, the combined area-valued matrix. Unlike an RMF, its rows
  are not normalized to one:

  $$
  \sum_j {\rm RSP}_{ij}={\rm ARF}_i.
  $$

- **Response metadata JSON**  
  Preserves information that is not fully represented by the OGIP tables,
  including the source direction or orientation interval, HDF5 checksum,
  calibration and resolution provenance.

The function returns the same `BTOResponseProduct`, after setting:

```python
response.arf_path
response.rmf_path
response.rsp_path
response.metadata_path
```

Writing the response before writing the BAK or PHA is recommended. Those later
OGIP products can then record the appropriate response filenames in their
headers.

Existing files are protected unless `overwrite=True`.

### 3. Reload a previously generated response

A combined RSP can be loaded directly:

```python
reloaded = load_bto_response(
    output_path / "bto1.rsp",
)
```

A separated RMF and ARF pair can be loaded as:

```python
reloaded = load_bto_response(
    output_path / "bto1.rmf",
    arf_file=output_path / "bto1.arf",
)
```

The returned object again satisfies

```python
reloaded.rsp_cm2 == (
    reloaded.arf_cm2[:, None] * reloaded.rmf_probability
)
```

The detector effects have already been applied to the OGIP product. Loading it
does **not** apply gain, resolution or threshold a second time.

OGIP response files do not contain the original Geant4 deposited-energy
matrix. Therefore, for a loaded response:

```python
reloaded.metadata["unsmeared_components_available"] == False
```

and the unsmeared deposit arrays are placeholders. Use the original HDF5 and
`get_bto_response` when deposit-level diagnostics, a different source/time
interval or a different calibration are required.

### 4. Regrid the BTO background

```python
background = get_bto_background(
    response,
    background_template=background_template_path,
)
```

#### Purpose

The source response defines the measured detector channels, but the reference
background template may have a different energy grid. `get_bto_background`
integrates the template by fractional energy-bin overlap and returns a mean
background rate on exactly the response grid.

The measured background is **not** passed through the source RMF. It is already
a measured-energy spectrum and is only regridded.

#### Input choices

Supply the response and at most one background source:

- `background_template=path`  
  Recommended for this tutorial. The FITS file contains both a `SPECTRUM`
  extension and its own `EBOUNDS`. `RATE`, or `COUNTS/EXPOSURE`, is converted
  into counts s$^{-1}$ per native bin and conservatively regridded.

- `rate_per_channel=array_or_scalar`  
  Explicit background rate in counts s$^{-1}$ channel$^{-1}$. An array must
  already match the response channel count. A scalar is repeated over channels.

- `total_rate_hz=value`  
  Total background rate distributed uniformly over channels that overlap the
  acquisition threshold/saturation interval. This is mainly a simplified test
  option, not the preferred scientific background model.

- No explicit background argument  
  Uses the project’s bundled 10–3000 keV default background template.

`background_pha`, `background_rmf` and
`background_energy_edges_keV` remain available for compatibility with older
background files that do not have a self-contained energy grid.

The template is not extrapolated beyond its stored energy range. Target
channels outside that range receive zero contribution.

#### Returned object

`BTOBackground` contains:

- `rate_per_channel[n_channel]` in counts s$^{-1}$ channel$^{-1}$;
- `measured_energy_edges_keV[n_channel + 1]`, identical to the response grid;
- provenance and total-rate metadata;
- `path`, which remains `None` until an OGIP BAK is written.

This object represents a mean rate. It does not yet contain a Poisson
realization.

### 5. Write the OGIP background

```python
write_ogip_background(
    background,
    response,
    output_path / "bto1_background.bak",
    exposure_s=exposure_s,
    overwrite=True,
)
```

The response and background must have the same measured-energy boundaries and
number of channels.

The output BAK stores:

- `CHANNEL`;
- mean `RATE` in counts s$^{-1}$;
- `STAT_ERR`;
- `QUALITY`, `GROUPING`, `AREASCAL` and `BACKSCAL`;
- detector and response metadata.

`exposure_s` is used to construct the Poisson-equivalent statistical error,

$$
{\rm STAT\_ ERR}_j
=\frac{\sqrt{R_j\,T}}{T},
$$

where $R_j$ is the mean rate in channel $j$ and $T$ is the supplied
exposure. It does not Poisson-randomize the mean background.

After writing:

```python
background.path
```

points to the generated BAK. This path can be passed directly to `OGIPLike`.

### 6. Fold a photon model and generate a fake spectrum

```python
fake = simulate_bto_spectrum(
    response,
    photon_model,
    exposure_s,
    background=background,
    poisson=True,
    random_seed=36512,
    output_pha=output_path / "bto1_fake.pha",
    overwrite=True,
)
```

#### Photon-model input

`photon_model` may be either:

- a callable returning the differential photon spectrum

  $$
  \frac{dN}{dE} \quad [{\rm ph\ cm^{-2}\ s^{-1}\ keV^{-1}}],
  $$

  evaluated at energies in keV; or

- an array containing the already integrated photon flux in each true-energy
  bin, in ph cm$^{-2}$ s$^{-1}$ per bin.

For a callable, the function numerically integrates the model inside every
true-energy bin. It does not simply multiply a bin-center value by the bin
width.

The expected source counts in measured channel $j$ are

$$
\mu_{{\rm src},j} = T\sum_i F_i\,{\rm RSP}_{ij},
$$

where $F_i$ is the photon flux integrated over true-energy bin $i$, and
$T$ is `exposure_s`.

If a `BTOBackground` is supplied,

$$
\mu_{{\rm bkg},j}=T R_{{\rm bkg},j},
$$

and

$$
\mu_{{\rm total},j} = \mu_{{\rm src},j} + \mu_{{\rm bkg},j}.
$$

With `poisson=True`, the simulated channel counts are drawn as

$$
N_j\sim{\rm Poisson}(\mu_{{\rm total},j}).
$$

`random_seed` makes this realization reproducible. No dead-time correction,
pile-up or additional electronics effect is silently applied.

#### Returned object

`BTOSpectrumSimulation` contains:

- `integrated_photon_flux[n_photon]`;
- `expected_source_counts[n_channel]`;
- `expected_background_counts[n_channel]`;
- `total_expected_counts[n_channel]`;
- `simulated_counts[n_channel]`;
- `exposure_s`;
- channel and measured-energy boundaries;
- a copy of the response metadata;
- `pha_path`, if a PHA was written.

These arrays are useful for quick-look plots and validation. In particular,
expected source, expected background and simulated total counts remain
separate; the PHA contains the simulated total observation.

When `output_pha` is given, the function writes an OGIP Type-I spectrum. It is
best to write the response and background first, so the PHA headers can refer
to their filenames.

### 7. Create the threeML BTO plugin

The separated ARF and RMF can be supplied to the installed threeML
`OGIPLike` interface:

```python
from threeML.plugins.OGIPLike import OGIPLike

bto_plugin = OGIPLike(
    "bto1",
    observation=fake.pha_path,
    background=background.path,
    response=str(response.rmf_path),
    arf_file=str(response.arf_path),
    verbose=False,
)

bto_plugin.set_active_measurements("30-2500")
bto_plugin.rebin_on_background(20)
```

Here:

- `observation` is the simulated or measured Type-I PHA;
- `background` is the response-grid BAK;
- `response` is the conditional RMF;
- `arf_file` supplies the recorded effective area;
- `set_active_measurements("30-2500")` restricts the likelihood to measured
  energies from 30 to 2500 keV;
- `rebin_on_background(20)` groups neighboring channels according to the
  expected background counts.

The active fitting interval is conceptually distinct from the detector
threshold and saturation:

- threshold/saturation are applied while constructing the response and
  determine which events can be recorded;
- `set_active_measurements` selects which already-recorded channels are used
  by the likelihood.

For a joint CC+BTO fit, construct one response, background, fake/observed PHA
and `OGIPLike` plugin separately for BTO1 and BTO2. Assign both BTO plugins to
the same astrophysical source used by the CC plugin:

```python
bto1_plugin.assign_to_source(source_name)
bto2_plugin.assign_to_source(source_name)
```

The source model parameters are then shared across CC, BTO1 and BTO2. Optional
BTO-to-CC cross-normalization parameters belong to the fitting layer, not to
the response matrix itself.

### Complete minimal recipe

```python
from pathlib import Path
import sys

from astropy.coordinates import SkyCoord
from astropy.time import Time
import astropy.units as u
from astromodels import Band
from threeML.plugins.OGIPLike import OGIPLike

project_path = Path(
    "/path/to/project/"
)
response_path = project_path / "response"
output_path = project_path / "notebook_output" / "example"
output_path.mkdir(parents=True, exist_ok=True)

sys.path.insert(0, str(project_path / "api"))

from bto_response import (
    BTOCalibration,
    Prototype_Resolution,
    get_bto_response,
    write_ogip_response,
    get_bto_background,
    write_ogip_background,
    simulate_bto_spectrum,
)

response_h5_path = response_path / "bto_response_order2.h5"
orientation_path = response_path / "20280301_3_month_with_orbital_info.fits"
background_template_path = (
    response_path / "BTO_background_template_10-3000keV.fits"
)

source_coord = SkyCoord(l=93.0*u.deg, b=-53.0*u.deg, frame="galactic")
tstart = Time("2028-05-22T08:36:50")
tstop = tstart + 40.0*u.s
exposure_s = 40.0

calibration = BTOCalibration(
    n_channels=4096,
    gain_quadratic_ch_per_keV2=0.0,
    gain_slope_ch_per_keV=1.0/0.8125,
    gain_offset_ch=0.0,
    lower_threshold_keV=30.0,
    upper_saturation_keV=3000.0,
)

# 1. Build the direction- and interval-specific detector response.
response = get_bto_response(
    response_h5_path,
    "BTO1",
    orientation_file=orientation_path,
    source_coord=source_coord,
    tstart=tstart,
    tstop=tstop,
    calibration=calibration,
    resolution=Prototype_Resolution,
    earth_occultation=False,
    write_ogip=False,
)

# 2. Write the response products before the PHA so that their paths are known.
write_ogip_response(
    response,
    output_path / "bto1",
    overwrite=True,
)

# 3. Regrid the mean background onto this response's measured channels.
background = get_bto_background(
    response,
    background_template=background_template_path,
)

# 4. Write a background file with the same channel grid and exposure.
write_ogip_background(
    background,
    response,
    output_path / "bto1_background.bak",
    exposure_s=exposure_s,
    overwrite=True,
)

# 5. Define the incident photon model.
photon_model = Band()
photon_model.beta.value = -3.0
photon_model.beta.min_value = -20.0
photon_model.beta.max_value = -2.01
photon_model.K.value = 0.030210450807
photon_model.alpha.value = -0.360
photon_model.xp.value = 472.34624
photon_model.beta.value = -11.921
photon_model.piv.value = 500.0

# 6. Fold the model, add background and draw a reproducible fake spectrum.
fake = simulate_bto_spectrum(
    response,
    photon_model,
    exposure_s,
    background=background,
    poisson=True,
    random_seed=36512,
    output_pha=output_path / "bto1_fake.pha",
    overwrite=True,
)

# 7. Load the products as a standard threeML OGIP plugin.
bto_plugin = OGIPLike(
    "bto1",
    observation=fake.pha_path,
    background=background.path,
    response=str(response.rmf_path),
    arf_file=str(response.arf_path),
    verbose=False,
)
bto_plugin.set_active_measurements("30-2500")
bto_plugin.rebin_on_background(20)
```

For BTO2, repeat the same sequence with `"BTO2"` and a different output prefix.
Do not reuse the BTO1 response object for BTO2, because the two detector
geometries and direction-dependent effective areas are stored separately in
the master HDF5.
