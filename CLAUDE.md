# CLAUDE.md

Commissioning and detector-characterization tools for the LCO (Las Cumbres Observatory) telescope network: CCD/CMOS noise & gain (photon transfer curves), crosstalk, bad-pixel masks, focus curves, flat-field comparison, and scripts that submit engineering observations to the LCO scheduler. Also hosts a long-running "noise/gain crawler" that analyses archive data nightly and publishes gain-history plots via a small Flask web app.

Most code is single-author research tooling. Only the scripts registered in `pyproject.toml` `[project.scripts]` (and those described in README.md) are maintained; many other modules are vestigial.

## Setup and commands

```sh
pip install -e .                 # build is pyproject.toml/setuptools (no setup.py since 2.2.1)
venv/bin/python -m pytest tests/test_Image.py tests/test_noisegain.py   # offline tests
```

- A local `venv/` (Python 3.12) exists; Docker image uses `python:3.12`; `requires-python >= 3.10`.
- Dependencies are hard-pinned in `pyproject.toml` (e.g. `sqlalchemy==1.3.16`, `astropy==8.0.0rc2`). Code uses SQLAlchemy 1.3 APIs (`declarative_base` from `sqlalchemy.ext.declarative`) — don't write 2.x-style ORM code.
- Tests depend on FITS files in `testdata/` that are **gitignored** (`*.fits`, `*.fits.fz`); they exist only on the developer's machine.
  - `tests/test_noisegain.py` runs the installed `noisegainmef` console script via subprocess on `testdata/ep60/` and compares output `.dat` tables against the tracked references in `testdata/ep60/ptc_data_*.dat`. Requires `pip install -e .` first.
  - `tests/test_archive.py` hits the live LCO archive/OpenSearch — network-dependent.
- Version is set in `pyproject.toml`; record releases in `Changenotes` (newest first).

## Layout

- `lcocommissioning/` — main package (flat layout, no `src/`, despite what `Changenotes` 2.2.1 says).
  - `common/` — shared library code (see below).
  - Top level: detector analysis CLIs — `noisegainrawmef.py` (the core noise/gain/PTC tool), `crosstalk.py`, `bpm.py`, `flatfield_comparison.py`, dark-current / linearity scripts, `submitXtalkObservation.py`, `submitNameModeTest.py`.
  - Per-instrument subpackages: `focus/`, `cdk/` (Delta Rho / CDK), `muscat/`, `nres/`, `floyds/`, `sophia/`, `sbig/`, `archon/`, `cmostest/` (QHY CMOS; vendors the QHYCCD SDK + `.so` in `qhyccdpython/`), `gpstiming/`, `guide/`.
- `noisegaincrawler/` — archive crawler (`crawl_noisegain.py`), history plotting + S3 upload (`analysegainhistory.py`), Flask web app (`webapp.py`, `templates/`). Driven by `crawlgain.sh`; this is what the Docker image / k8s deployment runs.
- `sciopstools/`, `allocationtools/` — science-operations helpers (proposal review spreadsheets, RTI block submission). `sciopstools` is packaged; `allocationtools` is not.
- `notebooks/`, `experimental/` — ad-hoc Jupyter analysis; not part of the package.
- Root shell scripts (`archonlabptc.sh`, `sophialabptc_binned.sh`, `crawlgain.sh`, …) are lab/ops wrappers around `noisegainmef`/`crawlnoisegain`.
- Root `*.sqlite`, `*.pickle`, `plots/`, `flatcache/`, `2026B_*.txt` are local working data, not source.

## Key concepts

**`common/Image.py`** — `Image(filename_or_hdulist, alreadyopenedhdu=..., overscancorrect=..., trim=...)` loads all SCI/SPECTRUM/COMPRESSED_IMAGE extensions of a (possibly `.fz`) MEF into a 3-D `self.data` array, handling overscan and DATASEC trimming (derived from BANZAI). For single-extension images, `getccddata(ext)` can simulate four "virtual quadrant" extensions via `common.quadrantboundaries()` (used for Sophia / single-amp cameras).

**Noise/gain pipeline** — `noisegainrawmef.do_noisegain_for_fileset(files, database, args, frameidtranslationtable)`:
1. `sortinputfitsfiles()` classifies inputs as bias (`b00`) vs flat (`f00`, `x00`, …) using headers (`CONFMODE`, `CCDSTEMP`/`CCDATEMP`), rejects mismatched readmode/temperature, and pairs flats by `--sortby exptime|filterlevel`.
2. `common/ccd_noisegain.py` computes per-extension (or per-quadrant with `--quadrants`) gain and read noise from bias pairs + flat pairs.
3. Results go to a `noisegaindb` (`common/noisegaindb_orm.py`, SQLAlchemy; any SQLAlchemy URL — sqlite locally, PostgreSQL in production) and optionally PNG plots / `ptc_data_*.dat` tables.
- With `--useaws`, frames are downloaded through the LCO archive API by frame id instead of read from disk. The crawler always passes the archive-query result table as `frameidtranslationtable`.
- `--noreprocessing` skips a file set if any file already appears in the DB.

**Archive access** — `common/lco_archive_utilities.py`: OpenSearch queries (`opensearch.lco.global`) to find calibration frames, `download_from_archive(frameid)`, and an on-disk crawler for `/archive/engineering`.

**Observation submission scripts** (`submit*` entry points) build LCO request groups and send them via `common.common.send_request_to_portal()` (observation portal, normal scheduling) or `submit_request_group()` (direct submission to a specific telescope). **They only submit when `--CONFIRM` is passed**; otherwise they log the request and exit. Never add `--CONFIRM` when running these yourself — it schedules real telescope time.

**Newer shared helpers (2026 refactor, partially adopted)** — `common/logging_config.py` (`setup_logging`), `common/request_builder.py` (`RequestBuilder`, `InstrumentConfigFactory`), `common/instrument_config.py` (site/instrument/constraint defaults), `common/date_time_utils.py`, `common/focus_fitting.py` (`FocusCurveFitter`). Prefer these over duplicating logic when touching submit or focus scripts. Planning notes are in `REFACTORING_*.md` and `MIGRATION_CHECKLIST.md` (Phase 1 done; Phases 2–5 open).

## Configuration (environment variables)

- `VALHALLA_URL`, `VALHALLA_API_TOKEN` — observation portal / scheduler (used by submit scripts).
- `ARCHIVE_API_TOKEN` — LCO science archive downloads.
- `AWS_S3_BUCKET`, `AWS_ACCESS_KEY_ID`, `AWS_SECRET_ACCESS_KEY`, `AWS_DEFAULT_REGION` — gain-history plot storage for the crawler/web app.
- `crawlgain.sh`: `DB_HOST`, `DB_PORT`, `DB_NAME`, `DB_USER`, `DB_PASS`, `NDAYS`, `OUTPUTDIR`.

## Deployment

- `Dockerfile` installs the package into `python:3.12`; GitHub Actions (`.github/workflows/publish-container-image.yaml`) builds with Nix + Skaffold and pushes to GHCR, then auto-PRs updated image tags into `k8s/base/components/ghcr-images`. `Jenkinsfile` is a legacy LCO Docker build.
- Local k8s dev: Nix flake / direnv shell (`./develop.sh`), then Skaffold — see `DEVELOPMENT_K8s.md`.

## Conventions

- CLIs use `argparse` with `ArgumentDefaultsHelpFormatter`, a `parseCommandLine()` function, and a `main()` referenced from `pyproject.toml`. Add new tools the same way and register them under `[project.scripts]`.
- Module loggers: `_log = logging.getLogger(__name__)` (some files use `_logger` / `log`).
- LCO naming: site codes (`lsc`, `cpt`, `coj`, `elp`, `ogg`, `tfn`, …), camera ids by family — `fa` Sinistro, `ep` MuSCAT, `sq` QHY CMOS, `fl`/`fs`/`kb` legacy; filenames like `coj1m003-fa19-20211018-0010-b00.fits.fz` (`<site><telescope>-<camera>-<DAY-OBS>-<seq>-<type>`; `b00` bias, `d00` dark, `f00`/`x00` flat/exposure).
