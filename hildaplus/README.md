
# hildaplus — HILDA+ Smoothing & Upscaling

# Description

LPJ-GUESS landcover input files derived from HILDA+ (smoothed and scaled up to 0.5deg).

This directory is part of the `landsymm_py` project and contains standalone scripts
for processing raw HILDA+ data into half-degree NetCDF/text files suitable for the
remapping stage.

The HILDA+ → LPJ-GUESS land-cover aggregation policy is **YAML-configurable**
via the files under `config/`. See `config/README.md` for the schema and
how to use named profiles or supply a custom mapping.

# Project Organization

```
hildaplus/
├── config/                     <- YAML-driven HILDA+ → LPJ-GUESS mapping config
│   ├── loader.py               <- Loader API used by the upscaling scripts
│   ├── lpjg_landcover_mapping.yaml   <- Default config (matches historical hardcoded mapping)
│   ├── profiles/               <- Named alternative profiles
│   │   ├── lpjg_v3_default.yaml
│   │   ├── lpjg_legacy_v1.yaml
│   │   └── lpjg_treecrops_as_forest.yaml
│   └── README.md               <- Config schema documentation
├── data_docs/                  <- Documentation for HILDA+ data
├── docs/                       <- Project documentation
│   └── DMP.md                  <- Data Management Plan
├── .gitignore
├── README.md                   <- This file
├── Readme_original.md          <- Original project README (historical)
├── reports/                    <- Reports (smoothing/upscaling diffs)
├── requirements.txt
└── scripts/                    <- Processing scripts
    ├── run_chain.sh            <- Pipeline orchestration (accepts --mapping-config / --mapping-profile)
    ├── hilda_smoothing/        <- Smoothing scripts
    └── hildaplus-upscaling/    <- Upscaling scripts (config-driven via hildaplus.config)
```

**Note**: HILDA+ input/output data is stored centrally at
`landsymm_py/data/geodata_py/HILDA+/data/` (not inside this directory).

# Data sharing and re-use 

This data can be reused by LEMG members.
Do not share this data or pass on the data without the author's approval.
Please contact the author, if you have questions (Martin Wittenbrink (martin.wittenbrink@kit.edu)).
