
# hildaplus — HILDA+ Smoothing & Upscaling

# Description

LPJ-GUESS landcover input files derived from HILDA+ (smoothed and scaled up to 0.5deg).

This directory is part of the `landsymm_py` project and contains standalone scripts
for processing raw HILDA+ data into half-degree NetCDF/text files suitable for the
remapping stage.

# Project Organization

```
hildaplus/
├── data_docs/                  <- Documentation for HILDA+ data
├── docs/                       <- Project documentation
│   └── DMP.md                  <- Data Management Plan
├── .gitignore
├── README.md                   <- This file
├── Readme_original.md          <- Original project README (historical)
├── reports/                    <- Reports (smoothing/upscaling diffs)
├── requirements.txt
└── scripts/                    <- Processing scripts
    ├── run_chain.sh            <- Pipeline orchestration
    ├── hilda_smoothing/        <- Smoothing scripts
    └── hildaplus-upscaling/    <- Upscaling scripts
```

**Note**: HILDA+ input/output data is stored centrally at
`landsymm_py/data/geodata_py/HILDA+/data/` (not inside this directory).

# Data sharing and re-use 

This data can be reused by LEMG members.
Do not share this data or pass on the data without the author's approval.
Please contact the author, if you have questions (Martin Wittenbrink (martin.wittenbrink@kit.edu)).
