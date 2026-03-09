# data

Data sets should be stored here and organized accordingly to the user needs.
E.g.
├── data
│   ├── climate
│   ├── soil
│   ├── land_use
    
or e.g.:
├── data
│   ├── raw            <- The original
│   ├── processed      <- Data that was processed for modeling.

Remember to never copy common input data here, but rather make symbolic links to the original path.
