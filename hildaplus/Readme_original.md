### Specifications of __hildap_vGLOB-2.0_geotiff_wgs84__ and __hildap_vGLOB-2.0_geotiff_eckert4__

Data set: Land Use/Cover states (1) and transitions (2)

Thematic classes

* Land Use/Cover categories (**__states__**)  
  00 ocean  
  11 urban  
  22 annual crops  
  23 tree crops  
  24 agroforestry  
  33 pasture/rangeland  
  40 forest (unknown/other)  
  41 forest (evergreen, needle leaf)  
  42 forest (evergreen, broad leaf)  
  43 forest (deciduous, needle leaf)  
  44 forest (deciduous, broad leaf)  
  45 forest (mixed)  
  55 unmanaged grass/shrubland  
  66 sparse/no vegetation  
  77 water 
  99 no data
* Transitions between Land Use/Cover categories (**__transitions__**) with 4-digit codes (XXYY, where XX is LULC category of previous year and YY is LULC category of reference year)

  e.g.  
  1111 urban (stable)  
  1122 urban to annual crops2  
  1123 urban to tree crops  
  1124 urban to agroforestry  
  1133 urban to pasture  
  1140 urban to forest (unknown/other)  
  1141 urban to forest (evergreen, needle leaf)  
  1142 urban to forest (evergreen, broad leaf)  
  1143 urban to forest (deciduous, needle leaf)  
  1144 urban to forest (deciduous, broad leaf)  
  1155 urban to shrub/grassland  
  1166 urban to other land  
  ...  
  etc



### __...eckert4__

**(equal-area projection; Unit: metres)**

Format: GeoTIFF  
Projection: EPSG:54012 - World_Eckert_IV  
Spatial coverage/extent: -18000000, 9000000, 18000000, -9000000  
Temporal coverage: 1960-2020  
Spatial resolution: 1 km  
Temporal resolution: 1 year

### __...wgs84__

**(latitude longitude; Unit: degree)**

Format: GeoTIFF  
Projection: EPSG:4326 - WGS 84  
Spatial coverage/extent: -180, 90, 180, -90  
Temporal coverage: 1960-2020  
Spatial resolution: 0.01°  
Temporal resolution: 1 year