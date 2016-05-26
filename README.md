# Digital-Forestry-Toolbox
The Digital Forestry Toolbox is collection of tools and tutorials for Matlab designed to help process datasets (e.g. remote sensing and field surveys) related to forests.


# Table of contents

+ [Scripts](#scripts)
  - [Input/Output](#scripts-1)
  - [Basic](#scripts-2)
  - [Grids](#scripts-3)
  - [Laser metrics](#scripts-4)
  - [Tree metrics](#scripts-5)
+ [Data](#data)
  - [Airborne Laser Scanning - 3D point clouds](#data-1)
  - [National Forest Inventories](#data-2)

## License

## Scripts <a id="scripts"></a>

#### Input/Output <a id="scripts-1"></a>

| Function        | Description                                                       | 
| --------------- | ----------------------------------------------------------------- | 
| [*LASread.m*](scripts/io/las/LASread.m)     | Read 3D point cloud from ASPRS LAS (1.0-1.4) file (.las)          | 
| [*LASwrite.m*](scripts/io/las/LASwrite.m)    | Write 3D point cloud to ASPRS LAS (1.0-1.4) file (.las)           | 
| [*HDRread.m*](scripts/io/envi/HDRread.m)    | Read ENVI header file (.hdr)                                      | 
| [*SBETread.m*](scripts/io/sbet/SBETread.m)    | Read Applanix SBET trajectory file (.sbet)                        | 
| [*TRJread.m*](scripts/io/trj/TRJread.m)     | Read Terrascan TRJ trajectory file (.trj)                         | 

#### Basic <a id="scripts-2"></a>

| Function         | Description                                                       | 
| ---------------- | ----------------------------------------------------------------- | 
| [*LASclip.m*](scripts/basic/LASclip.m)      | Clip 3D point cloud (LAS file) with a polygon               |
| [*LASmerge.m*](scripts/basic/LASmerge.m)     | Merge 3D point clouds (LAS files)                           | 
| [*LASextent.m*](scripts/basic/LASextent.m)    | Extract spatial extent of 3D point clouds (LAS files) and write result to ESRI shapefile (.shp) | 
| [*crossSection.m*](scripts/crossSection.m) | Extract a 2D cross-section from a 3D point cloud| 

#### Grids <a id="scripts-3"></a>

| Function        | Description                                                       | 
| --------------- | ----------------------------------------------------------------- | 
| [*rasterize.m*](scripts/grids/rasterize.m) | Convert a 3D point cloud to a 2D/3D raster | 
| [*elevationModels.m*](scripts/grids/elevationModels.m) | Compute elevation models (i.e. terrain, surface, height) from a 3D classified point cloud | 

#### Laser metrics <a id="scripts-4"></a>

| Function           | Description                                                       | 
| ------------------ | ----------------------------------------------------------------- | 
| [*laserEchoRatio.m*](scripts/laser metrics/laserEchoRatio.m) | Computes the echo ratio of a 3D point cloud                       | 
| [*laserPulses.m*](scripts/laser metrics/laserPulses.m)     | Determines the pulse number associated with the individual returns (sorted by acquisition GPS time) |
| [*laserTimeFormat.m*](scripts/laser metrics/laserTimeFormat.m) | Convert GPS week time to GPS satellite time, GPS adjusted time or UTC time |
| [*laserIntensityCorrection.m*](scripts/laser metrics/laserIntensityCorrection.m) | Corrects return intensity for the sensor to target range |


#### Tree metrics <a id="scripts-5"></a>

| Function           | Description                                                       | 
| ------------------ | ----------------------------------------------------------------- | 
| [*treeRegionGrowing.m*](scripts/tree metrics/treeRegionGrowing.m) |  Extracts individual tree crowns from a 3D point cloud using a modified version of the top down region growing algorithm described in Li et al. (2012) | 

## Data <a id="data"></a>

#### Airborne Laser Scanning - 3D point clouds <a id="data-1"></a>

| File | Land cover | Location | Sensor | Date |
| ------------------------- | ---------- | ----------------------------------------------------------------- | ------------ | ----------------- |
| [*6995_2710.laz* (external link)][1] | Mixed forest | Zürich, Switzerland (47.58384 N, 8.76594 E, ~380 m ASL), [see map][2] | Trimble AX60 | March 10-13, 2014 |
| [*6850_2475.laz* (external link)][3] | Urban, mixed forest | Zürich, Switzerland (47.37460 N, 8.56824 E, ~550 m ASL), [see map][4] | Trimble AX60 | April 2, 2014 |
| [*zh_2014_coniferous.laz*](data/measurements/vector/als/zh_2014_coniferous.laz) | Coniferous forest | Zürich, Switzerland (47.58378 N, 8.76275 E, ~380 m ASL), [see map][5] | Trimble AX60 | March 10-13, 2014 |

#### National Forest Inventories <a id="data-2"></a>

| File                      | Description                                                       |
| ------------------------- | ----------------------------------------------------------------- |
| [*swiss_nfi_flora.xlsx*](/data/reference/tabular/nfi/swiss_nfi_flora.xlsx) | Swiss National Forest Inventory (NFI), tree species Latin/vernacular (English, French, German, Italian) names and codes |


[1]: http://maps.zh.ch/download/hoehen/2014/lidar/6995_2710.laz
[2]: https://map.geo.admin.ch/?topic=ech&lang=fr&bgLayer=ch.swisstopo.swissimage&layers=ch.swisstopo.zeitreihen,ch.bfs.gebaeude_wohnungs_register,ch.bafu.wrz-wildruhezonen_portal,ch.swisstopo.swisstlm3d-wanderwege&layers_visibility=false,false,false,false&layers_timestamp=18641231,,,&X=271212&Y=699817&zoom=10&crosshair=marker                                           

[3]: http://maps.zh.ch/download/hoehen/2014/lidar/6850_2475.laz
[4]: https://map.geo.admin.ch/?topic=ech&lang=fr&bgLayer=ch.swisstopo.swissimage&layers=ch.swisstopo.zeitreihen,ch.bfs.gebaeude_wohnungs_register,ch.bafu.wrz-wildruhezonen_portal,ch.swisstopo.swisstlm3d-wanderwege&layers_visibility=false,false,false,false&layers_timestamp=18641231,,,&X=247761&Y=685308&zoom=9&crosshair=marker

[5]:https://map.geo.admin.ch/?topic=ech&lang=fr&bgLayer=ch.swisstopo.swissimage&layers=ch.swisstopo.zeitreihen,ch.bfs.gebaeude_wohnungs_register,ch.bav.haltestellen-oev,ch.swisstopo.swisstlm3d-wanderwege&layers_visibility=false,false,false,false&layers_timestamp=18641231,,,&X=271187&Y=699615&zoom=11&crosshair=marker
