# Digital-Forestry-Toolbox
A collection of digital forestry tools for Matlab.


## Scripts

#### Input/Output

| Function        | Description                                                       | 
| --------------- | ----------------------------------------------------------------- | 
| *LASread.m*     | Read 3D point cloud from ASPRS LAS (1.0-1.4) file (.las)          | 
| *LASwrite.m*    | Write 3D point cloud to ASPRS LAS (1.0-1.4) file (.las)           | 
| *HDRread.m*     | Read ENVI header file (.hdr)                                      | 
| *SBETread.m*    | Read Applanix SBET trajectory file (.sbet)                        | 
| *TRJread.m*     | Read Terrascan TRJ trajectory file (.trj)                         | 

#### Basic

| Function         | Description                                                       | 
| ---------------- | ----------------------------------------------------------------- | 
| *LASclip.m*      | Clip 3D point cloud (LAS file) with a polygon               |
| *LASmerge.m*     | Merge 3D point clouds (LAS files)                           | 
| *LASextent.m*    | Extract spatial extent of 3D point clouds (LAS files) and write result to ESRI shapefile (.shp)| 
| *crossSection.m* | Extract a 2D cross-section from a 3D point cloud| 

#### Grids

| Function        | Description                                                       | 
| --------------- | ----------------------------------------------------------------- | 
| *rasterize.m*   | Convert a 3D point cloud to a 2D/3D raster                        | 

#### Laser metrics

| Function           | Description                                                       | 
| ------------------ | ----------------------------------------------------------------- | 
| *laserEchoRatio.m* | Computes the echo ratio of a 3D point cloud                       | 
| *laserPulses.m*     | Determines the pulse number associated with the individual returns (sorted by acquisition GPS time) |
| *laserTimeFormat.m* | Convert GPS week time to GPS satellite time, GPS adjusted time or UTC time |
| *laserIntensityCorrection.m* | Corrects return intensity for the sensor to target range |


#### Individual tree crown segmentation

| Function           | Description                                                       | 
| ------------------ | ----------------------------------------------------------------- | 
| *treeRegionGrowing.m* |                                                                | 

## Data

#### Airborne Laser Scanning

| File (download link) | Land cover | Location | Sensor | Date |
| ------------------------- | ---------- | ----------------------------------------------------------------- | ------------ | ----------------- |
| [*6995_2710.laz*][1] | Mixed forest | Zürich, Switzerland (47.58384 N, 8.76594 E), [see map][2] | Trimble AX60 | March 10-13, 2014 |
| [*6850_2475.laz*][3] | Urban, mixed forest | Zürich, Switzerland (47.37460 N, 8.56824 E), [see map][4] | Trimble AX60 | April 2, 2014 |

#### National Forest Inventories

| File                      | Description                                                       |
| ------------------------- | ----------------------------------------------------------------- |
| *swiss_nfi_flora.xlsx*  | Swiss National Forest Inventory (NFI), tree species Latin/vernacular (English, French, German, Italian) names and codes |


[1]: http://maps.zh.ch/download/hoehen/2014/lidar/6995_2710.laz
[2]: https://map.geo.admin.ch/?topic=ech&lang=fr&bgLayer=ch.swisstopo.swissimage&layers=ch.swisstopo.zeitreihen,ch.bfs.gebaeude_wohnungs_register,ch.bafu.wrz-wildruhezonen_portal,ch.swisstopo.swisstlm3d-wanderwege&layers_visibility=false,false,false,false&layers_timestamp=18641231,,,&X=271212&Y=699817&zoom=10&crosshair=marker                                           

[3]: http://maps.zh.ch/download/hoehen/2014/lidar/6850_2475.laz
[4]: https://map.geo.admin.ch/?topic=ech&lang=fr&bgLayer=ch.swisstopo.swissimage&layers=ch.swisstopo.zeitreihen,ch.bfs.gebaeude_wohnungs_register,ch.bafu.wrz-wildruhezonen_portal,ch.swisstopo.swisstlm3d-wanderwege&layers_visibility=false,false,false,false&layers_timestamp=18641231,,,&X=247761&Y=685308&zoom=9&crosshair=marker
