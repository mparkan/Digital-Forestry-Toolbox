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
| *LASextent.m*    | Compute spatial extent of 3D point clouds (LAS files) and write result to ESRI shapefile (.shp)| 

#### Grids

| Function        | Description                                                       | 
| --------------- | ----------------------------------------------------------------- | 
| *rasterize.m*   | Convert a 3D point cloud to a 2D/3D raster                        | 


## Data

#### Airborne Laser Scanning

| File                      | Location                                                          | 
| ------------------------- | ----------------------------------------------------------------- | 
| *zh_6995_2710.las*        |                                                                   | 
