# enriching_3dBAG_with_building_metrics
3D building metrics such as roof surface orientation, volunme, rectangularity, hemisphericality, roughness index of buildings are calculated and added to the 3D BAG file.  We use CityJSON file format.

## Implementation
 - sample data is given in data folder
 - to test new 3D BAG file, add the file in data folder and modify the file name in line 530 of main.cpp in src folder
 - all the code files are in src directory
 - the updated enriched CityJSON file name can be defined in line 537 of main.cpp in src folder, after implementing the code the output could be found in output directory