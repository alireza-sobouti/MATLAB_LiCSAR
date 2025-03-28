MATLAB LiCSAR Data Processing Script
This MATLAB script is designed to process interferometric data from the LiCSAR project. It extracts data from LiCSAR frames, reads and processes metadata, and performs geometric calculations. The data is read and processed in patches for memory efficiency, especially when dealing with large geospatial datasets. The script outputs various results, including unwrapped interferograms, range values, and coherence.

Requirements
MATLAB (R2019a or later is recommended)

GeoTIFF file format support for geotiffread

Input data: LiCSAR interferogram data in .geo.*.tif format

Script Overview
The script performs the following main tasks:

Extracts the average height of the satellite:

Reads metadata and calculates the satellite's average height.

Reads and saves metadata files:

Extracts important metadata from files like .geo.E.tif, .geo.N.tif, .geo.U.tif, and the DEM file.

Performs geometric calculations:

Computes the incidence angle and azimuth look direction.

Calculates the distance of every pixel to the satellite.

Processes data in patches:

Divides the data into smaller patches to improve memory efficiency.

For each patch, it reads unwrapped interferograms, calculates phase, range, and coherence, and saves the results.

Handles bad interferograms:

If any interferograms are unreadable or corrupted, they are logged in a bad_interfero.txt file.

Variables and Parameters
foldername: Name of the folder where the output files will be saved (e.g., 'SBAS').

frame_ID: Identifier of the frame to process (e.g., '086A_06402_111313').

n_patch: Number of patches to divide the data into for memory efficiency (e.g., 40).

Data_Path: Path to the directory containing the LiCSAR data.

start_date: Start date for filtering interferograms (not used directly in the script but can be useful if uncommented).

patches: Number of pixels in each patch for efficient memory management.

LatitudeLimits and LongitudeLimits: Geospatial limits of the frame.

baseline: Perpendicular baseline information for each interferogram.

Phi and Lambda: Latitude and longitude values for each pixel.

How to Use
Prepare the Data:

Place your LiCSAR interferogram data files and metadata in the correct folder.

Ensure that the required .geo.*.tif files (e.g., .geo.E.tif, .geo.N.tif, etc.) are present for the frame you're processing.

Set the Parameters:

Modify the foldername and frame_ID variables to match your data and desired output folder.

Adjust the n_patch variable to change the number of patches (recommended default is 40).

Run the Script:

Open MATLAB and navigate to the folder containing the script.

Run the script in MATLAB by typing the following in the command window:

matlab
Copy
Edit
run('LiCSAR_data_processing.m')
View the Output:

The script will generate multiple output files saved in the specified foldername. These files include:

patchXcoordinate.mat: Contains the geographic coordinates for each patch.

patchXunw.mat: Contains the unwrapped interferogram for each patch.

patchXRsin_inc.mat: Contains the calculated Rsin(incidence_angle) values.

patchXcoh.mat: Contains the coherence values for each patch.

Additional files such as temp_base.mat, master_date.mat, and slave_date.mat will also be created to store metadata.

Handling Bad Interferograms:

If any interferograms are corrupted or unreadable, they will be logged in the bad_interfero.txt file. These interferograms will be skipped during processing.

Example Output Files
patch1coordinate.mat: Contains the geographic coordinates (latitude, longitude, and height) for the first patch.

patch1unw.mat: Contains the unwrapped phase for the first patch.

patch1Rsin_inc.mat: Contains the range times sine of the incidence angle for the first patch.

patch1coh.mat: Contains the coherence values for the first patch.

Notes
Memory Efficiency: The data is processed in patches to avoid memory overflow when dealing with large interferogram files. The size of each patch is controlled by the n_patch variable.

Coordinate Conversion: The script computes the geographic coordinates (latitude and longitude) for each pixel based on the metadata and geographic reference system.

Baseline Calculation: The script computes the perpendicular baseline for each interferogram, which is necessary for phase unwrapping.

Troubleshooting
Error in reading files: Ensure that the required .geo.*.tif files are present in the same folder as the script.

Memory issues: If the script runs out of memory, try reducing the number of patches by adjusting the n_patch variable.

Bad interferograms: If a specific interferogram cannot be read, check that the file exists and is not corrupted.

License
This script is released under the MIT License.

