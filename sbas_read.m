% *************************************************************************
% Extracting data from LiCSAR frames
% Author: Alireza Sobouti (alireza.sobouti@icloud.com)
% Before running the code, specify the desired folder in which to save the file.
% *************************************************************************

% Clear workspace and close all files/figures
clc; clear; close all; fclose('all');
format long g

% ------------------------- USER-DEFINED PARAMETERS -------------------------
foldername = 'SBAS';               % Name of the folder to save results
frame_ID = '020D_06347_151515';    % Frame ID of the LiCSAR dataset
n_patch = 40;                      % Number of patches for processing
% ---------------------------------------------------------------------------

% ----------------------------- PATH SETTINGS -------------------------------
Data_Path = strcat('../LiCSAR/', frame_ID); % Path to the data folder
cd(Data_Path)
cd metadata

% ----------------------- EXTRACT SATELLITE PARAMETERS ----------------------
% Extract the center range of the satellite from metadata
fm = fopen('metadata.txt');
while ~feof(fm)
    linef = fgetl(fm);
    if contains(linef, 'centre_range_m')
        centre_range_m = str2double(linef(16:end));
        break
    end
end
fclose(fm);

% --------------------- READ GEOSPATIAL DATA AND METADATA -------------------
[E, georef] = geotiffread(strcat(frame_ID, '.geo.E.tif'));  % East component
N = geotiffread(strcat(frame_ID, '.geo.N.tif'));            % North component
U = geotiffread(strcat(frame_ID, '.geo.U.tif'));            % Up component
H = geotiffread(strcat(frame_ID, '.geo.hgt.tif'));          % DEM data
baseline = dlmread('baselines');                            % Baseline data

% Extract georeferencing details
LatitudeLimits = georef.LatitudeLimits;
LongitudeLimits = georef.LongitudeLimits;
RasterSize = georef.RasterSize;
SampleSpacingInLatitude = georef.SampleSpacingInLatitude;
SampleSpacingInLongitude = georef.SampleSpacingInLongitude;

% ------------------------ CALCULATE GEOMETRIC PARAMETERS -------------------
% Calculate incidence angle and azimuth look direction
incidence_angle = acotd(U ./ sqrt(E.^2 + N.^2));
ALD_angle = mod(atan2d(E, N) + 360, 360);

% Calculate range (distance to satellite) for each pixel
centre_theta_m = incidence_angle(round(end/2), round(end/2)); % Center incidence angle
height = centre_range_m * cosd(centre_theta_m);
range = height ./ cosd(incidence_angle);

% Calculate number of pixels
[r, c] = size(E);
pixels = r * c;

% -------------------------- PATCH PARAMETERS -------------------------------
% Divide data into patches for processing
fprintf('Number of lines per patch: %d\n', round(c / n_patch));
patch = r * round(c / n_patch);
patches = [1; (1:n_patch)' * patch];
patches(end) = patch - (patch * n_patch - pixels);

% ------------------- PREPARE OUTPUT DIRECTORY ------------------------------
cd ..
if ~exist(foldername, 'dir')
    mkdir(foldername);
end
cd(foldername)

% ----------------------- GENERATE COORDINATE DATA --------------------------
phi0 = LatitudeLimits(2);
lambda0 = LongitudeLimits(1);
phi = phi0 - rem(1:pixels, r) .* SampleSpacingInLatitude;
lambda = lambda0 + floor((1:pixels) ./ r) * SampleSpacingInLongitude;

% ------------------- PROCESS INTERFEROGRAMS PER PATCH ----------------------
cd ../interferograms
dirData = dir;
names = char(dirData(3:end).name);

% Filter interferograms based on start date
valid_idx = startsWith(names(:, 1:8), start_date);
names = names(valid_idx, :);

% Initialize variables
ints = size(names, 1);
bad_interfero = [];
fid = fopen('bad_interfero.txt', 'a+');

for p = 1:n_patch
    patch_idx = patches(p):patches(p+1) - 1; % Indices for the current patch
    patch_coordinate = nan(numel(patch_idx), 3); % Coordinate storage
    patch_data = nan(numel(patch_idx), 3 * ints); % Data storage
    
    for ii = 1:ints
        cd(names(ii, :));
        try
            % Read unwrapped phase and coherence data
            unw = geotiffread(strcat(names(ii, :), '.geo.unw.tif'));
            cc = geotiffread(strcat(names(ii, :), '.geo.cc.tif'));
        catch
            % Handle missing or invalid interferograms
            bad_interfero = [bad_interfero; ii];
            fprintf(fid, '%i\n', ii);
            cd ..;
            fprintf('Bad interferogram: %i\n', ii);
            continue;
        end
        
        % Calculate perpendicular baseline and dates
        ind_m = find(baseline(:, 2) == str2double(names(ii, 1:8)));
        ind_s = find(baseline(:, 2) == str2double(names(ii, 10:17)));
        perp_base = baseline(ind_s, 3) - baseline(ind_m, 3);

        % Store results for the current patch
        row = [patch_idx', unw(patch_idx)', range(patch_idx)' .* sind(incidence_angle(patch_idx))', cc(patch_idx)'];
        row(isnan(row(:, 2)), :) = [];
        patch_coordinate(:, :) = [phi(row(:, 1)), lambda(row(:, 1)), H(row(:, 1))];
        patch_data(:, :) = [row(:, 2), row(:, 3), row(:, 4)];
        
        cd ..;
    end
    
    % Save patch data
    save(strcat('patch', num2str(p), 'coordinate.mat'), 'patch_coordinate');
    save(strcat('patch', num2str(p), 'unw.mat'), 'patch_data(:, 1:3:end)');
    save(strcat('patch', num2str(p), 'Rsin_inc.mat'), 'patch_data(:, 2:3:end)');
    save(strcat('patch', num2str(p), 'coh.mat'), 'patch_data(:, 3:3:end)');
end

fclose(fid);

% ------------------------- HELPER FUNCTION ---------------------------------
function numDaysBetween = my_daysact(date0, date1)
    dForm = 'mm/dd/yyyy';
    numDaysBetween = abs(datenum(date0, dForm) - datenum(date1, dForm));
end
