% *************************************************************************
% Extracting data from LICSAR frames
% Author: Alireza Sobouti(alireza.sobouti@icloud.com)
% Upadate Log: 
% 1/11/2021 : Second Column of Matrix C replaced by Calculation of range*sin(inc) 
% Before runing the code name the desired folder in which to save the file
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% Explaining Variables:
% foldername:
% frame_ID: 
% avg_height: 
% *************************************************************************
clc, clear, close all,fclose('all');
format long g
foldername = 'SBAS';
% start_date = '20141016';
frame_ID = '086A_06402_111313';
%159A_06375_101213
n_patch = 40; 

% ########## Extracting average height of satellite ##########
Data_Path = strcat('../LiCSAR/',frame_ID);
cd(Data_Path)
cd metadata
fm = fopen('metadata.txt');
while ~feof(fm)
    linef = fgetl(fm);
    if strcmp(linef(1:14),'centre_range_m')
        centre_range_m = str2double(linef(16:end));
        break
    end
end
clear fm line 
% ########## Raeding and saving information from metadate files ##########
[E,georef] = geotiffread(strcat(frame_ID,'.geo.E.tif'));
N = geotiffread(strcat(frame_ID,'.geo.N.tif'));
U = geotiffread(strcat(frame_ID,'.geo.U.tif'));
H = geotiffread(strcat(frame_ID,'.geo.hgt.tif')); % Reading DEM file
baseline = dlmread('baselines'); % Reading perpendicular baseline file
LatitudeLimits = georef.LatitudeLimits;
LongitudeLimits = georef.LongitudeLimits;
RasterSize = georef.RasterSize;
SampleSpacingInLatitude = georef.SampleSpacingInLatitude;
SampleSpacingInLongitude = georef.SampleSpacingInLongitude;
% ########################################################
% ########## Caculation of Geometric Parameters ##########
% ########################################################
incidence_angle = acotd(U./sqrt(E.^2+N.^2)); % Calculating incidence angle
ALD_angle = rem(atan2d(E,N)+360,360); % Calculating Azimuth Look Direction
 %Calculating distance of every pixel to satellite
centre_teta_m = incidence_angle(round(size(incidence_angle,1)/2),round(size(incidence_angle,2)/2));
height = centre_range_m*cosd(centre_teta_m);
range = height./cosd(incidence_angle);
r = RasterSize(1);
c = RasterSize(2);
pixels = r*c; % number of total pixels
clear avg_height RasterSize
% Get user input of the number of patches
fprintf('number of lines for every patch: %d\n',round(c/n_patch))
patch = r*(round(c/n_patch)); % number of pixels in every patch
patches = patch*ones(n_patch+1,1);
patches(1) = 1;
patches(end) = patch - (patch*n_patch-pixels); % number of pixels in last patch
cd ..
if exist(foldername,'dir') ~= 7
[~,msg,~] = rmdir(foldername,'s');
disp(msg)
[~,msg,~] = mkdir(foldername);
disp(msg)
end
cd(foldername)
phi0 = georef.LatitudeLimits(2);
lambda0 = georef.LongitudeLimits(1);
phi = phi0 - rem(1:pixels,r).*georef.SampleSpacingInLatitude;
lambda = lambda0 + floor((1:pixels)./r)*georef.SampleSpacingInLongitude;
fid = fopen('bad_interfero.txt','a+');
cd ../interferograms
d = dir;
n = char(d(3:end,1).name);
% s = sum([n(:,1:8) == start_date]');
% s = find(s==8);
% n = n(s(1):end,:);
ints = size(n,1);
patch_temp = nan(ints,1);
master_date = nan(ints,1);
slave_date = nan(ints,1);
interferos = 1:ints;
bad_interfero = [];
for p = 1:n_patch
    pp = sum(patches(1:p)):sum(patches(2:p+1));
    patch_coordinate = nan(length(pp),3);
    patch_data = nan(length(pp),3*ints);
    for jj = 1:length(interferos)
        ii = interferos(jj);
        cd(n(ii,:))
        try
            unw = geotiffread(strcat(n(ii,:),'.geo.unw.tif'));
            cc = geotiffread(strcat(n(ii,:),'.geo.cc.tif'));
            %diff_pha = geotiffread(strcat(n(ii,:),'.geo.diff_pha.tif'));
        catch
            bad_interfero = [bad_interfero;ii];
            fprintf(fid,'%i\n',ii);
            cd ..
            fprintf('bad interferogram: %i\n',ii)
            continue
        end
        m_year = str2double(n(ii,1:4));
        m_month = str2double(n(ii,5:6));
        m_day = str2double(n(ii,7:8));
        s_year = str2double(n(ii,10:13));
        s_month = str2double(n(ii,14:15));
        s_day = str2double(n(ii,16:17));
        ind_m = find(baseline(:,2)==str2double(n(ii,1:8)));
        ind_s = find(baseline(:,2)==str2double(n(ii,10:17)));
        perp_base = baseline(ind_s,3) - baseline(ind_m,3);
        if p == 1
            temp_base = my_daysact(strcat(num2str(m_month),'/',num2str(m_day),'/',num2str(m_year)),strcat(num2str(s_month),'/',num2str(s_day),'/',num2str(s_year))); %days
            %temp_base = daysact(datetime([m_year,m_month,m_day,0,0,0]),datetime([s_year,s_month,s_day,0,0,0])); %days
            patch_temp(ii,1) = temp_base;
            master_date(ii,1) = str2double(n(ii,1:8));
            slave_date(ii,1) = str2double(n(ii,10:17));
        end
        row = [pp',unw(pp)',range(pp').*sind(incidence_angle(pp')),double(cc(pp')),[1:length(pp)]'];
        row(isnan(row(:,2)),:) = [];
        row(row(:,4)==0,:) = [];
        patch_coordinate(row(:,5),1) = phi(row(:,1));
        patch_coordinate(row(:,5),2) = lambda(row(:,1));
        patch_coordinate(row(:,5),3) = H(row(:,1));
        patch_data(row(:,5),3*ii-2) = row(:,2);
        patch_data(row(:,5),3*ii-1) = row(:,3);
        patch_data(row(:,5),3*ii) = row(:,4);
        fprintf('Patch %i of %i & Interferogram %i of %i\n',[p n_patch ii length(n)])
        clear unw cc
        cd ..
    end
    cd ..
    cd(foldername)
    if p == 1
        save('temp_base.mat','patch_temp');
        save('master_date.mat','master_date');
        save('slave_date.mat','slave_date');
    end   
    patch_data(isnan(patch_coordinate(:,1)),:) = [];
    patch_coordinate(isnan(patch_coordinate(:,1)),:) = [];
    save(strcat('patch',num2str(p),'coordinate.mat'),'patch_coordinate');
    unw_phase = patch_data(:,1:3:end);
    save(strcat('patch',num2str(p),'unw.mat'),'unw_phase');
    Rsin_inc = patch_data(:,2:3:end);
    save(strcat('patch',num2str(p),'Rsin_inc.mat'),'Rsin_inc');
    coh = patch_data(:,3:3:end);
    save(strcat('patch',num2str(p),'coh.mat'),'coh');
    cd ../interferograms
    clear patch_coordinate patch_data
    interferos(bad_interfero) = [];
end
fclose('all');





function numDaysBetween = my_daysact(date0,date1)

dForm = 'mm/dd/yyyy'; % Date Format String, datenum does not necessarily need this
numDaysBetween = abs(datenum(date0,dForm) - datenum(date1,dForm));

end