%% ************************************************************************************************
%%               Time Series and Velocity Estimation from LiCSAR Unwrapped Interferograms
%% ************************************************************************************************
% This script performs time series and velocity estimation from LiCSAR unwrapped interferograms.
% Key features include:
% - Coherence filtering applied at the interferogram level
% - Options for SBAS network plotting
% - Support for plate motion models and estimation of residual effects
% Changes:
% - Filtering based on network size removed (2022/10/6)
% - Option for SBAS network plotting added (2022/10/6)

%% Initialization
clc; clear; close all;
format long g;
fclose('all');

% -----------------------------
% User-defined parameters
% -----------------------------
frame_ID = '020D_06347_151515';    % Frame ID for the LiCSAR data
patches = 40;                     % Number of patches in SBAS folder
process_name = 'ION_10';          % Name for the process folder
foldername = 'velocity';          % Subfolder to save outputs
refpoint_id = 10000;              % Reference point ID

% Plate motion model coefficients
lat_coef = -0.509070127346254;
lon_coef = -2.3066985484587;
cons_coef = 148.161299774462; 

%% Path Setup
Data_Path = strcat('../LiCSAR/', frame_ID); % Path to LiCSAR data
cd(Data_Path);
cd metadata;
load baselines;  % Load baseline data
cd ..;
mkdir(process_name);
cd(process_name);

% Create output folder if it doesn't exist
if exist(foldername, 'dir') ~= 7
    [~, msg, ~] = mkdir(foldername);
    disp(msg);
end

%% Load SBAS Data
cd ../SBAS;
load('temp_base.mat');    % Temporal baselines of SBAS network
load('master_date.mat');  % Dates of master images
load('slave_date.mat');   % Dates of slave images
load('SBAS2SM.mat');      % SBAS to single-master inversion matrix
load('sm_perp_base.mat'); % Perpendicular baselines
load('sm_temp_base.mat'); % Temporal baselines for single-master

% Wavelength of Sentinel-1 SAR imagery (in meters)
lambda = 5.5465763e-2;
coherence_threshold = 0;  % Coherence threshold for filtering
net_size = floor(1.5 * (size(A, 2) + 1)); % Filtering based on network size
N = length(sm_temp_base);

%% Model Definitions
% Model 1: Linear Trend
A_sm_1 = zeros(N, 2);
A_sm_1(:, 1) = -(4 * pi / lambda) * (sm_temp_base / 365); % Linear trend (yearly)
A_sm_1(:, 2) = -(4 * pi / lambda);                       % Constant term

% Model 2: Trend + Topographic Residual + Ionospheric Effects
A_sm_2 = zeros(N, 3);
A_sm_2(:, 1) = -(4 * pi / lambda) * (sm_temp_base / 365); % Linear trend
A_sm_2(:, 3) = -(4 * pi / lambda);                       % Constant term

% Model 3: Trend + Topographic Residual + Harmonic Signal
A_sm_3 = zeros(N, 4);
A_sm_3(:, 1) = -(4 * pi / lambda) * (sm_temp_base / 365); % Linear trend
A_sm_3(:, 3) = -(4 * pi / lambda);                       % Constant term
A_sm_3(:, 4) = -(4 * pi / lambda) * cos((2 * pi / 11) * ((1 / 365) * (sm_temp_base + 198))); % Harmonic signal

%% Process Each Patch
for patch_number = 1:patches
    tic;
    
    % Load patch-specific data
    load(strcat('patch', num2str(patch_number), 'unw.mat')); % Unwrapped phase
    load(strcat('patch', num2str(patch_number), 'Rsin_inc.mat')); % Incident angle matrix
    load(strcat('patch', num2str(patch_number), 'coh.mat')); % Coherence data
    load(strcat('patch', num2str(patch_number), 'coordinate.mat')); % Coordinate matrix

    fprintf('Processing Patch %i of %i\n', patch_number, patches);

    % Coherence filtering
    c = coherence_filtering(coh, coherence_threshold);
    unw_phase(c) = nan; % Mark low-coherence pixels as NaN

    % Initialize output arrays
    pixels = size(patch_coordinate, 1); % Number of valid pixels after filtering
    patch_vel = nan(pixels, 6); % Velocity vector initialization
    patch_ts = nan(pixels, size(B, 2)); % Time-series vector initialization
    
    for pixel_number = 1:pixels
        % Remove NaN observations
        L_sbas = unw_phase(pixel_number, :)';
        flag = sbas_plot(master_date, slave_date, L_sbas, baselines, flase);
        flag = flag(2:end);
        pix_coh = coh(pixel_number, ~isnan(L_sbas))';
        B = A(~isnan(L_sbas), :);
        B(:, ~flag) = [];
        Rsin = Rsin_inc(pixel_number, ~isnan(L_sbas))';
        L_sbas(isnan(L_sbas)) = [];

        % SBAS to Single-Master Conversion
        if rank(B) == size(B, 2) && ~isempty(L_sbas)
            A_sm_2(:, 2) = -(4 * pi / lambda) * sm_perp_base / Rsin(1);
            A_sm_3(:, 2) = -(4 * pi / lambda) * sm_perp_base / Rsin(1);
            newA_1 = A_sm_1(flag, :);
            newA_2 = A_sm_2(flag, :);
            newA_3 = A_sm_3(flag, :);

            % Time-series estimation
            ts_cap_L2 = (B' * B) \ (B' * L_sbas);
            [ts_cap_L1, W1] = IRLS(B, ts_cap_L2, L_sbas);

            % Residual estimation
            r_cap_L1 = (lambda / (4 * pi)) * (B * ts_cap_L1 - L_sbas);
            patch_ts(pixel_number, flag) = ts_cap_L1;

            % Parameter estimation
            x_cap_1 = (newA_1' * newA_1) \ (newA_1' * ts_cap_L1);
            x_cap_2 = (newA_2' * newA_2) \ (newA_2' * ts_cap_L1);
            x_cap_3 = (newA_3' * newA_3) \ (newA_3' * ts_cap_L1);
        end
    end
    toc;
end


function c = coherence_filtering(coh,t)
c = coh/256 < t;
end

function f = netfiltering(unw_phase,t)
f = false(size(unw_phase,1),1);
for ii = 1:size(unw_phase,1)
    not_nan = ~isnan(unw_phase(ii,:));
    n = sum(not_nan);
    f(ii,1) = (n < t);
end
end

function [m1,W]  = IRLS(G,m2,d)
r_cap = abs(G*m2 - d);
r_cap(r_cap < 1e-6) = 1e-6;
W = diag(1./r_cap);
m = zeros(size(m2));
m1 = m2;
while norm(m1 - m) / (1+norm(m1)) > 1e-4
    m = m1;
    m1 = (G'*W*G)\(G'*W*d);
    r_cap = abs(G*m2 - d);
    r_cap(r_cap < 1e-6) = 1e-6;
    W = diag(1./r_cap);
end
end


function dates_flag = sbas_plot(mdate,sdate,sbas,baselines,plot_flag)
close all
dates = sort(unique([mdate;sdate]));
I = [];
for ii = 1:length(dates)
     I = [I;find(baselines(:,2) == dates(ii))];
end
baselines = baselines(I,:);
baseline = baselines(:,3);
myear = floor(mdate/10000);
syear = floor(sdate/10000);
mmonth = floor((mdate/1e4 - myear)*1e2);
smonth = floor((sdate/1e4 - syear)*1e2);
mday = round(((mdate/1e4 - myear)*1e2 - mmonth)*1e2);
sday = round(((sdate/1e4 - syear)*1e2 - smonth)*1e2);
dates = baselines(:,2);
year = floor(dates/10000);
month = floor((dates/1e4 - year)*1e2);
day = round(((dates/1e4 - year)*1e2 - month)*1e2);
x = year + month/12 + day/365;
xm = myear + mmonth/12 + mday/365;
xs = syear + smonth/12 + sday/365;
% 
% 
dates_flag = false(size(dates));
% sdate_flag = false(size(sdate));
if plot_flag
    figure('WindowState','maximized')
    plot(x,baseline,'.r','MarkerSize',5)
    hold on
end
gmt = [];
for ii = 1:length(mdate)
    if ~isnan(sbas(ii))
        jj = find(mdate(ii,1)==dates);
        kk = find(sdate(ii,1)==dates);
        if plot_flag
            plot([xm(ii),xs(ii)],[baseline(jj),baseline(kk)],'-','Color',[0.8 0.8 0.8],'LineWidth',0.05)
            gmt = [gmt;xm(ii) baseline(jj);xs(ii) baseline(kk)];
        end
        dates_flag(jj) = true;
        dates_flag(kk) = true;
    end
end

if plot_flag
  plot(x(~dates_flag),baseline(~dates_flag),'.b')
  daspect([1 200 1])
  xlabel('Acquisition Date')
  ylabel('Perp. Baseline(m)')
end

[~,ix] = sort(x);
dates_flag = dates_flag(ix);



end

function I = crop_patch(B,C)
y  = C(:,1);
x = C(:,2);
x1 = B(1,1);
x2 = B(1,2);
y1 = B(2,1);
y2 = B(2,2);

I = find( x <= x2 & x >= x1 & y <= y2 & y >= y1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = winsinc(x,dt,fc,ns,win,ftype,fo)
    switch nargin
        case 0
            n = 500;
            t = linspace(0,20,n);
            dt = diff(t(1:2));
            s1 = 3*sin(2*pi*t/1);
            s2 = 6*sin(2*pi*t/3);
            s3 = 4*sin(2*pi*t*3);
            x = s1 + s2 + s3 + 1*randn(1,n);
            fc = 1;
            ns = 80;
            win = 'welch';
            win = 'kaiser';
            win = 'lanczos';
            %         win = 'gauss';
            %         win = 'boxcar';
            ftype = 'low';
            %         fc = 0.09;
            %         ftype = 'pass';
            %         fo = 1;
            winsinc(x,dt,fc,ns,win,ftype);
            return
        case {1,2,3}
            error('Faltan argumentos de entrada')
        case 4
            ftype = 'low'; %pasa baja
        case {5,6}
            ftype = lower(ftype);
            switch ftype
                case {'low','high'}
                    if nargin==6
                        warning('El filtrado pasa-baja y alta no requieren fo')
                    end %if
                case 'pass'
                    if nargin==5
                        error('Bandpass filter requires fo')
                    end %if
                otherwise
                    error('Invalid Filter Type');
            end %switch
        otherwise
            error('Too many Inputs')
    end %switch
    
    N = length(x);
    fn = (1/(2*dt)); %Nyquist frequency
    df = 2*fn/(N-1); %delta frec
    fcc = fc;
    %Frecuencia de corte para pasa alta
    if strcmpi(ftype,'high')
        fc = fn - fc;
    end
    delta = 4*fn/(2*ns+1);
    d = delta/2;
    if d>fc
        warning('Cutoff Frequency too small')
        fc = d;
    end
    n = -ns:ns;
    h = zeros(1,2*ns+1);
    hh = h;
    warning off MATLAB:divideByZero
    w = zeros(size(n));
    an = abs(n);
    ic = (an<ns);
    switch lower(win(1:3))
        case 'wel' %welch
            w(ic) = 1 - (n(ic)/ns).^2
            w(~ic) = 0;
        case 'par' %parzen
            ic1 = (an>=0 & an<1);
            w(ic1) = 1/4*(4 - 6*an(ic1).^2  + 3*an(ic1).^3);
            ic2 = (an>=1 & an<2);
            w(ic2) = 1/4*(2 - an(ic2)).^3
            ic3 = (~ic1 & ~ic2); %otherwise
            w(ic3) = 0;
        case 'han' %hanning
            alpha = 1/2;
            w(ic) = alpha + (1 - alpha)*cos(pi*n(ic)/ns);
            w(~ic) = 0;
        case 'ham' %hamming
            alpha = 0.54;
            w(ic) = alpha + (1 - alpha)*cos(pi*n(ic)/ns);
            w(~ic) = 0;
        case 'bla' %blackman
            w(ic) = 0.42 + 0.5*cos(pi*n(ic)/ns) ...
                + 0.08*cos(2*pi*n(ic)/ns);
            w(~ic) = 0;
        case 'lan' %lanczos
            w(ic) = sin(pi*n(ic)/ns)./(pi*n(ic)/ns);
            w(~ic) = 0;
        case 'kai' %kaiser
            alpha = 0.01;
            y = alpha*sqrt(1 - (n(ic)/ns).^2);
            w(ic) = besselj(0,y)./besselj(0,alpha);
            w(~ic) = 0;
        case 'gau' %gauss
            sigma = 30;
            w(ic) = exp(-(n(ic)/sigma).^2);
            w(~ic) = 0;
        case 'box' %boxcar
            w = ones(size(n));
        otherwise
            error('Invalid Window')
            %w = ones(ns,1)/ns;
    end %switch
    %Basic Filter (sinc)
    a = fc/fn; %0<a<1
    c = sin(pi*n*a)./(pi*n);
    %Obtain filter Coefficients
    %Multiplication is time domain
    %is a convolution in frequency domain
    h = w.*c;
    ind = find(n==0);
    h(ind) = a;
    %pasa alta
    if strcmpi(ftype,'high')
        h(2:2:length(h)) = -1*h(2:2:length(h));
    end
    hh = h;
    hh(ind) = 0;
    %Filter Frequencies
    Fh = 0:df:fn;
    
    switch ftype
        case {'low','high'}
            %Respuesta en Frecuencia del Filtro
            for k=1:length(Fh);
                Rh(k) = a + sum(hh.*cos(2*pi*n*Fh(k)/(2*fn)));
            end
            %Filtrado por convolucion entre la Serie de Datos
            %y los Coeficientes del Filtro
            xf = conv(x,h);
            xf = wkeep(xf,N);
        case 'pass'
            %Respuesta en Frecuencia del Filtro
            for k=1:length(Fh);
                Rh(k) = a + sum(hh.*cos(2*pi*n*(Fh(k)-fo)/(2*fn)));
            end
            %Corrimiento de los pesos del filtro
            th = n*dt;
            i = sqrt(-1); %Redundant but instructive
            hpri = h.*exp(2.0*i*pi*fo.*th);
            %Filtrado en tiempo: Convolucion entre la Serie
            %de Datos y los Coeficientes del Filtro
            xfc = conv(x,hpri);
            xf = 2.0*real(xfc);
            xf = wkeep(xf,N);
    end
    xf = xf(:);
    if nargout==0
        subplot(3,1,1)
        plot(x,'b')
        hold on
        plot(xf,'r')
        hold off
        legend('Original','Filtered')
    
        subplot(3,2,3)
        plot(h)
        %     hold on
        %     plot(hh,'r')
        ylabel('Filter Coef.')
        set(gca,'xlim',[0,length(h)])
    
        subplot(3,2,4)
        [psd,f] = spectral(x,dt,'hanning');
        loglog(f,psd)
        hold on
        [psd,f] = spectral(xf,dt,'hanning');
        loglog(f,psd,'r')
        set(gca,'xlim',[f(1),f(end)])
        set(gca,'xticklabel',get(gca,'xtick'))
        hold off
        ylabel('PSD(f)')
    
        subplot(3,1,3)
        semilogx(Fh,Rh)
        hold on
        lim = axis;
        plot([fcc,fcc],[lim(3),lim(4)],'g')
        if exist('fo') & fo>Fh(1)
            plot([fo,fo],[lim(3),lim(4)],'r')
        end
        set(gca,'xlim',[Fh(1),Fh(end)],'ylim',[-0.1,1.1])
        set(gca,'xticklabel',get(gca,'xtick'))
        hold off
        xlabel('Frequency')
        %     xlabel('Period')
        ylabel('Response')
    
    end
    if nargout>0
        varargout{1} = xf;
    end
    if nargout>1
        varargout{2} = Rh;
    end
    if nargout>2
        varargout{3} = Fh;
    end
    warning on MATLAB:divideByZero
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Decomposing InSAR time_series using FIR filter  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ts_low,ts_high] = time_series_decomp(ts_res,sm_temp_base)
 
    I = find(~isnan(ts_res));
    ts = ts_res(I)';
    tb = sm_temp_base(I);
    t = [-tb(end:-1:1)',0,tb'];
    ns = length(t);
    win = 'lanczos';
    ftype = 'low';
    [ts_low,~,~] = winsinc(ts,mean(diff(tb)),1/(1.5*365),ns,win,ftype);
    ts_high = ts - ts_low;
    
end
