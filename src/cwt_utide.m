function [cd,OUT] = cwt_utide(dataIn,window,incr,good_pct,rayleigh,maxFLength,statWindow)

dt = seconds(nanmedian(diff(dataIn.dates)));

t1_orig = dataIn.datenums(1)*24*3600; % + dt * (maxFLength-1) / 2;

decFact = round(incr*24/(dt/3600));

nDataIn=length(dataIn.wl);
nDataOdd = nDataIn-mod(nDataIn+1,2);       % modified this from t_tide

tPts = dt*([1:nDataOdd]-ceil(nDataOdd/2));

times_alt = dataIn.dtime;

dt_orig = tPts(1)-t1_orig;

timesAnalyzed = tPts(floor((maxFLength-1)/2):(incr*24):nDataIn-floor((maxFLength-1)/2));      % times without edges

% t1_orig = timesAnalyzed(1);


if mod(length(timesAnalyzed),2) == 1
    odd = 1;
else
    timesAnalyzed = timesAnalyzed(1:end-1);
    odd = 0;
end

nDec = length(timesAnalyzed);

%% filtering data in same way as CWT analysis does
% defining filter criteria

co_stop = (maxFLength*3600)^-1;
co_pass = ((maxFLength+2*24)*3600)^-1;     % i.e., pass band of 1 cyc/day

dpass_pct_lo = 0.1; % 5 % passband ripple
dstop_db_lo = 0.005; % 47db stopband attenuation

% defining optimal filters
[co_n,wn,beta,ftype] = kaiserord([co_pass,co_stop],[1,0],[dpass_pct_lo,dstop_db_lo],dt^-1);
co_filt_lo = fir1(co_n,wn,ftype,kaiser(co_n+1,beta),'noscale');

% take care of NaNs
nanI = find(isnan(dataIn.wl));
dataIn.wl(nanI) = nanmean(dataIn.wl);

% applying filters
dataLo = conv(dataIn.wl,co_filt_lo,'same');
dataHi = dataIn.wl - dataLo;
dataIn.wl = dataHi;

dataIn.wl(nanI) = NaN;
% 
% p=figure();
% plot(dataIn.dates,dataIn.wl)

%% doing moving window harmonic analysis

[coeffs,cd] = moving_UTide(dataIn,window,incr,good_pct,rayleigh,t1_orig,dt_orig);

cd.wlHi = dataIn.wl;

% use coeffs{1}.name to see which constits we have...and their indices!!

%% reconstruction


% dt = 1/24;      % dt for reconstruction in days
t_int = dataIn.datenums;
% t_int = tPts;


% note that both t_int(i.e. cdnew.dates) and cd.dates are in the form of
% "hours since a very far away date"

% [res,OUT,cdnew] = prepcoef(coeffs,cd.dates,t_int,t_int);    % here's where recon happens
[res,OUT,cdnew] = prepcoef(coeffs,cd.times_analyzed*24*3600,t_int*24*3600);    % here's where recon happens


t_in = cd.dates;
t_out = t_int;
a = interp1(t_in,res.int0(:,1),t_out,'pchip');

OUT.wlAll = OUT.wl+dataLo.';

OUT.dtimes = dataIn.dates;

% statistics
[~,statStart] = min(abs(statWindow(1)-dataIn.dates));
[~,statEnd]   = min(abs(statWindow(2)-dataIn.dates));

[~,statStartDec] = min(abs(statWindow(1)-dataIn.dates(1:decFact:end)));
[~,statEndDec]   = min(abs(statWindow(2)-dataIn.dates(1:decFact:end)));

if statEndDec > length(OUT.(OUT.names{1}).phases)
    statEndDec = length(OUT.(OUT.names{1}).phases);
end

for i=1:length(OUT.names)
    nombre = invalidNameIn(OUT.names{i});
    cd.(nombre).avgSNR = nanmean(cd.(nombre).snr(statStartDec:statEndDec));
end

OUT.dtimesTrim = OUT.dtimes(statStart:statEnd);

OUT.resid = dataIn.wl(statStart:statEnd) - OUT.wl(statStart:statEnd).';

OUT.bias = mean(OUT.resid(~isnan(OUT.resid)));

OUT.rmse = sqrt(mean(OUT.resid(~isnan(OUT.resid)).^2));

% adding power spectrum part, I guess :,)

c = OUT.resid;

bad_vals = find(isnan(c));
for i=1:10
    c(bad_vals) = 0;
    c = c-mean(c);
end

d = floor(length(c)/2);         % so we get same-sized spectrum as residuals below

nfft = 2^12;

[pxx,f] = pwelch(c,hanning(d),round(d/2),nfft,1/dt);

OUT.resid_pwr=pxx/3600/24;
OUT.resid_freqs=f*3600*24;      % frequency in cpd


end

function nameOut = invalidNameIn(nameIn)
    % only goes up to name beginning with '9'
    if nameIn(1)=='1'
        nameOut = strcat('One',nameIn(2:end));
    elseif nameIn(1)=='2'
        nameOut = strcat('Two',nameIn(2:end));
    elseif nameIn(1)=='3'
        nameOut = strcat('Three',nameIn(2:end));
    elseif nameIn(1)=='4'
        nameOut = strcat('Four',nameIn(2:end));
    elseif nameIn(1)=='5'
        nameOut = strcat('Five',nameIn(2:end));
    elseif nameIn(1)=='6'
        nameOut = strcat('Six',nameIn(2:end));
    elseif nameIn(1)=='7'
        nameOut = strcat('Seven',nameIn(2:end));
    elseif nameIn(1)=='8'
        nameOut = strcat('Eight',nameIn(2:end));
    elseif nameIn(1)=='9'
        nameOut = strcat('Nine',nameIn(2:end));
    else
        nameOut = nameIn;
    end
end