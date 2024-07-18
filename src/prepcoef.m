function [results,OUT,cdnew] = prepcoef(coef,t_in,t_out)
% Whole goal here is to build new coef variable that includes
% amp/phase info for interpolated constituents.
% 
% Reconstruction is also done here now, as well.
%
% *** INPUTS ***
%
% coef: output from ut_solv
% t_in: times associated with moving harmonic analyses
% t_out: desired times for reconstructed tide data; t_in and t_out cannot
% be vastly different!
% t_or: the time series of the original data, it is used for the
% parameter OUT.decwl (see below)
% 
% *** OUTPUTS ***
% 
% results: shows intermediate values for amps, phase, complex values, etc. 
% OUT: contains all data pertaining to reconstructed water levels
% cdnew: contains all data needed for reconstruction, i.e. amplitudes and
% phases corresponding to t_out, as well as fields ut_reconstr needs
%
% see end of script for details of OUT

N = length(coef);
n = length(coef{1}.name);
ampsin = zeros(N,n);
phasesin = zeros(N,n);
rel_phase = zeros(N,n);
rel_phase_alt = zeros(N,n);
rel_time  = zeros(N,n);
int0 = zeros(N,n);
intout = zeros(n,length(t_out));
reftimein = zeros(N,1);
meanin = zeros(N,1);
slopein = zeros(N,1);

lat = coef{1}.aux.lat;
lind = coef{1}.aux.lind;
names = coef{1}.name;

for k=1:N
    ch = coef{k};
    if ~isempty(ch)
        for j=1:n
            ampsin(k,j) = ch.A(j); 
            phasesin(k,j) = ch.g(j);
            rel_time(k,j) = ch.aux.start_time;
            
            comp = ampsin(k,j)*exp(-1i*ch.g(j)*(pi/180));
            
            % ch.aux.reftime is average of first and last input times for each
            % window analysis
            
            rel_rad = 2 * pi *(ch.aux.reftime*24*3600) * (ch.aux.frq(j)/3600.0) ; % SAME AS ONE BELOW
            % rel_rad = 2 * pi *(t_in(k)) * (ch.aux.frq(j)/3600.0) ; SAME AS ONE ABOVE
            rel_phase(k,j) = mod(rel_rad - angle(comp),2*pi);
            
            rel_rad_alt = 2 * pi *(ch.time_alt) * (ch.aux.frq(j)/3600.0) ; %
            rel_phase_alt(k,j) = mod(rel_rad_alt - angle(comp),2*pi);            
            
            int0(k,j) = ch.A(j)*exp(-1i*rel_phase(k,j));
        end
        reftimein(k) = ch.aux.reftime;
        meanin(k) = ch.mean;
        slopein(k) = ch.slope;
    else 
    
    end
end

% for k=1:n   % for each constituent
%     ii = find(int0(:,k));
%     intout(k,:) = interp1(t_in(ii),int0(ii,k),t_out,'linear'); % bad; times align...
% 
% end

%% Alternate way  to interpolate complex solution, taken from CWT code
for k=1:n
    ii = find(int0(:,k));
    solnFunc{k} = spline(t_in(ii),int0(ii,k));
end

nData = length(t_out);

reconMat = zeros(n,nData);

for j=1:nData
    for k=1:n
        omega = coef{1}.aux.frq(k)*(2*pi)/3600;
        reconMat(k,j)= real(ppvalFast(solnFunc{k},t_out(j)) * exp(1i* omega * t_out(j)));
    end
end

ampLimit = 1000;
ampFloor = 1e-10;

% reconWts = cellit(@(it2,it3) (abs(reconMat(it2,it3)) < ampLimit) * 1, 1:nRecon,1:nData);

reconHi = cell2mat(cellit(@(it2) sum(reconMat(:,it2)),1:nData));

% reconAll = reconHi + dataInLo.' + meanI;

%% continuing old code

ii = find(reftimein);
% refout = interp1(t_in(ii),reftimein(ii),t_out);
% meanout = interp1(t_in(ii),meanin(ii),t_out,'linear');
% slopeout = interp1(t_in(ii),slopein(ii),t_out,'linear');
% 
% % now let's reconstruct
% ampsout = zeros(length(t_out),n);
% phasesout = zeros(length(t_out),n);


cdnew = cell(length(t_out),1);

% for k=1:length(t_out)
%     for j=1:n
%         im = intout(j,k);
%         ampsout(k,j) = abs(im);
%         phasesout(k,j) = mod(atan2d(imag(im),real(im)),360);
%     end
%     
%     % fields needed for ut_reconstr
%     cdnew{k}.A = ampsout(k,:)';
%     cdnew{k}.g = phasesout(k,:)';
%     cdnew{k}.aux.frq = coef{1}.aux.frq;
%     cdnew{k}.aux.reftime = refout(k);
%     cdnew{k}.mean = meanout(k);
%     cdnew{k}.aux.lat = lat;
%     cdnew{k}.aux.lind = lind;
%     cdnew{k}.slope = slopeout(k);
%     cdnew{k}.aux.opt.notrend = 0;
%     cdnew{k}.aux.opt.prefilt = [];
%     cdnew{k}.aux.opt.nodsatlint = 0;
%     cdnew{k}.aux.opt.nodsatnone = 1;
%     cdnew{k}.aux.opt.gwchlint = 0;
%     cdnew{k}.aux.opt.gwchnone = 1;
%     cdnew{k}.aux.opt.twodim=0;
%     cdnew{k}.name = names;
%  %   cdnew{k}.aux.opt.notrend = 1;
%     % these are dummy variables for now
%     cdnew{k}.A_ci = 0.001 * ones(1,n)';
%     cdnew{k}.g_ci = 0.001 * ones(1,n)';
% end

% 
% sl = zeros(length(t_out),1);
% 
% % think about band-limited interpolation here, this needs to be at dt=1hr
% for k = 1:length(t_out)
%     sl(k) = ut_reconstr(t_out(k),cdnew{k});
% end

%% once we have reconstruction...

% OUT.wl = sl ;           % water level reconstructed from moving_HA params and t_out
OUT.wl = reconHi;
OUT.dates = t_out';     % times corresponding to reconstructed data
OUT.lat = lat;
% OUT.amps = ampsout;
OUT.phase = phasesin;
OUT.names = names;
OUT.rel_phase = rel_phase*(180/pi);
OUT.rel_time = rel_time;
OUT.rel_phase_alt = rel_phase_alt*(180/pi);
% OUT.decwl = interp1(t_out,OUT.wl,t_or);
% OUT.decmean = interp1(t_out,meanout,t_or,'linear');
% OUT.mean = meanout';    % subtract this from OUT.wl to "high-pass" your data

results.int = intout;
results.int0 = int0;
results.amps = ampsin;
results.mean = meanin;
results.phases = phasesin;

for k=1:length(names)
    nombre = invalidNameIn(names{k});
    
    OUT.(nombre).phases = rel_phase(:,k)*(180/pi);
    OUT.(nombre).amps   = ampsin(:,k);
    
    OUT.(nombre).phases_alt = rel_phase_alt(:,k)*(180/pi);
end



end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function structOut = allocateMem(num_runs)
    % allocate memory for all UTide analyses
    structOut.amp = zeros(num_runs,1);
    structOut.ampci = zeros(num_runs,1);
    structOut.phase = zeros(num_runs,1);
    structOut.phaseci = zeros(num_runs,1);
    structOut.sig = zeros(num_runs,1);
end

function structIn = assignValues(structIn,coeffIn,name,num)
    % build structure for every constituent that UTide puts out
    badName = revertToBadName(name);
    aa = startsWith(coeffIn.name,badName);
    bb = find(aa==1); structIn.(name).amp(num) = coeffIn.A(bb);
    structIn.(name).ampci(num) = coeffIn.A_ci(bb); structIn.(name).phase(num) = coeffIn.g(bb); 
    structIn.(name).phaseci(num) = coeffIn.g_ci(bb);
    structIn.(name).sig(num) = coeffIn.diagn.SNR(bb);
end

function structIn = assignNans(structIn,name,num)
     structIn.(name).amp(num) = nan; structIn.(name).phase(num) = nan;  % not enough good data, set structIn.(name) measred to nan
     structIn.(name).ampci(num) = nan; structIn.(name).phaseci(num) = nan; structIn.(name).sig(num) = nan;
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

function badName = revertToBadName(nameIn)
    % only goes up to name beginning with 'nine'
    if all(nameIn(1:2)=='On')
        badName = strcat('1',nameIn(4:end));
    elseif all(nameIn(1:2)=='Tw')
        badName = strcat('2',nameIn(4:end));
    elseif all(nameIn(1:2)=='Th')
        badName = strcat('3',nameIn(6:end));
    elseif all(nameIn(1:2)=='Fo')
        badName = strcat('4',nameIn(5:end));
    elseif all(nameIn(1:2)=='Fi')
        badName = strcat('5',nameIn(5:end));
    elseif all(nameIn(1:2)=='Si')
        badName = strcat('6',nameIn(4:end));
    elseif all(nameIn(1:2)=='Se')
        badName = strcat('7',nameIn(6:end));
    elseif all(nameIn(1:2)=='Ei')
        badName = strcat('8',nameIn(6:end));
    elseif all(nameIn(1:2)=='Ni')
        badName = strcat('9',nameIn(5:end));
    else
        badName = nameIn;
    end
end



