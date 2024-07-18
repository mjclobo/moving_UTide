function [coeffs,CD] = moving_UTide(DATA,window,incr,good_pct,rayleigh,t1_orig,dt_orig)
% DATA: data structure containing time in datenumber format (datenums) and heights (dataHi)
% window: window for the harmonic analysis (e.g.:  32 days)
% incr: number of days to increment forward with each HA loop
% plot_flag: plot if ==1
% good_pct: what percent of data in window must be non-NaN
% rayleigh: Rayleigh criterion, implies cnstit arg is auto

% datenums = datenum(DATA.dates);
datenums = DATA.datenums;
% dtimes   = convertTo(DATA.dtime,'epochtime','Epoch','1970-01-01');
dtimes   = convertTo(DATA.dtime,'epochtime','Epoch',DATA.dtime(1));


dt = diff(datenums);
dt = median(dt,'omitnan')*24*3600;

dataHi = DATA.wl;

%%
num = 1;
tb = datenums(1) + window/2;

memBool=1;

possible_numdata  = (window*24*3600)/dt;
timeSteps = tb:incr:datenums(end)-window/2;
num_runs = length(timeSteps);

ref_times = zeros(num_runs,1);

CD.M4_M2sqr = zeros(num_runs,1);
t_HA = zeros(num_runs,1);

coeffs = cell(1,num_runs);

for num=1:num_runs
    
    tb = timeSteps(num);
    
    if rem(num,30)==1
        disp(['Loop ',num2str(num)]);
    end

    dk = find(datenums>=(tb -window/2) & datenums < (tb+window/2));
    kk = find(isfinite(dataHi(dk)));  % only use non-nan height data in utide
    jj = find(isfinite(datenums(dk(kk)))); % only use non-nan time data in utide
    
    if length(kk(jj))> (good_pct*possible_numdata)  % require a minimum amount of data before doing processing loop
        % the combination of this threshold and window and increment
        % affects which constits can be resolved
        
        % Utide solver
%         coef = ut_solv(datenums(dk(kk(jj))), dataHi(dk(kk(jj))),[],DATA.lat,'auto','GwchNone','NodsatNone','RunTimeDisp',...
%             'nnn','Rmin',rayleigh,'White','OLS','OrderCnstit','frq');
        
         coef = ut_solv(datenums(dk), dataHi(dk),[],DATA.lat,'auto','GwchNone','NodsatNone','RunTimeDisp',...
            'nnn','Rmin',rayleigh,'White','OLS','OrderCnstit','frq');
        if memBool==1
            for i=1:length(coef.name)
                nombre = invalidNameIn(coef.name{i});
                CD.(nombre) = allocateMem(num_runs);
            end
            memBool=0;
        end
        
        if num==1
            names_all = coef.name; % YOUR FIRST WINDOW MUST NOT HAVE ANY MISSING VALUES
        end
        
        for i=1:length(names_all)
            nombre = invalidNameIn(names_all{i});
            CD = assignValues(CD,coef,nombre,num);
        end
        
%         CD.M4_M2sqr(num) = CD.M4.amp(num)./(CD.M2.amp(num)*CD.M2.amp(num));
        % adjust phase to be relative to whole time series' start time
        t = dtimes(dk(kk(jj)));
        coef.aux.start_time = tb;  % start time of this analysis; using original time axis
        
        ref_times(num) = coef.aux.reftime;
        coef.time_alt = mean(t); %(t(1)+t(end))/2;
        
        coef.aux.t1_orig = t1_orig;         % initial time of the time series; on original time axis
        
        % save ut_solv() complete output
        coeffs{num} = coef;
%         for i=1:length(names_all)
%             bad_name = revertToBadName(names_all{i});
%             aa2 = startsWith(coef.name,bad_name);
%             bb2 = find(aa2==1);
%             if isempty(bb2)
%                 coeffs{num}
%             end
%         end

    else  % if not enough data points.....
        
        if memBool==0
            for i=1:length(coef.name)
                nombre = invalidNameIn(coef.name{i});
                CD = assignNans(CD,nombre,num);
            end
%             CD.M4_M2sqr(num)=nan;
            coeffs{num} = [];
        end
    
    end

t_HA(num) = tb;

end



CD.datetimes = datetime(t_HA,'ConvertFrom','datenum');
CD.dates = t_HA;

% these two times should be same thing..
CD.times_analyzed = timeSteps;
CD.ref_times = ref_times; 

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
    structOut.snr = zeros(num_runs,1);
end

function structIn = assignValues(structIn,coeffIn,name,num)
    % build structure for every constituent that UTide puts out
    badName = revertToBadName(name);
    aa = startsWith(coeffIn.name,badName);
    bb = find(aa==1);
    if isempty(bb)
        structIn.(name).amp(num) = 0.;
        structIn.(name).ampci(num) = 0.; structIn.(name).phase(num) = 0.; 
        structIn.(name).phaseci(num) = 0.;
        structIn.(name).snr(num) = NaN;
    else
        structIn.(name).amp(num) = coeffIn.A(bb);
        structIn.(name).ampci(num) = coeffIn.A_ci(bb); structIn.(name).phase(num) = coeffIn.g(bb); 
        structIn.(name).phaseci(num) = coeffIn.g_ci(bb);
        structIn.(name).snr(num) = coeffIn.diagn.SNR(bb);
    end
    
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

