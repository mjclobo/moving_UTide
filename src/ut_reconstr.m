function [u,v] = ut_reconstr(tin,coef,varargin)
% UT_RECONSTR() Reconstruct superposed harmonics. 
%   Generates hindcast and/or forecast/prediction at user-specified times
%       using output from companion harmonic analysis function UT_SOLV().
%
% OVERVIEW 
% 
% Syntax for two-dimensional raw input, such as velocities:
%   [ u_fit, v_fit ] = UT_RECONSTR ( t_fit, coef , {options} );
% 
% Syntax for one-dimensional raw input, such as sea level:
%   [ sl_fit, ~ ] = UT_RECONSTR ( t_fit, coef , {options} ); 
% 
% Input/output parameters (more detailed explanations below): 
%   * t_fit = arbitrary times for reconstructed superposed harmonics
%   * coef = results output structure from UT_SOLV()
%   * u_fit & v_fit, or sl_fit = reconstructed superposed harmonics
%   * {options} explained below
%
% Analysis of groups of records
%   The descriptions that follow next for INPUTS, DEFAULTS, OPTIONS, and 
%   OUTPUTS are for treatment of a single record. Following that, 
%   explanations are given for modifications that enable treating a group 
%   of records with a single execution.
%
% INPUTS
% 
% t_fit 
%
%   * Column vector of arbitrary times, datenum UCT/GMT, distributed
%       uniformly or irregularly.
%   * May include NaNs and if so the outputs (u_fit/v_fit or sl_fit) will 
%       have corresponding NaNs. 
%
% coef
%
%   * Output generated by call to UT_SOLV().
%
% {options}
%
%   * See section below.
%
% DEFAULTS
%
%   * Constituents with SNR > 2 are included in the superposition. 
%   * All other aspects of the calculation are determined by the 
%       information in coef, such that they match the calculation done 
%       by UT_SOLV(); to change them another run of UT_SOLV() is needed.
%
% OPTIONS
%
% Option flags are not case-sensitive but cannot be abbreviated. The order
%   of the option flags is not important but they must be passed in 
%   after all other arguments. See report for more complete explanations.
%
%   �MinSNR�, MinSNR 
%
%       * Only include constituents with SNR >= MinSNR. Default MinSNR=2.
%
%   �MinPE�, MinPE 
%
%       * Only include constituents with PE > MinPE. Default MinPE=0. 
%
% If both of �MinSNR� and �MinPE� are selected no constituent with either 
%   SNR or PE values lower than the specified thresholds will be included.
%
%   �Cnstit�, Cnstit
%
%       * Include only those constituents named in Cnstit, which must be 
%           selected from those in coef.name.
%               Cnstit = cell array of 4-character strings.
%       * If �Cnstit� is used then MinSNR and MinPE are ignored.
% 
% OUTPUTS
%
% u_fit & v_fit, or sl_fit
%
%   * Column vector(s) (same size as t_fit) containing the superposed
%       harmonics and the mean (and trend if included in model).
%
% GROUPS OF RECORDS
%
% Inputs are as in the single-record case except that t_fit can be either
%   a single n_t x 1 vector of times, to be used for all time sequences
%   in the group, or an n_t x n1 x n2 x n3 ... n-n_d array of times, which
%   specifies a different set of n_t times for each record. See comments
%   for UT_SOLV() for explanations of the parameters n1, n2, ..., n-n_d.
%
% Output is the same as in the single-record case except that u_fit & 
%   v_fit, or sl_fit, are n_t x n1 x n2 x n3 ... n-n_d arrays (the same 
%   size as the u_raw & v_raw, or sl_raw, inputs to the UT_SOLV() run
%   that created coef). 
%
% For more information see:
%   Codiga, D.L., 2011. Unified Tidal Analysis and Prediction Using the 
%       UTide Matlab Functions. Technical Report 2011-01. Graduate School 
%       of Oceanography, University of Rhode Island, Narragansett, RI. 
%       59pp. ftp://www.po.gso.uri.edu/pub/downloads/codiga/pubs/
%       2011Codiga-UTide-Report.pdf
%
% UTide v1p0 9/2011 d.codiga@gso.uri.edu

if isequal(size(coef.g,2),1) % single record
    [u,v] = ut_reconstr1(tin,coef,varargin{:});
else % group of records
    n_t = size(tin,1); 
    % alln = [n_1 n_2 n_3 ... n_n_d]
    alln = size(coef.g);
    alln = alln(2:end);
    % check inputs
    if (~isequal(size(tin),[n_t 1]) && ~isequal(size(tin),[n_t alln]))
        error(['ut_reconstr(group): input times have neither acceptable'...
            ' size.']);
    end
    % n_d = number of dimensions of group
    n_d = length(alln);
    % n_s = number of time sequences in group
    n_s = prod(alln); 
    % nallc = number of all constituents
    nallc = size(coef.g,1);
    % reshape t to a vector
    if size(tin,2)>1
        tin = reshape(tin,n_t,n_s);
    end
    % reshape fields of coef to vectors
    if n_d > 1
        coef.g = reshape(coef.g,nallc,n_s);
        coef.g_ci = reshape(coef.g_ci,nallc,n_s);
        if coef.aux.opt.twodim 
            coef.Lsmaj = reshape(coef.Lsmaj,nallc,n_s);
            coef.Lsmaj_ci = reshape(coef.Lsmaj_ci,nallc,n_s);
            coef.Lsmin = reshape(coef.Lsmin,nallc,n_s);
            coef.Lsmin_ci = reshape(coef.Lsmin_ci,nallc,n_s);
            coef.theta = reshape(coef.theta,nallc,n_s);
            coef.theta_ci = reshape(coef.theta_ci,nallc,n_s);
            coef.umean = reshape(coef.umean,1,n_s);
            coef.vmean = reshape(coef.vmean,1,n_s);
            if ~coef.aux.opt.notrend
                coef.uslope = reshape(coef.uslope,1,n_s);
                coef.vslope = reshape(coef.vslope,1,n_s);
            end
        else
            coef.A = reshape(coef.A,nallc,n_s);
            coef.A_ci = reshape(coef.A_ci,nallc,n_s);
            coef.mean = reshape(coef.mean,1,n_s);
            if ~coef.aux.opt.notrend
                coef.slope = reshape(coef.slope,1,n_s);
            end
        end
        coef.aux.reftime = reshape(coef.aux.reftime,1,n_s);
        if size(coef.aux.lat,2)>1
            coef.aux.lat = reshape(coef.aux.lat,1,n_s);
        end
    end
    % initialize storage
    u = nan*ones(n_t,n_s);
    if coef.aux.opt.twodim
        v = u;
    else
        v = [];
    end
    % main loop
    coef1.aux.frq = coef.aux.frq;
    coef1.aux.lind = coef.aux.lind;
    coef1.aux.opt = coef.aux.opt;
    for is = 1:n_s
        if size(tin,2) > 1
            tin1 = tin(:,is);
        else
            tin1 = tin;
        end
        % select coef results for the is'th record
        coef1.name = coef.name;
        coef1.g = coef.g(:,is);
        coef1.g_ci = coef.g_ci(:,is);
        if coef.aux.opt.twodim
            coef1.Lsmaj = coef.Lsmaj(:,is);
            coef1.Lsmaj_ci = coef.Lsmaj_ci(:,is);
            coef1.Lsmin = coef.Lsmin(:,is);
            coef1.Lsmin_ci = coef.Lsmin_ci(:,is);
            coef1.theta = coef.theta(:,is);
            coef1.theta_ci = coef.theta_ci(:,is);
            coef1.umean = coef.umean(is);
            coef1.vmean = coef.vmean(is);
            if ~coef.aux.opt.notrend
                coef1.uslope = coef.uslope(is);
                coef1.vslope = coef.vslope(is);
            end
        else
            coef1.A = coef.A(:,is);
            coef1.A_ci = coef.A_ci(:,is);
            coef1.mean = coef.mean(is);
            if ~coef.aux.opt.notrend
                coef1.slope = coef.slope(is);
            end
        end
        if size(coef.aux.lat,2)>1
            coef1.aux.lat = coef.aux.lat(is);
        else
            coef1.aux.lat = coef.aux.lat;
        end
        coef1.aux.reftime = coef.aux.reftime(is);
        % execute reconstruct for one record
        fprintf('[%d/%d]',is,n_s);
        [u1,v1] = ut_reconstr1(tin1,coef1,varargin{:});
        % store results
        u(:,is) = u1;
        if coef.aux.opt.twodim
            v(:,is) = v1;
        end
    end
    % reshape back to original array sizes
    if n_d > 1
        u = reshape(u,[n_t alln]);
        if coef.aux.opt.twodim
            v = reshape(v,[n_t alln]);
        end
    end
end

%%--------------------------------------------------------- 
function [u,v] = ut_reconstr1(tin,coef,varargin)
% UT_RECONSTR1()
% Reconstruction for a single record. See comments for UT_RECONSTR().
% UTide v1p0 9/2011 d.codiga@gso.uri.edu 

fprintf('ut_reconstr: ');

% parse inputs and options
[t,opt] = ut_rcninit(tin,varargin);

% determine constituents to include
if ~isempty(opt.cnstit)
    [~,ind] = ismember(cellstr(opt.cnstit),coef.name);
    if ~isequal(length(ind),length(cellstr(opt.cnstit)))
        error(['ut_reconstr: one or more of input constituents Cnstit '...
            'not found in coef.name']);
    end
else
    ind = 1:length(coef.aux.frq);
    if coef.aux.opt.twodim
        SNR = (coef.Lsmaj.^2 +coef.Lsmin.^2)./...
            ((coef.Lsmaj_ci/1.96).^2 + (coef.Lsmin_ci/1.96).^2);
        PE = sum(coef.Lsmaj.^2 + coef.Lsmin.^2);
        PE = 100*(coef.Lsmaj.^2 + coef.Lsmin.^2)/PE;
    else
        SNR = (coef.A.^2)./((coef.A_ci/1.96).^2);
        PE = 100*coef.A.^2/sum(coef.A.^2);
    end

    ind = ind(SNR(ind)>=opt.minsnr & PE(ind)>=opt.minpe);
    
end

% complex coefficients
rpd = pi/180;
if coef.aux.opt.twodim
    ap = 0.5*(coef.Lsmaj(ind) + coef.Lsmin(ind)) .* ...
        exp(1i*(coef.theta(ind) - coef.g(ind))*rpd);
    am = 0.5*(coef.Lsmaj(ind) - coef.Lsmin(ind)) .* ...
        exp(1i*(coef.theta(ind) + coef.g(ind))*rpd);
else
    ap = 0.5*coef.A(ind).*exp(-1i*coef.g(ind)*rpd);
    am = conj(ap);
end

% exponentials
ngflgs = [coef.aux.opt.nodsatlint coef.aux.opt.nodsatnone ...
    coef.aux.opt.gwchlint coef.aux.opt.gwchnone];
% fprintf('prep/calcs ... ');
E = ut_E(t,coef.aux.reftime,coef.aux.frq(ind),coef.aux.lind(ind),...
    coef.aux.lat,ngflgs,coef.aux.opt.prefilt);

% fit
fit = E*ap + conj(E)*am;

% mean (& trend)
u = nan*ones(size(tin));
whr = ~isnan(tin);
if coef.aux.opt.twodim
    v = u;
    if coef.aux.opt.notrend
        u(whr) = real(fit) + coef.umean;
        v(whr) = imag(fit) + coef.vmean;
    else
        u(whr) = real(fit) + coef.umean + ...
            coef.uslope*(t-coef.aux.reftime);
        v(whr) = imag(fit) + coef.vmean + ...
            coef.vslope*(t-coef.aux.reftime);
    end
else
    if coef.aux.opt.notrend
        u(whr) = real(fit) + coef.mean;
    else
        u(whr) = real(fit) + coef.mean + ...
            coef.slope*(t-coef.aux.reftime);
    end
    v = [];
end
fprintf('done.\n');

%%--------------------------------------------------------- 
function [t,opt] = ut_rcninit(tin,args)
% UT_RCNINIT()
% parse inputs and options for UT_RECONSTR1()
% UTide v1p0 9/2011 d.codiga@gso.uri.edu

t = tin(:);
t(isnan(t)) = [];
opt.cnstit = [];
opt.minsnr = 2;
opt.minpe = 0;
while ~isempty(args)
    switch(lower(args{1}))
        case 'cnstit'
            opt.cnstit = args{2};    
            args(1:2) = [];
        case 'minsnr'
            opt.minsnr = args{2};
            args(1:2) = [];
        case 'minpe'
            opt.minpe = args{2};
            args(1:2) = [];
        otherwise 
            error(['ut_reconstr: unrecognized input: ' args{1}]);
    end
end

%%--------------------------------------------------------- 
function E = ut_E(t,tref,frq,lind,lat,ngflgs,prefilt)
% UT_E()
% compute complex exponential basis function
% inputs
%   t = times [datenum UTC] (nt x 1)
%   tref = reference time [datenum UTC] (1 x 1)
%   frq = frequencies [cph] (nc x 1)
%   lind = list indices of constituents in ut_constants.mat (nc x 1)
%   lat = latitude [deg N] (1 x 1)
%   ngflgs = [NodsatLint NodsatNone GwchLint GwchNone] each 0/1
%       ([0 1 0 1] case not allowed, and not needed, in ut_E)
%   prefilt = 'prefilt' input to ut_solv
% output
%   E = complex exponential basis function [unitless] (nt x nc)
% UTide v1p0 9/2011 d.codiga@gso.uri.edu 

nt = length(t);
nc = length(lind);
if ngflgs(2) && ngflgs(4)
    F = ones(nt,nc);
    U = zeros(nt,nc);
    V = 24*(t-tref)*frq';
else
    [F,U,V] = ut_FUV(t,tref,lind,lat,ngflgs);
end
E = F.*exp(1i*(U+V)*2*pi);
if ~isempty(prefilt)
    P=interp1(prefilt.frq,prefilt.P,frq)';
    P( P>max(prefilt.rng) | P<min(prefilt.rng) | isnan(P) )=1;
    E = E.*P(ones(nt,1),:);
end

%%--------------------------------------------------------- 
function [F,U,V] = ut_FUV(t,tref,lind,lat,ngflgs)
% UT_FUV()
% compute nodal/satellite correction factors and astronomical argument
% inputs
%   t = times [datenum UTC] (nt x 1)
%   tref = reference time [datenum UTC] (1 x 1)
%   lind = list indices of constituents in ut_constants.mat (nc x 1)
%   lat = latitude [deg N] (1 x 1)
%   ngflgs = [NodsatLint NodsatNone GwchLint GwchNone] each 0/1
% output
%   F = real nodsat correction to amplitude [unitless] (nt x nc)
%   U = nodsat correction to phase [cycles] (nt x nc)
%   V = astronomical argument [cycles] (nt x nc)
% UTide v1p0 9/2011 d.codiga@gso.uri.edu 
% (uses parts of t_vuf.m from t_tide, Pawlowicz et al 2002)

nt = length(t);
nc = length(lind);
%% nodsat
if ngflgs(2) % none
    F = ones(nt,nc);
    U = zeros(nt,nc);
else
    if ngflgs(1) % linearized times
        tt = tref;
    else         % exact times
        tt = t;
    end
    ntt = length(tt);
    load('ut_constants.mat');
    [astro,~]=ut_astron(tt');
    if abs(lat)<5 
        lat=sign(lat).*5; 
    end
    slat=sin(pi*lat/180);
    rr=sat.amprat;
    j=find(sat.ilatfac==1);
    rr(j)=rr(j).*0.36309.*(1.0-5.0.*slat.*slat)./slat;
    j=find(sat.ilatfac==2);
    rr(j)=rr(j).*2.59808.*slat; 
    uu=rem( sat.deldood*astro(4:6,:)+sat.phcorr(:,ones(1,ntt)), 1);
    nfreq=length(const.isat); %#ok
    mat = rr(:,ones(1,ntt)).*exp(1i*2*pi*uu);
    F = ones(nfreq,ntt);
    ind = unique(sat.iconst);
    for i = 1:length(ind)
        F(ind(i),:) = 1+sum(mat(sat.iconst==ind(i),:),1);
    end
    U = imag(log(F))/(2*pi); % faster than angle(F)
    F=abs(F);
    for k=find(isfinite(const.ishallow))'
        ik=const.ishallow(k)+(0:const.nshallow(k)-1);
        j = shallow.iname(ik);
        exp1 = shallow.coef(ik);
        exp2 = abs(exp1);
        F(k,:)=prod(F(j,:).^exp2(:,ones(ntt,1)),1);
        U(k,:)=sum(U(j,:).*exp1(:,ones(ntt,1)),1);
    end
    F=F(lind,:)';
    U=U(lind,:)';
    if ngflgs(1) % nodal/satellite with linearized times
        F = F(ones(nt,1),:);
        U = U(ones(nt,1),:);
    end
end
%% gwch (astron arg)
if ngflgs(4) % none (raw phase lags not greenwich phase lags)
    if ~exist('const','var')
        load('ut_constants.mat','const');
    end
    [~,ader] = ut_astron(tref);
    ii=isfinite(const.ishallow); 
    const.freq(~ii) = (const.doodson(~ii,:)*ader)/(24);
    for k=find(ii)'
        ik=const.ishallow(k)+(0:const.nshallow(k)-1);
        const.freq(k)=sum(const.freq(shallow.iname(ik)).*shallow.coef(ik));
    end
    V = 24*(t-tref)*const.freq(lind)';
else 
    if ngflgs(3)  % linearized times
        tt = tref;
    else 
        tt = t;   % exact times
    end
    ntt = length(tt);
    if exist('astro','var')
        if ~isequal(size(astro,2),ntt)
            [astro,~]=ut_astron(tt');
        end        
    else
        [astro,~]=ut_astron(tt');
    end
    if ~exist('const','var')
        load('ut_constants.mat');
    end
    V=rem( const.doodson*astro+const.semi(:,ones(1,ntt)), 1);
    for k=find(isfinite(const.ishallow))'
        ik=const.ishallow(k)+(0:const.nshallow(k)-1);
        j = shallow.iname(ik);
        exp1 = shallow.coef(ik);
        V(k,:) = sum(V(j,:).*exp1(:,ones(ntt,1)),1);
    end
    V=V(lind,:)';
    if ngflgs(3)    % linearized times
        [~,ader] = ut_astron(tref);
        ii=isfinite(const.ishallow);
        const.freq(~ii) = (const.doodson(~ii,:)*ader)/(24);
        for k=find(ii)'
            ik=const.ishallow(k)+(0:const.nshallow(k)-1);
            const.freq(k)=sum( const.freq(shallow.iname(ik)).* ...
                shallow.coef(ik) );
        end
        V = V(ones(1,nt),:) + 24*(t-tref)*const.freq(lind)';
    end
end

%%--------------------------------------------------------- 
function [astro,ader] = ut_astron(jd)
% UT_ASTRON()
% calculate astronomical constants
% input
%   jd = time [datenum UTC] (1 x nt)
% outputs
%   astro = matrix [tau s h p np pp]T, units are [cycles] (6 x nt)
%   ader = matrix of derivatives of astro [cycles/day] (6 x nt)
% UTide v1p0 9/2011 d.codiga@gso.uri.edu 
% (copy of t_astron.m from t_tide, Pawlowicz et al 2002)

d=jd(:)'-datenum(1899,12,31,12,0,0);
D=d/10000;
args=[ones(size(jd));
      d;
      D.*D;
      D.^3];
sc= [ 270.434164,13.1763965268,-0.0000850, 0.000000039];
hc= [ 279.696678, 0.9856473354, 0.00002267,0.000000000];
pc= [ 334.329556, 0.1114040803,-0.0007739,-0.00000026];
npc=[-259.183275, 0.0529539222,-0.0001557,-0.000000050];
ppc=[ 281.220844, 0.0000470684, 0.0000339, 0.000000070];
astro=rem( [sc;hc;pc;npc;ppc]*args./360.0 ,1);
tau=rem(jd(:)',1)+astro(2,:)-astro(1,:);
astro=[tau;astro];
dargs=[zeros(size(jd));
       ones(size(jd));
       2.0e-4.*D;
       3.0e-4.*D.*D];
ader=[sc;hc;pc;npc;ppc]*dargs./360.0;
dtau=1.0+ader(2,:)-ader(1,:);
ader=[dtau;ader];

