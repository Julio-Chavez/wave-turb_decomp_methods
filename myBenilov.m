function [B]= myBenilov(U,V,W,P,nfft,doffp,SurfType,fs,fcutoff,rho, f_swell_low, f_swell_high)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Benilov and Fuyiskin, 1970 method
% for Reynolds Stress computation with strong waves
%
% INPUT:
% U         [m/s] u velocity signal
% V         [m/s] v velocity signal
% P         [m] or [dBar] pressure/elevation signal
% nfft      points to use for FFT, recommend power of 2
%            =5*60/dt_sec;% approximate 5 minute window for ig waves
%            =1*60/dt_sec, appoximate 1 minute window for swell waves
% doffp     [m] height of p measurement off bottom
% fs        [Hz] sampling frequency
% fcutoff   [Hz] high frequency cutoff, depends on fs and instrument depth
% rho       [kg/m^3] density
%
% OUTPUT: WaveStats structure file
% Structure file B with following variables
% B.See = Seet;
% B.Suu = Suut;
% B.Svv = Svvt;
% B.Sww = Swwt;
% B.Suv = Suvt;
% B.Suw = Suwt;
% B.Svw = Suvt;
% B.Sue = Suet;
% B.Sve = Svet;
% B.Swe = Swet;
% 
% B.S_uwave_wwave = S_uwave_wwave;
% B.S_vwave_wwave = S_vwave_wwave;
% B.Supwp = Supwp;
% B.Svpwp = Svpwp;
% B.upwp_bar = upwp_bar;
% B.vpwp_Bar = vpwp_bar;
% % fmt       (f) [Hz] frequency vector

% df        [Hz] frequency resolution is:  = fs/(nfft-1); fs is sampling frequency

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%
g = 9.81;    % gravity, m/s^2
max_attenuation_correction = 3; % normally the maximum attenuation correction should not be higher than 5

% check for vector direction, need to be colunm vectors
if ~iscolumn(U) || ~iscolumn(W) || ~iscolumn(P)
    error("The input signals, U, W, P, should be column vectors.")
end

% Check whether there are NaN values in array
if (all(isnan(U)) || all(isnan(W)) || all(isnan(P)) )
    error("The input signals, U, W, P, should have finite value.")
end

% interpolate out nans
P = naninterp(P);
U = naninterp(U);
V = naninterp(V);
W = naninterp(W);

% Detrend signals
P = P - mean(P);
U = U - mean(U);
V = V - mean(V);
W = W - mean(W);

if strcmp(SurfType,'pressure')
    % Determine the depth time series using preassure measurements
    P = P*1E4/(rho*g);
    depth = mean(P,'omitnan')*1E4/(rho*g) + doffp;  % This is the averaged water depth
elseif strcmp(SurfType,'elevation')
    % Determine the depth time series using preassure measurements
    depth = mean(P,'omitnan') + doffp;  % This is the averaged water depth
else
    error("Usage: SurfType is either''pressure'' or ''elevation''.")
end

if (depth <= 0.0)
    error("Warning! Inconsistent depth values (depth < 0).")
end 

%%%%%%%%%%%%%%%%%%%%%%
% correlation
df = fs/(nfft-1);   % frequency resolution
noverlap = nfft/2;
% ==================

% Autospectra
[Suu, fm] = cpsd(U,U,hann(nfft),noverlap,[],fs); 
[Svv, ~] = cpsd(V,V,hann(nfft),noverlap,[],fs);  
[Sww, ~] = cpsd(W,W,hann(nfft),noverlap,[],fs);  
[Spp, ~] = cpsd(P,P,hann(nfft),noverlap,[],fs); 
%Cross spectra
[Suv, ~] = cpsd(U,V,hann(nfft),noverlap,[],fs);
[Suw, ~] = cpsd(U,W,hann(nfft),noverlap,[],fs); 
[Svw, ~] = cpsd(V,W,hann(nfft),noverlap,[],fs); 
[Spu, ~] = cpsd(U,P,hann(nfft),noverlap,[],fs);  
[Spv, ~] = cpsd(V,P,hann(nfft),noverlap,[],fs); 
[Spw, ~] = cpsd(W,P,hann(nfft),noverlap,[],fs); 
% The result is now that  nansum(Suu)*df = var(U) - actually < var(U) because of detrending

% DEPTH CORRECTION & SPECTRAL WEIGHTED AVERAGES / STATS--------------------
% find correction/conversion rcoefs at each f-band 
% to calc sea surface elevation and convert velocities to pressure units 


omega = 2*pi*fm;   % This is radian frequency

warning off MATLAB:divideByZero
% k = get_wavenumber( omega , depth);
k = get_wavenumber( omega , depth);
warning on MATLAB:divideByZero


correction = ones(1,length(fm));

ii = find(fm <= fcutoff);
correction(ii) = 1E4 *cosh(k(ii)*depth)./(rho * g * cosh(k(ii)*doffp));
correction(fm>fcutoff) = 1;	% attenuation for frequencies less or equal fcutoff only
correction(correction < 1/max_attenuation_correction) = 1/max_attenuation_correction;
correction(correction > max_attenuation_correction) = max_attenuation_correction;

See = Spp .* (correction'.^2) ; % correct pressure for attentuation to get sea-surface elevation
Sue = Spu .* correction';
Sve = Spv .* correction';
Swe = Spw .* correction';

% find indices of freq bands
i_swell = find(fm>f_swell_low & fm<f_swell_high);

%Full raw Reynolds stresses
uu = sum(real(Suu)*df, 'omitnan');
uu_check = mean(U.*U); %should match
uv = sum(real(Suv)*df, 'omitnan');
uv_check = mean(U.*V); %should match
uw = sum(real(Suw)*df, 'omitnan');
uw_check = mean(U.*W); %should match
vv = sum(real(Svv)*df, 'omitnan');
vw = sum(real(Svw)*df, 'omitnan');
ww = sum(real(Sww)*df, 'omitnan');
ee = sum(real(See)*df, 'omitnan');
ue = sum(real(Sue)*df, 'omitnan');
ve = sum(real(Sve)*df, 'omitnan');
we = sum(real(Swe)*df, 'omitnan');

% Benilov and Filyushkin 1970 --------------------------------------------------------
% Su_wave_w_wave = Sue(f)(Swe(f)*)/See(f)
% Su’w’ = Suw – Su_wave_w_wave
% <u’w’> = integral(real(Su’w’(f)))df

% wave results
% auto
S_uwave_uwave = Sue .* conj(Sue) ./ See;
S_vwave_vwave = Sve .* conj(Sve) ./ See;
S_wwave_wwave = Swe .* conj(Swe) ./ See;
% cross
S_uwave_wwave = Sue .* conj(Swe) ./ See;
S_vwave_wwave = Sve .* conj(Swe) ./ See;
S_uwave_vwave = Sue .* conj(Sve) ./ See;

% turbulence results
% auto
Supup = abs(Suu) - abs(S_uwave_uwave);
Svpvp = abs(Svv) - abs(S_vwave_vwave);
Swpwp = abs(Sww) - abs(S_wwave_wwave);
% cross
Supwp = abs(Suw) - abs(S_uwave_wwave);
Svpwp = abs(Svw) - abs(S_vwave_wwave);
Supvp = abs(Suv) - abs(S_uwave_vwave);

% integrate wave over spectra, only in wave band
% auto
uwave_uwave_bar = nansum( real(S_uwave_uwave(i_swell)) .* df);
vwave_vwave_bar = nansum( real(S_vwave_vwave(i_swell)) .* df);
wwave_wwave_bar = nansum( real(S_wwave_wwave(i_swell)) .* df);
% cross
uwave_wwave_bar = nansum( real(S_uwave_wwave(i_swell)) .* df); % i think this is right (coherence not quadrature)
vwave_wwave_bar = nansum( real(S_vwave_wwave(i_swell)) .* df);
uwave_vwave_bar = nansum( real(S_uwave_vwave(i_swell)) .* df);

% integrate turbulence over entire spectra
% the 2 comes from a double sided spectra
% auto
upup_bar = 2*nansum( real(Supup) .* df);
vpvp_bar = 2*nansum( real(Svpvp) .* df);
wpwp_bar = 2*nansum( real(Swpwp) .* df);
% cross
upwp_bar = 2*nansum( real(Supwp) .* df); % i think this is right (coherence not quadrature)
% upwp_bar = nansum( abs(Supwp) .* df); I think this is incorrect
vpwp_bar = 2*nansum( real(Svpwp) .* df);
upvp_bar = 2*nansum( real(Supvp) .* df);

% Turbulent kinetic energy
TKEturb = 0.5*(upup_bar + vpvp_bar + wpwp_bar);

% condense results
B.See = See;
B.Suu = Suu;
B.Svv = Svv;
B.Sww = Sww;
B.Suv = Suv;
B.Suw = Suw;
B.Svw = Svw;
B.Sue = Sue;
B.Sve = Sve;
B.Swe = Swe;

B.S_uwave_uwave = S_uwave_uwave;
B.S_vwave_vwave = S_vwave_vwave;
B.S_wwave_wwave = S_wwave_wwave;
B.S_uwave_wwave = S_uwave_wwave;
B.S_vwave_wwave = S_vwave_wwave;
B.S_uwave_vwave = S_uwave_vwave;
B.Supup = Supup;
B.Svpvp = Svpvp;
B.Swpwp = Swpwp;
B.Supwp = Supwp;
B.Svpwp = Svpwp;
B.Supvp = Supvp;

% wave tensor
B.uwave_uwave_bar = uwave_uwave_bar;
B.vwave_vwave_bar = vwave_vwave_bar;
B.wwave_wwave_bar = wwave_wwave_bar;

B.uwave_vwave_bar = uwave_vwave_bar;
B.uwave_wwave_bar = uwave_wwave_bar;
B.vwave_wwave_bar = vwave_wwave_bar;

% reynolds stress tensor
B.upup_bar = upup_bar;
B.vpvp_bar = vpvp_bar;
B.wpwp_bar = wpwp_bar;

B.upwp_bar = upwp_bar;
B.vpwp_bar = vpwp_bar;
B.upvp_bar = upvp_bar;

% Raw Reynolds stresses
B.uu = uu;
B.uv = uv;
B.uw = uw; 
B.vv = vv;
B.vw = vw;
B.ww = ww;
B.ee = ee;
B.ue = ue;
B.ve = ve;
B.we = we;
B.uu_check = uu_check;
B.uv_check = uv_check;
B.uw_check = uw_check;

B.TKEturb = TKEturb;

% other info
B.fm = fm;
B.df = df;
B.nfft = nfft;
B.fcutoff = fcutoff;
B. depth = depth;
B.fs = fs;

%%
% subplot(1,2,1)
% loglog(fmt,abs(Supwp),fmt,abs(Suwt),fmt,abs(S_uwave_wwave))
% legend('<u''w''>','<uw>','<uwave wwave>','location','ne')
% ylabel('S'), xlabel('Hz')
% 
% subplot(1,2,2)
% loglog(fmt,abs(Svpwp),fmt,abs(Svwt),fmt,abs(S_vwave_wwave))
% legend('<v''w''>','<vw>','<vwave wwave>','location','ne')
% ylabel('S'), xlabel('Hz')



end

function X = naninterp(X)
% Interpolate over NaNs
% See INTERP1 for more info
X(isnan(X)) = interp1(find(~isnan(X)), X(~isnan(X)), find(isnan(X)));
end

function [k,err] = get_wavenumber(omega,depth,tolerance,num_iterations)
% VECTWAVENUM Linear dispersion relation wavenumber solver
% 
%	k = vectWavenum(omega,depth) solves for the wave-numbers 'k' for a
%	given depth that satisfy the dispersion relation:
%
%   	omega^2 = g*k*tanh(k*depth)
%	
%   using the Newton-Raphson method.  The variables omega and depth can be
%	arrays.
% 
%	k = vectWavenum(omega,depth,tolerance) assigns a tolerance value for
%	this calculation. If 'tolerance' is omitted, or entered as [], then a
%	default value of 1e-6 is used.  The solver continues iterating until
%	all values of 'k' meet this tolerance level or until it reaches the
%	maximum number of iterations.
% 
%	k = vectWavenum(omega,depth,tolerance,num_iterations) computes 'k'
%   within a given tolerance using a maximum number of interations:
%   'num_iterations'. If 'num_iterations' is omitted, or entered as [],
%   then a default value of 10 is used.
% 
%	[k,err] = vectWavenum(omega,depth,...) returns the approximation error
%   associated with each value 'k': err = abs( g*k*tanh(k*d) - omega^2 )
% 
%	S.D.Brenner, April 2018

%% Error checking:
narginchk(2,4);
nargoutchk(1,2);

% omega and depth should be either scalars, or arrays of the same size:
if all(size(omega) ~= size(depth)) && length(omega(:)) ~= 1 && length(depth(:)) ~= 1
    error('''omega'' and ''depth'' must be arrays of the same size or be of length 1');
end

% If tolerance and num_iterations are included, they should be scalars:
if nargin >= 3 &&  length(tolerance) >1
   error('If ''tolerance'' is included, it should be length 1');
end
if nargin == 4 && length(num_iterations) >1
    error('If ''num_iterations'' is included, it should be length 1');
end

%% Set optional argument values
g = 9.81;
args = {1e-6,30}; % Default values
if nargin >= 3 && ~isempty(tolerance)
    args{1} = tolerance;
end
if nargin == 4 && ~isempty(num_iterations)
    args{2} = num_iterations;
end
tol = args{1};
num_iter = args{2};

% If one of omega or depth are scalars, they should be converted to vectors
% of equal size:
if length(omega(:)) == 1
    omega = omega* ones(size(depth));
end
if length(depth(:)) == 1
    depth = depth* ones(size(omega));
end


%% Initial guess:

% from CEM II-1-11:
% kn = (w.^2/g) .* ( tanh( (w.^2/g).*d ) ).^(-1/2)

% from Fenton and McKee, 1990:
k0 = omega.^2/g; % deep water wavenumber
kn = k0 .* ( tanh((k0.*depth).^(3/4)) ).^(-2/3);

% Error for initial guess
err = abs( omega.^2 - g*kn.*tanh(kn.*depth) ); 

%% Use Newton's method: k_(n+1) = k_n + f(k_n)/f'(k_n)

converged = 0; % flag to identify if the iteration has converged with tolerance
for n = 1:num_iter
    
    % on each loop, re-compute only the k's that don't meet the tolerance
    idx = find(err > tol); 
    
    fkn = g*kn(idx).*tanh(kn(idx).*depth(idx)) - omega(idx).^2;
    dfkn = g*tanh(kn(idx).*depth(idx)) + g*depth(idx).*kn(idx).*sech(kn(idx).*depth(idx)).^2;
    kn(idx) = kn(idx) - fkn./dfkn;
    
    % Estimation error
    err = abs( omega.^2 - g*kn.*tanh(kn.*depth) );
    
    % Tolerance check
    if max(err(:)) <= tol 
        converged = 1;
        break;
    end
end
if converged == 0
    warning('Method failed to a tolerance of %2.2d converge in %g iterations.  Maximum error of %2.2d recorded. Considering decreasing ''tolerance'' or increasing ''num_iterations''',tol, num_iter, max(err(:)));
end

% Return result:
k = kn;

% For any places with zero depth, set k to zero:
k(depth==0) = 0;


end

function k = waveno(omega,depth)

%
% Matlab function to return the wavenumber (m^-1) of a surface-gravity wave 
%      given the inputs of angular frequency (rad/s) and depth (m)
%
%      k = wavenumber( omega, depth )
%
% linear wave theory dispersion relation:  omega^2 = g*k * tanh(k*h)
% solved using the approximate explicit solution of J. Guo (2001)

g = 9.81;   % gravity in m/s/s
% Deep water limit
kdeep = omega.^2./g; 

k = kdeep .* ( 1 - exp( -(kdeep.*depth).^1.25 ) ).^(-.4);
end

% function A = calculate_fft(X,nfft,fs)
%   - Falk Feddersen
%   - returns matrix A(num,nfft)  of fourier coefficients
%   - this uses 75% overlapping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% notes Justin Rogers 2/27/13
% 
% Input:
% X     data input
% nfft  length of data to use for fft
% 
% Output:
% A     fourier coefficients A(num,nfft), where num = floor(4*n/nfft)-3;
% 
% 1. revised to allow row vector input for X
%

function output = VelSpectra(U, V, W, P, fs, nfft)
Amu = calculate_fft2(U,nfft);    % this is w/ 75% overlapping
Amv = calculate_fft2(V,nfft);
Amw = calculate_fft2(W,nfft);
Amp = calculate_fft2(P,nfft);

df = fs/(nfft-1);   % This is the frequency resolution
nnyq = nfft/2 +1;

output.fm = (0:nnyq-1)*df;

Suu = mean( Amu .* conj(Amu), 'omitnan' ) / (nnyq * df);  % need to normalize by freq resolution 
output.Suu = Suu(1:nnyq);

Svv = mean( Amv .* conj(Amv), 'omitnan' ) / (nnyq * df);  % need to normalize by freq resolution 
output.Svv = Svv(1:nnyq);

Sww = mean( Amw .* conj(Amw), 'omitnan' ) / (nnyq * df);  % need to normalize by freq resolution 
output.Sww = Sww(1:nnyq);

Suv = mean( Amu .* conj(Amv), 'omitnan' ) / (nnyq * df);  % need to normalize by freq resolution 
output.Suv = Suv(1:nnyq);

Suw = mean( Amu .* conj(Amw), 'omitnan' ) / (nnyq * df);  % need to normalize by freq resolution 
output.Suw = Suw(1:nnyq);



Svw = mean( Amv .* conj(Amw), 'omitnan' ) / (nnyq * df);  % need to normalize by freq resolution 
output.Svw = Svw(1:nnyq);

Spp = mean( Amp .* conj(Amp), 'omitnan' ) / (nnyq * df);  % need to normalize by freq resolution 
output.Spp = Spp(1:nnyq);

Spu = mean( Amp .* conj(Amu), 'omitnan' ) / (nnyq * df);  % need to normalize by freq resolution 
output.Spu = Spu(1:nnyq);

Spv = mean( Amp .* conj(Amv), 'omitnan' ) / (nnyq * df);  % need to normalize by freq resolution 
output.Spv = Spv(1:nnyq);

Spw = mean( Amp .* conj(Amw), 'omitnan' ) / (nnyq * df);  % need to normalize by freq resolution 
output.Spw = Spw(1:nnyq);
end

function A = calculate_fft2(X,nfft)

[n,m]=size(X);

if m>n
   X=X';
   [n,~]=size(X);
end
    
num = floor(4*n/nfft)-3;

%[n nfft num]

X = X-mean(X);
% sumXt = (X'*X)/n;

%WIN = hanningwindow(@hamming,nfft);
jj = (0:nfft-1)';
WIN = 0.5 * ( 1 - cos(2*pi*jj/(nfft-1)));


A = zeros(num,nfft);

% set it up so that SQR(|A(i,:)|^2) = sum(X^2)

varXwindtot = 0;

for i=1:num
      istart = ceil((i-1)*(nfft/4)+1);
      istop = floor(istart+nfft-1);
      Xwind = X(istart:istop);
      lenX = length(Xwind);
      Xwind = Xwind - mean(Xwind);  % demean.   detrend?
      varXwind =( Xwind'*Xwind)/lenX;            
      Xwind = detrend(Xwind);

      varXwindtot = varXwindtot + varXwind;
      Xwind = Xwind .* WIN;
      tmp = ( Xwind'*Xwind)/lenX;
      if (tmp == 0)
	 Xwind = Xwind * 0.0;
      else
	 Xwind = Xwind*sqrt(varXwind/tmp);
      end
      A(i,:) = fft(Xwind')/sqrt(nfft);
%       meanA2 = mean( A(i,:) .* conj(A(i,:)));
%      [meanA2 varXwind  meanA2/varXwind]
%      A(i,:) = A(i,:);    %  no longer necessary  * (varXwind / meanA2);   % make so above condition holds
end

%meanA2 = mean( mean( A .* conj(A) ));

%disp('--here is total---')
%[sumXt varXwindtot/num  meanA2 ]
end