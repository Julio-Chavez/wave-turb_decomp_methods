function [Bricker] = myBricker(fs,nfft,inertialRange,waveRange,Data)
%%=========================================================================
% This function provides code to separate wave and turbulence components of 
% Reynold Stresses by the phase method [Bricker and Monismith, 2007]

% This function was adapted from C.E. Reimers for the manuscript published 
% in JGR-Oceans by Reimers and Fogaren (2021). 

%*** INPUTS:
% fs = sampling frequency (Hz) 
% nfft = number of fft points
% inertialRange = start and end of inertial range (Hz)
% waveRange = start and end of wave peak (Hz)
%   Data = data after rotation and time lag correction for oxygen 
%        col 1: horizontal velocity in dominate direction of waves (cm/s)
%        col 2: horizontal velocity perpendicular to u (cm/s)
%        col 3: vertical velocity after rotation to minimize wave           
%            components of w (cm/s)


% *** OUTPUTS:
% Bricker: structure with the Reynold stresses summary 
% TKEturb: turbulent kinetic energy (cm2 s-2); calculated from variances
%   which removed wave components using the phase method 

%%=========================================================================
% Check that the frequency is a scalar greater than zero
if ~isscalar(fs) && (fs <= 0); error("Hz is a scalar value greater than zero."); end

dt = 1/fs; % sampling interval (s)

% Time series of velocity and oxygen measurements
u = Data(:,1); % horizontal velocity in dominate direction of waves, cm/s units
v = Data(:,2); % horizontal velocity perpendicular to u
w = Data(:,3); % vertical velocity after rotation to minimize wave components of w

% Fluctuations velocity
u_nonan = naninterp(u);
v_nonan = naninterp(v);
w_nonan = naninterp(w);

% Time series of detrended velocity and oxygen measurements
up = u_nonan - mean(u_nonan); 
vp = v_nonan - mean(v_nonan);
wp = w_nonan - mean(w_nonan);


% dfold = 1/(dt*(length(wp)-1));
% nt = size(wp,1);
% 
% f = (0:dfold:dfold*(length(wp)-1)); % frequency vector
% nny = ceil((length(wp)./2));
nny = ceil(nfft./2);

% Fourier transforms of up, vp, and wp. Velocities changed from cm/s to m/s
% Amu = fft(up)./(nt-1); 
% Amv = fft(vp)./(nt-1);
% Amw = fft(wp)./(nt-1);

Amu = amplitudeSpectrumAveraged(up, nfft);
Amv = amplitudeSpectrumAveraged(vp, nfft);
Amw = amplitudeSpectrumAveraged(wp, nfft);
df = 1/(dt*(nfft-1));
% Compute power densities
% Sww  =2.*abs(Amw(1:nny)).^2/df;
% Suu = 2.*abs(Amu(1:nny)).^2/df;
% Svv = 2.*abs(Amv(1:nny)).^2/df;
% fm = f(1:nny);
noverlap = nfft/2;
[Suu, fm] = cpsd(up,up,hann(nfft),noverlap,[],fs);
[Svv, ~] = cpsd(vp,vp,hann(nfft),noverlap,[],fs);  
[Sww, ~] = cpsd(wp,wp,hann(nfft),noverlap,[],fs);

% Compute phase, output in radians 
Uph = atan2(imag(Amu),real(Amu));  
Vph = atan2(imag(Amv),real(Amv));
Wph = atan2(imag(Amw),real(Amw));

% Computing the cross spectra
Suv = 2.*real(Amu(1:nny).*conj(Amv(1:nny)))/(df);
Suw = 2.*real(Amu(1:nny).*conj(Amw(1:nny)))/(df);
Svw = 2.*real(Amv(1:nny).*conj(Amw(1:nny)))/(df);

waverange = find(fm > waveRange(1) & fm < waveRange(2)); %indices
fWave = fm(waverange); %frequencies for waves

% Calculate range of freq for linear interpretation (set limits):
%   must adjust for each data set

%range for turb slope, adjust depending on where wave peak sits
interprange = find(fm >waveRange(2) & fm < inertialRange(2));  
anchor = find(fm < waveRange(1) & fm > inertialRange(1));
interprange = cat(1,anchor,interprange);
fInter = fm(interprange');
Suu_inter = Suu(interprange);
Suu_wave = Suu(waverange);

%check on peak assignments
figure(1)
subplot(1,2,1)
loglog(fm,Suu,fWave,Suu_wave,'r-',fInter,Suu_inter','g-x');  

% Linear interpolation of turbulent spectra beneath wave peak
F = log10(fInter)'; 
Sl = log10(Suu_inter);
subplot(1,2,2)
plot(F,Sl,'x');
[P(1:2),s] = polyfit(F,Sl,1);
    slope = P(1);
    interc = P(2);
logy = (slope*(log10(fWave)) + interc);   
y = 10.^(logy);

figure(2)
subplot(2,2,1)
loglog(fm,Suu,'k-',fWave, Suu_wave,'r-',fWave,y,'b--')
title('Suu')
ylabel('S_u_u (m^2 s ^-^2 Hz^-^1)');
xlabel('f (Hz)');

[g,h] = size(Suu_wave);
Suu_wavecomp = zeros(g,h);
y = y';
for b = 1:g   
   Suu_wavecomp(b) = (Suu_wave(b) - y(b));  
end
%wave fourier component
Amuu_wave = sqrt((Suu_wavecomp+0j).*(df));
% pause

%now for vv    
Svv_inter = Svv(interprange);
Svv_wave = Svv(waverange);

%Linear interpolation over turbulent spectra
Sl = log10(Svv_inter);
[P(1:2),s] = polyfit(F,Sl,1);
    slope = P(1);
    interc = P(2);
logy = (slope*(log10(fWave)) + interc);   
y = 10.^(logy);

subplot(2,2,2)
loglog(fm,Svv,'k-',fWave,Svv_wave,'r-',fWave,y,'b-')
title('Svv')
ylabel('S_u_u (m^2 s ^-^2 Hz^-^1)');
xlabel('f (Hz)');
[g,h] = size(Svv_wave);
Svv_wavecomp = zeros(g,h);
y = y';
for b=1:g
    Svv_wavecomp(b) = (Svv_wave(b) - y(b));
end
%wave fourier component
Amvv_wave = sqrt((Svv_wavecomp+0j).*(df));      
% pause

%now for ww    
Sww_inter = Sww(interprange);
Sww_wave = Sww(waverange);

%Linear interpolation over turbulent spectra
Sl = log10(Sww_inter);
[P(1:2),s] = polyfit(F,Sl,1);
    slope = P(1);
    interc = P(2);
logy = (slope*(log10(fWave)) + interc);   
y = 10.^(logy);

subplot(2,2,3)
loglog(fm,Sww,'k-',fWave,Sww_wave,'r-',fWave,y,'b-')
title('Sww')
ylabel('S_w_w (m^2 s ^-^2 Hz^-^1)');
xlabel('f (Hz)');

[g,h] = size(Sww_wave);
Sww_wavecomp = zeros(g,h);
y = y';
for b=1:g   
   Sww_wavecomp(b) = (Sww_wave(b) - y(b));  
end
%wave fourier component
Amww_wave = sqrt((Sww_wavecomp+0j).*(df));

% pause   %use to reset ranges if needed

%Wave Magnitudes
Um_wave = sqrt(real(Amuu_wave).^2 + imag(Amuu_wave).^2);
Vm_wave = sqrt(real(Amvv_wave).^2 + imag(Amvv_wave).^2);
Wm_wave = sqrt(real(Amww_wave).^2 + imag(Amww_wave).^2);
   
%Wave reynolds stresses
uw_wave = sum(Um_wave.*Wm_wave.*cos(Wph(waverange) - Uph(waverange)), 'omitnan');
uv_wave =  sum(Um_wave.*Vm_wave.*cos(Vph(waverange) - Uph(waverange)), 'omitnan');
vw_wave = sum(Vm_wave.*Wm_wave.*cos(Wph(waverange) - Vph(waverange)), 'omitnan');

%sum wave components
uu_wave = sum(Suu_wavecomp*df, 'omitnan'); 
vv_wave = sum(Svv_wavecomp*df, 'omitnan');
ww_wave = sum(Sww_wavecomp*df, 'omitnan');
 
%Full Reynolds stresses
uu = sum(real(Suu)*df, 'omitnan');
uu_check = mean(up.*up); %should match
uv = sum(real(Suv)*df, 'omitnan');
uv_check = mean(up.*vp); 
uw = sum(real(Suw)*df, 'omitnan');
uw_check = mean(up.*wp); 
vv = sum(real(Svv)*df, 'omitnan');
vw = sum(real(Svw)*df, 'omitnan');
ww = sum(real(Sww)*df, 'omitnan');

%Turbulent reynolds stresses corrected
upup = uu - uu_wave;
vpvp = vv - vv_wave;
wpwp = ww - ww_wave;
upwp = uw - uw_wave;
upvp = uv - uv_wave;
vpwp = vw - vw_wave;

% fullfluxes = [uu vv ww uw uv vw];
% tfluxes = [upup vpvp wpwp upwp upvp vpwp];

TKEturb = 0.5*(upup + vpvp + wpwp); % in units cm2 s-2 

% condensing results in structure
Bricker.uu_wave = uu_wave;
Bricker.vv_wave = vv_wave;
Bricker.ww_wave = ww_wave;
Bricker.uv_wave = uv_wave;
Bricker.uw_wave = uw_wave;
Bricker.vw_wave = vw_wave;
Bricker.uu = uu;
Bricker.vv = vv;
Bricker.ww = ww;
Bricker.uv = uv;
Bricker.uw = uw;
Bricker.vw = vw;
Bricker.uu_check = uu_check;
Bricker.uv_check = uv_check;
Bricker.uw_check = uw_check;
Bricker.upup = upup;
Bricker.vpvp = vpvp;
Bricker.wpwp = wpwp;
Bricker.upvp = upvp;
Bricker.upwp = upwp;
Bricker.vpwp = vpwp;
Bricker.TKEturb = TKEturb;

end

function Am = amplitudeSpectrumAveraged(Signal, nfft)

% signal divisions
Ndiv = floor(size(Signal,1)/nfft*2) - 1; %1; %
% window and compute fft over each division
hannwin = hann(nfft);
for j = 1:Ndiv
    idx_start = (j-1)*nfft/2 + 1;
    idx_end = idx_start + nfft - 1;
    Am = fft(hannwin.*Signal( idx_start:idx_end, :))./(nfft-1);
end

end
function X = naninterp(X)
% Interpolate over NaNs
% See INTERP1 for more info
X(isnan(X)) = interp1(find(~isnan(X)), X(~isnan(X)), find(isnan(X)));
end