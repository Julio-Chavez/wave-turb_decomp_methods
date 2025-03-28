% This MATLAB script performs wave-turbulence decomposition on raw flow
% measurement data. It includes DMD, Coherence (Benilov), Phase (Bricker),
% and EEMD (Huang).
% It uses a synthetic dataset composed of trubulence and linear waves from
% Jonswap spectrum. % the turbulence signal we used to construct the synthetic 
% dataset only had horizontal component so we couldn't calculate the 
% Reynolds stresses but all the functions work with several signals from
% ADV (i.e. U, V, W, P time series)

% Input:
%   - Raw flow measurement data loaded from 'data.mat'
%   - Sampling frequency, fs
%   - Wave frequency range
%   - For spectra: number of fft points, nfft and window size, win

% Output:
%   - Various data and plots.

% Dependencies
% This script needs MATLAB signal processing toolbox 23.2 or later
%   - Functions in DMD folder
%   - Functions in EEMD folder
%   - myBenilov.m
%   - myBricker.m


%% Load data
close
clear
clc
%==========================================================================
nfft = 2048; win = 1000;
load('SynthData1.mat', 'u_tot'); % 1xk time series of flow measurements (ideally in m/s)
DataRaw = u_tot;
fs = 10; % sampling frequency (Hz)
waveRange = [1.17e-1 6e-1]; % [low-freq_end hig-freq_end]
%==========================================================================
dt = 1/fs; 
[numSignals, numObs] = size(DataRaw);
t = (0:numObs-1)/fs;


%% plot spectra of fluctuation and adjust waveRange accordingly

DataFluc = DataRaw - mean(DataRaw);
% spectra of fluctuations
[Sxx, fm] = cpsd(DataFluc', DataFluc',hann(win),win/2,nfft,fs);

figure(4);
loglog(fm, Sxx, 'Color', 'k', 'LineWidth',1)
hold off
xlim([fm(1) fm(end)])
grid on
xlabel('Frequency (Hz)')
ylabel('PSD (m$^2$s$^{-2}$/Hz)')

% section of time series
figure(5);
ts = 500; ss = 1500;
plot(t(1,ts:ts+ss),DataFluc(1,ts:ts+ss), 'Color', 'k', 'LineWidth',1)
xlabel('Time (s)', 'Interpreter','latex')
ylabel('Velocity (m/s)', 'Interpreter','latex')
xlim([t(ts) t(ts+ss)])


%% Decomposition
clc
% ========================================================================
rawInputSignal = DataFluc; % Should not have mean
r =[1 131]; % range of modes to use in truncation
matParam.Lag = 1; % number of obs to lag
matParam.NumObsLag = 1500; % number of columns in time-delayed matrix
% ========================================================================

% POD spectrum

[X1, ~] = HankelMatrix(rawInputSignal,matParam.NumObsLag ,numObs - matParam.NumObsLag, 'column');
% Decompose time-delayed matrix using SVD
[U, S, V] = svd(X1, 'econ');
V_win = round(size(V, 1)/3);
Vspec = nan(size(V, 1), nfft/2 + 1);
for ii = 1:size(V, 1)
    % get power spectrum of each row of V (or each column of V')
    Vspec(ii,:) = cpsd(V(:,ii),V(:,ii),hann(V_win),round(V_win/2),nfft,fs);
end
figure(12);
[X,Y] = meshgrid(fm,1:size(V, 1));
surf(X,Y,Vspec, 'EdgeColor', 'none')
% colormap('sky')
colormap(flipud(gray))
set(gca,'xscale','log')
set(gca,'yscale','log')
view(2); % Set the viewing angle to show only the x and y axes
box on
grid off
xline(waveRange, ':')
yline(r,'--')
xlabel('Frequency (Hz)', 'Interpreter','latex')
ylabel('Mode $j$', 'Interpreter','latex')
xlim([fm(1) fm(end)])

%% DMD
[X1, X2, Out] = DMD_processing(rawInputSignal, r, matParam, dt);

% DMD spectrum
[P2, f2, P1, f1] = DMD_spectrum(Out.b, Out.lambda, dt);
% find the indices of values within 'waveRange' in the input vector 'f1' 
% for 1-sided spectra and 'f2' for two-sided spectra.
[in1,out1, ~, ~] = find_in_range(f1,waveRange);
[in2,out2, inPos2, outPos2] = find_in_range(f2,waveRange);

figure(101)
loglog(f1(in1), P1(in1), 'xr'); hold on
loglog(f1(out1), P1(out1), 'xk');
xline(waveRange, '--')
box on
ylabel('DMD spectrum', 'Interpreter','latex'); xlabel('Frequency (Hz)', 'Interpreter','latex')


[u_time_dynamics_wave, u_Xdmd_wave] = DMD_reconstruct(X1,Out.Phi(:,in2),Out.omega(in2), Out.b(in2),dt);

u_wave_dmd = [real(u_Xdmd_wave(1,:)) real(u_Xdmd_wave(2:end,end))'];

u_turb_dmd = rawInputSignal(1:end-1) - u_wave_dmd;

% Plotting wave time series (DMD)
figure(3);
hold on
s1 = 2000; s2 = 2800;
plot(t(1:s2-s1+1), rawInputSignal(1,s1:s2), 'Color', '#444444', 'LineStyle','-', 'LineWidth', 1)
plot(t(1:s2-s1+1), u_wave_dmd(s1:s2), 'Color', '#1982c4', 'LineStyle','-', 'LineWidth', 1)
plot(t(1:s2-s1+1), u_turb_dmd(s1:s2), 'Color', '#ff595e', 'LineStyle','-', 'LineWidth', 1)
hold off
xlim([t(1) t(s2-s1+1)])
ylabel('$u$ (m/s)', 'Interpreter','latex'); xlabel('Time (seconds)', 'Interpreter','latex')
legend('Raw signal', 'DMD wave', 'DMD turbulence', 'Interpreter','latex')
box on


% Main DMD spectra
figure(4);
[SxxXs, fmXs] = cpsd(rawInputSignal, rawInputSignal,hann(win),win/2,nfft,fs);
[SxxXfs, fmXfs] = cpsd(u_wave_dmd, u_wave_dmd,hann(win),win/2,nfft,fs);
[SxxXnfs, fmXnfs] = cpsd(u_turb_dmd, u_turb_dmd,hann(win),win/2,nfft,fs);

loglog(fmXs, SxxXs, 'k', 'LineWidth', 1.5); hold on
loglog(fmXfs, SxxXfs, 'Color', '#1982c4', 'LineStyle','-', 'LineWidth', 1)
loglog(fmXnfs, SxxXnfs, 'Color', '#ff595e', 'LineStyle','-', 'LineWidth', 1)

xline(waveRange, ':')
xlim([fmXs(1) fmXs(end)])
hold off
xlabel('Frequency (Hz)', 'Interpreter','latex')
ylabel('PSD (m$^2$s$^{-2}$/Hz)', 'Interpreter','latex')
legend('Raw','DMD wave', 'DMD turbulence', 'Interpreter','latex')


%% Check eigenvalues
[lambdaCompsIn, lambdaAngleIn, ~, omegaCompsIn] = eigvals_lambda_omega(Out.lambda(in2), Out.omega(in2));
[lambdaCompsOut, lambdaAngleOut, circleArcOut, omegaCompsOut] = eigvals_lambda_omega(Out.lambda(out2), Out.omega(out2));

[~, lambdaAngleInPos, ~, ~] = eigvals_lambda_omega(Out.lambda(inPos2), Out.omega(inPos2));
[~, lambdaAngleOutPos, ~, ~] = eigvals_lambda_omega(Out.lambda(outPos2), Out.omega(outPos2));

figure(10);
subplot(1,2,1)
plot(circleArcOut(1,:),circleArcOut(2,:),'k') % plot unit circle
hold on, grid on
scatter(lambdaCompsIn(1,:),lambdaCompsIn(2,:),'or'); hold on
scatter(lambdaCompsOut(1,:),lambdaCompsOut(2,:),'ok')
xlabel('Real($\lambda$)', 'Interpreter','latex'); ylabel('Imaginary($\lambda$)', 'Interpreter','latex')
axis equal
subplot(1,2,2)
hold on, grid on
scatter(omegaCompsIn(1,:),omegaCompsIn(2,:),'+r')
scatter(omegaCompsOut(1,:),omegaCompsOut(2,:),'+k')
xline(0); yline(0)
box on
xlabel('Real($\omega$)', 'Interpreter','latex'); ylabel('Imaginary($\omega$)', 'Interpreter','latex')


figure(12);
subplot(1,2,1)
plot(in2, lambdaAngleIn, '.r'); hold on
plot(out2, lambdaAngleOut, '.k');
xlabel('mode', 'Interpreter','latex'); ylabel('$|$angle($\lambda$)$|$ [$0$ $\pi$]', 'Interpreter','latex')
subplot(1,2,2)
plot(f2(in2), lambdaAngleIn, '.r'); hold on
plot(f2(out2), lambdaAngleOut, '.k');
xlabel('Frequency (Hz)', 'Interpreter','latex'); ylabel('$|$angle($\lambda$)$|$ [$0$ $\pi$]', 'Interpreter','latex')



%% Benilov
% the turbulence signal we used to construct the synthetic dataset only had
% horizontal component so we couldn't calculate the Reynolds stresses
load('SynthData1.mat', 'u_tot','eta', 'u_turb'); % 1xk time series of flow measurements (in m/s)
DataRaw(1,:) = u_tot;
DataRaw(2,:) = eta;
DataFluc(1,:) = DataRaw (1,:)- mean(DataRaw(1,:));
DataFluc(2,:) = DataRaw (2,:)- mean(DataRaw(2,:));

% ---------- Benilov ----------

doffp = 0.3; SurfType = "elevation"; fcutoff = 1; rho = 1030;
f_swell_low = .039; f_swell_high = 0.148; % for u
% NOTE: the second and third inputs in myBenilov can be replaced by the v'
% and w' components, respectively
[Ben]= myBenilov(DataFluc(1,:)',zeros(length(DataFluc(1,:)'),1), zeros(length(DataFluc(1,:)'),1),DataFluc(2,:)',...
    nfft,doffp,SurfType,fs,fcutoff,rho, waveRange(1), waveRange(2)); 

%  ---------- Bricker ----------
inertialRange = [0.0234 1.4375];
% NOTE: the second and third row in hte last input in myBricker can be replaced by the v'
% and w' components, respectively
[Bri] = myBricker(fs,nfft,inertialRange,waveRange, [DataFluc(1,:)',zeros(length(DataFluc(1,:)'),1), zeros(length(DataFluc(1,:)'),1),DataFluc(2,:)']);


[SuuRaw, fmRaw] = cpsd(DataFluc(1,:)', DataFluc(1,:)',hann(win),win/2,nfft,fs);
[SxxT, fm_T] = cpsd(u_turb - mean(u_turb), u_turb - mean(u_turb),hann(win),win/2,nfft,fs);



% --------- EEMD ------------
noisy_signal = DataFluc'; 
Nstd = .8; %.2
noisy_signal_std = std(noisy_signal);
noisy_gauss_std = Nstd * noisy_signal_std;
noise = noisy_gauss_std.*randn(size(noisy_signal));
Nstd = noisy_gauss_std ./ noisy_signal_std;
inData = noisy_signal + noise;
rslt=eemd(inData(:,1),Nstd(1),100);

EEMDwave = sum(rslt(:,5:6),2); % 5:6(.8)
EEMDturb = noisy_signal(:,1) - EEMDwave;

[SuuEEMD, fmEEMD] = cpsd(EEMDturb, EEMDturb,hann(win),win/2,nfft,fs);



%% plots
clr = ["#e6a176", "#00678a", "#984464", "#5eccab", "#c0affb", "#56641a","#85754d"];
fontsz = 20;
figure
loglog(fmRaw, abs(SuuRaw), 'Color','k', 'LineWidth',3); hold on
loglog(fm_T, SxxT, 'Color', '#0000ff', 'LineStyle','-', 'LineWidth',1.5); hold on
loglog(Ben.fm, abs(Ben.Supup),'-', 'Color',clr(7), 'LineWidth',1.5); hold on
loglog(fmEEMD, abs(SuuEEMD),'-', 'Color',clr(4), 'LineWidth',1.5); hold on
loglog(fmXnfs, SxxXnfs, 'Color','r', 'LineWidth',1.5); hold on
legend("$u' + \tilde{u}$","$u'$", "$u'$ Coherence","$u'$ EEMD", "$u'$ DMD", 'FontSize', fontsz, 'Location', 'southwest')
xlim([0 fs/2])
% ylim([1e-6 2e0])
grid on
xlabel('Frequency (Hz)','FontSize',fontsz)
ylabel('PSD (m$^2$s$^{-2}$/Hz)','FontSize',fontsz)
hold off
ax = gca;
ax.FontSize = fontsz;



























































