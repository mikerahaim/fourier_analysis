% *********************************************************************** %
% Michael Rahaim
% Modified Examples from http://www.gaussianwaves.com/
%
% ------------------------------------------------------------------------
% Note on Power: 
%       Amplitude = 1/2 Peak to Peak
%       RMS of a sine = P-P / 2*sqrt(2)
%       For a real signal,
%           Power at f = Power at -f = (RMS)^2 / 2
%       Power calculation also accounts for the % of the FFT that the
%           signal accounts for since fft(x,NFFT) fills x with zeros is
%           length of x < NFFT.
% ------------------------------------------------------------------------
% Note on DFT:
%       When calculating the FFT and iFFT, Matlab's fft and ifft functions
%       do NOT use the 1/N and N factor that are found in text definitions
%       of the complex DFT.
% ------------------------------------------------------------------------
% Suppress warngngs:
%#ok<*UNRCH>
% *********************************************************************** %
clear all %#ok<CLALL>
close all

%% Parameters 
A             = sqrt(2);    % Amplitude
f             = 10;         % frequency of sine wave
overSampRate  = 16;         % oversampling rate
fs = overSampRate*f;        % sampling frequency

phase         = 0;%-1/2*pi; % desired phase shift in radians
nCyl          = 4;      % to generate five cycles of sine wave
NFFT          = 64;   % NFFT-point DFT	 	 

tolerance     = 1e-8;               % To avoid FP rounding errors
noise_var     = 0;                  % Noise variance
noise_std     = sqrt(noise_var);    % Noise standard deviation

%% Figure Settings
fig_c1 = 10;
fig_r1 = 700;
fig_w  = 460;
fig_h  = 400;
del_w  = 15;
del_h  = 100;
fig_c  = [fig_c1, fig_c1+fig_w+del_w, fig_c1+2*(fig_w+del_w), fig_c1+3*(fig_w+del_w)];
fig_r  = [fig_r1, fig_r1-(fig_h+del_h)];

STEM = 1; % 0 for regular plots, 1 for stem plots

%% Generate Signal / Info
t = 0:1/fs:(nCyl*1/f-1/fs);         % time base
x = A*cos(2*pi*f*t+phase);          % time domain signal
%x = x + 0.5*A*cos(2*pi*2*f*t+2*phase);
x = x + noise_std*randn(size(x));   % signal plus noise
L = length(x);                      % signal length


%% Generate DFT
X1     = fft(x,NFFT);               % Compute DFT using FFT	 	 
X2     = fftshift(X1);              % Shift the FFT values
Px     = X2.*conj(X2)/(NFFT*L);     % Power of each freq components	 	 
Px2    = X1.*conj(X1)/(NFFT*L);     % Single sided PSD

nVals1 = 0:NFFT-1;                  % DFT Sample points	 	 
nVals2 = (0:NFFT-1)/NFFT;           % Normalized DFT Sample points	 	 
nVals3 = (-NFFT/2:NFFT/2-1)/NFFT;   % DFT Sample points (Shifted)
fVals1 = fs*nVals3;                 % DFT frequencies
fVals2 = fs*(0:NFFT/2-1)/NFFT;	 	% Single sided DFT frequencies

% Get rid of values below the threshold tolerance (occur because of the
% rounding errors in floating point nums. Causes phase calculation errors)
% X1(abs(X1)<tolerance) = 0;
% X2(abs(X2)<tolerance) = 0;
X1(imag(X1) < tolerance & imag(X1) > -tolerance) = ...
    real(X1(imag(X1) < tolerance & imag(X1) > -tolerance));
X1(real(X1) < tolerance & real(X1) > -tolerance) = ...
    imag(X1(real(X1) < tolerance & real(X1) > -tolerance));
X2(imag(X2) < tolerance & imag(X2) > -tolerance) = ...
    real(X2(imag(X2) < tolerance & imag(X2) > -tolerance));
X2(real(X2) < tolerance & real(X2) > -tolerance) = ...
    imag(X2(real(X2) < tolerance & real(X2) > -tolerance));


%% Figure - Signal
figure('Position',[fig_c(1),fig_r(1),fig_w,fig_h]);
if (STEM) 
    stem(t,x);
else
    plot(t,x);
end
title('Input Signal');
%title(['Sine Wave f=', num2str(f), 'Hz']);
xlabel('Time(s)');
ylabel('Amplitude');


%% Figure - Real and Imaginary components
figure('Position',[fig_c(2),fig_r(1),fig_w,fig_h]);
subplot(2,1,1);
if (STEM) 
    stem(nVals1,real(X1));
else
    plot(nVals1,real(X1));	 	 
end
title('Double Sided FFT');	 	 
ylabel('Real DFT Values');
xlim([nVals1(1),nVals1(NFFT-1)]);
subplot(2,1,2);
if (STEM) 
    stem(nVals1,imag(X1));
else
    plot(nVals1,imag(X1)); 
end
%title('Double Sided FFT - Imaginary Compononts');	 	 
xlabel(['Sample points (', num2str(NFFT), '-point DFT)']);
ylabel('Imag DFT Values');
xlim([nVals1(1),nVals1(NFFT-1)]);

%% Figure - Magnitude
figure('Position',[fig_c(3),fig_r(1),fig_w,fig_h]);
if (STEM) 
    stem(nVals1,abs(X1));
else
    plot(nVals1,abs(X1));
end
title('Double Sided FFT - without FFTShift');	 	 
xlabel(['Sample points (', num2str(NFFT), '-point DFT)']);
ylabel('|DFT Values|');
xlim([0,NFFT-1]);

figure('Position',[fig_c(4),fig_r(1),fig_w,fig_h]);
subplot(2,1,1);
if (STEM) 
    stem(nVals2,abs(X1));
else
    plot(nVals2,abs(X1));
end
title('Double Sided FFT - without FFTShift');	 	 
%xlabel('Normalized Frequency');
ylabel('|DFT Values|');
subplot(2,1,2);
if (STEM) 
    stem(nVals3,abs(X2));
else
    plot(nVals3,abs(X2));
end
title('Double Sided FFT - with FFTShift');
xlabel('Normalized Frequency');
ylabel('|DFT Values|');
%xlim([nVals3(1),nVals3(NFFT-1)]);


%% Figure - Magnitude and Phase
figure('Position',[fig_c(1),fig_r(2),fig_w,fig_h]);
subplot(2,1,1)
if (STEM) 
    stem(fVals1,abs(X2),'b');
else
    plot(fVals1,abs(X2),'b');
end
title('Double Sided FFT - with FFTShift');	 	 
%xlabel('Frequency (Hz)');
ylabel('|DFT Values|');
xlim([fVals1(1),fVals1(NFFT-1)]);
subplot(2,1,2)
if (STEM) 
    stem(fVals1,(180/pi)*atan2(imag(X2),real(X2)),'b');
else
    plot(fVals1,(180/pi)*atan2(imag(X2),real(X2)),'b');
end
xlabel('Frequency (Hz)');
ylabel('Phase (deg)');
xlim([fVals1(1),fVals1(NFFT-1)]);
ylim([-180, 180]);
ax = gca;
ax.YTick = [-180,-90,0,90,180];


%% Figure - Power Spectral Density
figure('Position',[fig_c(2),fig_r(2),fig_w,fig_h]);
if (STEM) 
    stem(fVals1,Px,'b');
else
    plot(fVals1,Px,'b');
end
title('Power Spectral Density');	 	 
xlabel('Frequency (Hz)');
ylabel('Power');
xlim([-fs/2,fs/2]);

figure('Position',[fig_c(3),fig_r(2),fig_w,fig_h]);
plot(fVals1,10*log10(Px),'b'); % Stem plots look funny in dB...
title('Power Spectral Density');	 	 
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
xlim([-fs/2,fs/2]);


%% Figure - Single Sided Power Spectral Density
figure('Position',[fig_c(4),fig_r(2),fig_w,fig_h]);
if (STEM) 
    stem(fVals2,Px2(1:NFFT/2),'b','LineWidth',1);
else
    plot(fVals2,Px2(1:NFFT/2),'b','LineWidth',1);	 	 
end
title('One Sided Power Spectral Density');	 	 
xlabel('Frequency (Hz)');
ylabel('PSD');