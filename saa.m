% © 2018 Digital Signal Processing Laboratory, University of the Philippines Diliman

function A = new_saa(WAV_FILE, Freq, lf, hf)

[sig, Fs] = audioread(WAV_FILE);
disp(['Loading ' WAV_FILE]);
sig = medfilt1(sig,5,'truncate');   % Signal smoothing for sudden spikes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%Signal Filtering to divide Freq%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Freq = [80 100 125 160 200 250 315 400 500 630 :: 800 1000 1250 1600 2000 2500 3150 4000 5000];
% Divide the freq into two

% Shift signal so it is centered at zero 
S_buff = sig;
S_min = min(S_buff);
S_max = max(S_buff);
add = abs(S_min + S_max)/2;
sig = sig + add;

% Normalization of signal limit from 1 to -1
S_buff = sig;
S_min = min(S_buff);
S_max = max(S_buff);
sig = sig/max(S_min,S_max);        
sig=sig-mean(sig(:));
min(sig);
max(sig);
[thresh_start, thresh_off] = stability(sig);    % Taking the stable value(dB SPL) of the signal
S_time = sig;

% Divide freq into two
freq1 = Freq(1:10);
freq2 = Freq(11:19);

% Apply filters
% Note high freq has lower energy levels, hence higher threshold than 25dB
fc1 = 750;
fc2 = 700;
[b1, a1] = butter(10, fc1/(Fs/2));
[b2, a2] = butter(10, fc2/(Fs/2), 'high');
S_time1 = filter(b1,a1,S_time);      %Low Freq
S_time2 = filter(b2,a2,S_time);      %High Freq

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Starting and End point determination %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Calculating points of interest');

% Initialization
s = 32; % mV/Pa
p0 = 20e-6;
jmp = 0.1;
x = 0.1;		% Start of signal
i = 0;
j = 0;

% Find the point where stable state will start
while(i~=1)
    j = j + 1;
	S_buff = (rms(S_time(round(x*44100):round(((x*44100)+(0.1*44100)))))^2)*1e3;	%100 ms increments
	S1 = 20*log10(S_buff/(p0*s));
	x = x + jmp;
		if S1>thresh_start  %Through tests 80
			i = 1;
% 			ansx = sprintf('Starting point of generation at %d seconds with %d dB', x, S1);
% 			disp(ansx)
		else
			i = 0;
		end
	end
	
i = 0;
jmp = 0.01;

% Find the point where stable where pink noise generation will turn off
while(i~=1)
	S_buff = (rms(S_time(round(x*44100):round(((x*44100)+(0.1*44100)))))^2)*1e3;	%100 ms increments
	S1 = 20*log10(S_buff/(p0*s));
% 	A(j,1) = S1;
% 	j = j + 1;
	x = x + jmp;
		if S1<thresh_off  %Through tests 77
			i = 1;
% 			ansx = sprintf('Pink Noise generation will stop at %d with %d dB', x, S1);
% 			disp(ansx)
		else
			i = 0;
		end
end

i = 0;
% S_buff = (rms(S_time(round(x*44100):round(((x*44100)+(0.1*44100)))))^2)*1e3;	%100 ms increments
% S1 = 20*log10(S_buff/(p0*s));
jmp = 0.001;
n0 = round(x*Fs);
x1 = x;
x2 = x;

% Find the points where the low/high frequency will reach their certain threshold
while(i~=1)
	x1 = x1 + jmp;
	S_buff = (rms(S_time1(round(x1*44100):round(((x1*44100)+(0.1*44100)))))^2)*1e3;	%100 ms increments
	S21 = 20*log10(S_buff/(p0*s));
		if S1-S21>=lf   %Through tests
			i = 1;
		else
			i = 0;
		end
end
i = 0;
n1 = round(x1*Fs);
% ansx = sprintf('End of measurement for low frequencies at point %d with %d dB', x1, S21);
% disp(ansx)

while(i~=1)
	x2 = x2 + jmp;
	S_buff = (rms(S_time2(round(x2*44100):round(((x2*44100)+(0.1*44100)))))^2)*1e3;	%100 ms increments
	S22 = 20*log10(S_buff/(p0*s));
		if S1-S22>=hf   %Through tests
			i = 1;
		else
			i = 0;
		end
end
n2 = round(x2*Fs);
% ansx = sprintf('End of measurement for high frequencies at point %d with %d dB', x2, S22);
% disp(ansx)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%Short-time Fourier Analysis%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WINDOW = hamming(20);
NOVERLAP = round(length(WINDOW)*0.6);
[S_dB1,~,T1] = spectrogram(S_time1(n0:n1),WINDOW,NOVERLAP,freq1,Fs);
% S_dB1 = 20*log10(abs(S_dB1)+1);
[S_dB2,~,T2] = spectrogram(S_time2(n0:n2),WINDOW,NOVERLAP,freq2,Fs);
% S_dB2 = 20*log10(abs(S_dB2)+1);
% S_dB = log(hilbert(abs(S_time(round(n1):round(n2)))));
% T = 0:1:length(S_dB);
% T = T(1:length(S_dB));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%Linear Approximation Code%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[m,~] = size(S_dB1);
[n,~] = size(S_dB2);
d = zeros(m+n,1);
for i=1:m+n
    if i<=m
        env = abs(hilbert(S_dB1(i,:)));     % Hilbert transform is said to be effective in taking the decay rate of a signal
        buffx = polyfit(T1,log(env),1);     % https://www.mathworks.com/matlabcentral/answers/282827-how-to-find-the-decay-rate-of-a-decaying-signal
    else
        env = abs(hilbert(S_dB2(i-m,:)));
        buffx = polyfit(T2,log(env),1);
    end
    d(i,1) = -buffx(1);
end

%%%% Sabine Formula %%%%
% disp('Computing Absorption Coefficient...');

V = 58*3;
T = 22;						%Room tempt is 22 degree Celcius
c = 20.047*sqrt(273.15+T);
A = 0.9210*((V*d(:,1))/c);
