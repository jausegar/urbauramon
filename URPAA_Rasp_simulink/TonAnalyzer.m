%close all, clear all, clc;
% Note: It is mandatory to generate H matrix before launching tonality_cal
%[H,~]=GenerateFilters(fk,fs);
function [T,ff,LX,WT,Wgr]=TonAnalyzer(xr,fso,reH16k,imH16k) %#codegen
%load('H_16k.mat');
H16k=reH16k+1i*imH16k;

%% Tonalidad de un sonido 
%
% Siguiendo el modelo de Aures descrito en: 
% Berechnungsverfahren für den sensorischen Wohlklang beliebiger
% Schallsignale. Von W. Aures (Ver pag 135 Fig. 8).
%
% Extracción de componentes tonales usando Terhardt: 
% Algorithm for extraction of pitch and pitch salience from complex tonal signals
%
% Resumido en la tesis de Zhong Zhang y Merina Shrestha:
% Sound Quality User-defined Cursor Reading Control-Tonality Metric
%

%% 1) READ OR CREATE THE SOUND
%
%________________________________________________
% [x,fso] = audioread('fluteA.wav');
% 
% % new sampling frequency
 fs = 16000.;
 if(fso~=fs)
     xs = resample1(xr,fs,fso);
 else
     xs = xr(:);
 end

%________________________________________________


%% 2) COMPUTE THE SPECTRUM ANALYSIS
% Using a Hamming window
winlen = 80e-3; %80 ms recomended
wlen   = winlen*fs;
hop    = wlen/2;
K      = wlen/2+1;

% Como ejemplo trabajamos con la ventana 'l'
l = 5;

xf = enframe(xs,hanning(wlen),hop);
Xf = fft(xf');
Xf = Xf(1:wlen/2+1,:);%<--- A partir de aquí trabajamos sólo con 1 ventana

freq = (0:wlen-1)*fs/wlen;
freq = freq(1:wlen/2+1);

%fig_a = figure(1);
%fig_a.Position = ([150, 450, 700, 500]);
%subplot(4,1,1)
%plot(freq, abs(Xf(:,l)));
%ylabel('FFT')
%result_text = uicontrol('Style', 'text', 'String', 'TONALITY', ...
%    'Position', [80 15 280 20], 'FontSize', 14); 
 
% SPL mapping
% Tranform FFT values to SPL using normalization factor
dftmax = wlen/4;
Xf_spl = 96 + 20*log10(abs(Xf)/dftmax);

% Ejemplo
fk = (0:wlen-1)*fs/wlen;
fk = fk(1:wlen/2+1);

% Como ejemplo trabajamos con la ventana 20
%l = 20;
Xw = Xf_spl(:,l);


%% 3) EXTRACTION OF TONAL COMPONENTS

th = 7; % Threshold de 7dB
%
aux1 = circshift(Xw,2,1);
aux2 = circshift(Xw,3,1);
aux3 = circshift(Xw,-2,1);
aux4 = circshift(Xw,-3,1);
sel_ini = (Xw - aux1)>=th & (Xw - aux2)>=th & (Xw - aux3)>=th & (Xw - aux4)>=th;



% Cálculo de las frecuencias tonales "exactas"
% ft = fi + 0.46(Hz/dB)(Li+1 - Li-1)
% auxf1 = circshift(sel_ini,1,1);
% auxf2 = circshift(sel_ini,-1,1);
% 
% fk(sel_ini) = fk(sel_ini).' + 0.46*(Xw(auxf1) - Xw(auxf2)); %(3)

%% 4) Sound Pressure Level Excess
%
N = sum(sel_ini); % Número de elementos tonales
PeakIndex = find(sel_ini>0); % Vector que contiene los indices (número) de los elementos tonales
LX = zeros(1,N);

fkkhz = fk/1000;

% Calculo de la señal en intensidad.
XI = (10^-12) * 10.^(Xw/10); %

z = 13 * atan(0.76*fkkhz) + 3.5*atan(fkkhz/7.5).^2; %(6) Todas las frecuencias en barks
    
for u = 1: N % Para cada elemento tonal (fu)
    uu = PeakIndex(u);
    % Cálculo del nivel de excitacion LEv(fu)
    LE = zeros(1,N);
    for v = 1:N
        if (u == v)
            continue % If u==v dont compute LEv
        end

        vv = PeakIndex(v);
        
        if (fk(uu) <= fk(vv))
            s = 27; % (7a)
        else
            s = -24 - (0.23/fkkhz(vv))+(0.2*Xw(vv)); % (dB/Bark) (7b)
        end
        LE(v) = Xw(vv) - s * (z(vv) - z(uu)); % (5)
    end
    
    % Cálculo de la intensidad del ruido INu
    fmin = find(z >= z(uu)-0.5,1,'first');
    fmax = find(z <= z(uu)+0.5,1,'last');
    ff = [fmin:uu-3,uu+3:fmax];
    IN = sum(XI(ff));
    
    % Cálculo del humbral de audición
    LTH = 3.64*(fkkhz(uu))^-0.8 - 6.5*exp(-0.6*(fkkhz(uu)-3.3)^2) + 1e-3*(fkkhz(uu))^4; % (8)
    
    % Cálculo del exceso de SPL
    LX(u) = Xw(uu) - 10*log10( (sum(10.^(LE/20)))^2 + IN + 10^(LTH/10) ); %(4)
end

% Quitamos aquellas SPL que quedan debajo de 0
index0 = LX>0; % Calculamos los indices que necesitamos
PeakIndex0 = PeakIndex(index0);


%% 5) Tonal Weighting
ff = fk(PeakIndex0);
LX = LX(index0);
ff = ff(:);
LX = LX(:);
%%
% 5.1) W1
%flower = fk(PeakIndex0-2);
%fupper = fk(PeakIndex0+2);
%Zupper=13*atan(0.76*fupper/1000)+3.5*atan((fupper/7500).^2);
%Zlower=13*atan(0.76*flower/1000)+3.5*atan((flower/7500).^2);
%dzCorrect=(Zupper-Zlower);
zz = z(PeakIndex0+1) - z(PeakIndex0-1);

%zz=zz*0.05; %% REVISAR!!
w1 = 0.13./(zz+0.13); %(7)
w1=w1(:);

%%
% 5.2) W2
x = (1 + 0.2* (ff./0.7e3 + 0.7e3./ff ).^2).^(1/2);
w2 = (1./x).^0.29; % (9)

%%
% 5.3) W3
w3 = (1-exp(-LX./15)).^0.29; % (10)

%%
% 5.4) WT
%w1 = w1.^(1/0.29);
w2 = w2.^(1/0.29);
w3 = w3.^(1/0.29);

WT = sqrt(sum(w1.*w2.*w3).^2); % (11)

%fprintf('W_T: %.4e\n',WT);

%% 6) Cálculo del Loudness del ruido
%
NOISE = Xw(:)';
auxf1 = circshift(sel_ini,-2,1);
auxf2 = circshift(sel_ini,-1,1);
auxf3 = circshift(sel_ini,1,1);
auxf4 = circshift(sel_ini,2,1);

NOISE(sel_ini) = min(NOISE);
NOISE(auxf1) = min(NOISE);
NOISE(auxf2) = min(NOISE);
NOISE(auxf3) = min(NOISE);
NOISE(auxf4) = min(NOISE);

% NOISE(sel_ini) = -30;
% NOISE(auxf1) = -30;
% NOISE(auxf2) = -30;
% NOISE(auxf3) = -30;
% NOISE(auxf4) = -30;

% figure()
% plot(fk,NOISE);


MaxFilt=28;
% Note: It is mandatory to generate H matrix before launching tonality_cal
%[H,~]=GenerateFilters(fk,fs);
%load H_16k.mat;
%H=0.9639*H;

NoiseZinput=zeros(1,MaxFilt);
for ink=1:MaxFilt
    NoiseZinput(ink)=10*log10(sum((10.^(NOISE/10)).*(abs(H16k(ink,:).^2))));
end
field = 'f';
NoiseZinput(NoiseZinput<-70) = -70;
[NoiseN,NoiseNs,~]=DIN45631(NoiseZinput,field);


%% 7) Cálculo del Loudness de la señal
YdBZinput=zeros(1,MaxFilt);
for ink=1:MaxFilt
    YdBZinput(ink)=10*log10(sum((10.^(Xw'/10)).*(abs(H16k(ink,:).^2))));
end
YdBZinput(YdBZinput<-70) = -70;
[YdBN,YdBNs,~]=DIN45631(YdBZinput,field);


%% 8) Cálculo del Peso Tonal del Loudness
% 
Wgr=1-NoiseN/YdBN; 
%fprintf('W_Gr: %.4e\n',Wgr);

%% 9) Cálculo de la Tonalidad
%
%c = 1.0129;
c = 1.397; % Ajustado para que un tono puro de 1KHz 60 dB SPL de tonalidad de 1.00
T = c * (WT^0.29)*(Wgr^0.79);
%disp('_________');
%fprintf('TONALITY: %.5f\n',T);
%disp('_________');


%result_text.String = ['TONALITY = ', num2str(T)]; 
