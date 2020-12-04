dir = './dataExample/';
filename = ['DopplerSignal'];
[y,Fs] = audioread([dir filename '.wav']);

% filtre passe-haut pour eliminer les artefacts
fcArt = 150;% [Hz] Low-Pass filter to remove artefacts

fin=length(y);%3*60/(1/Fs);%
Xreal      =  y(1:fin,1);
Ximag      = -y(1:fin,2);
X_complex  = Xreal + 1i*Ximag;

EmboleAdapt = Study0_L3(X_complex, Fs, 1);


%% Spectrogram
Lt      = 64;
Lf      = 256;
Yf1     = spectrogram(X_complex,Lt,round(Lt*0.8),Lf,Fs);
timeYf1 = linspace(0,length(X_complex)/Fs,length(X_complex)/round(Lt*0.8)-1);
freq    = (0:Lf-1)/Lf*Fs;

fig = figure(1);

ax(1) = subplot(1,1,1);
imagesc(timeYf1, freq, abs(Yf1(1:Lf/2,:))), axis xy
caxis([0 3.5545])
%
hold on
% Plot of detections
if ~isempty(EmboleAdapt.pos); plot(repmat(EmboleAdapt.pos,length(freq),1), repmat(freq,length(EmboleAdapt.pos),1)','w--','linewidth',2);end
hold off
xlabel('Temps (s)'), ylabel('Frequence (Hz)')
xlim([0 timeYf1(end)])
ylim([0 Fs])