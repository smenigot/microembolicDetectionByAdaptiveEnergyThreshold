function Embole = Study0_L3(X, Fs, display)
% X analytical signal
% Fs sampling frequency
% display : 0 no / 1 yes

%% Setting
Tolerance = 5;
Lwind     = 64;


%% Extraction of a priori information
OverLap = 20/100;
fcArt = 150;% [Hz] Low-Pass filter to remove artefacts

[Mpv_apriori, Med_apriori] = passe1(X, Lwind, OverLap, fcArt, Fs);


%% loop pass2
Lf      = 256;
OverLap = 80/100;
Tobs    = 10; % [s]
NObs    = Tobs*Fs;% 10 seconds extracted from the data

% Init
TotalNemb=0;
DynamicSNR               = zeros(round(length(X)/NObs),1);
SimilarityThresholdRatio = zeros(round(length(X)/NObs),1);
for k=1:fix(length(X)/NObs)
    EmboleTroncon(k).troncon = [];
    EmboleTroncon(k).Amp     = [];
    EmboleTroncon(k).posEch  = [];
    EmboleTroncon(k).pos     = [];
    EmboleTroncon(k).length  = [];
end

for k=1:fix(length(X)/NObs)
    x = X((k-1)*NObs+1:k*NObs,:);
    if std(x) > 0.01
%%  Information used for detection (Energy)
        [ESubNorm,ESubNormF,VarThreshold,std_flux] = Energies_L1(x, NObs, Lwind, OverLap, fcArt, Fs);
        
%% Indicator of Quality
        [DynamicSNR(k), Mpv_pos] = Quality2(ESubNorm,ESubNormF);
        
%% HITS detection
        [ThresholdHits,SimilarityThresholdRatio(k),PositionHits,DurationHits,ZHits] =...
            HITS_detection2(ESubNorm, ESubNormF, VarThreshold, std_flux(1), Mpv_pos);

%% Artefact detection
        [PositionArts,ThresholdArts,DurationArts,ZArts] = ...
            DetectionArts(ESubNorm,ESubNormF, Mpv_apriori(2));

%% Embolus detection
        [NEmb, EmbPos, EmbEner, EmbDur] = ...
            EmbDetector(PositionHits,PositionArts, DurationHits, DurationArts, ZHits,ZArts,ESubNorm,ESubNormF,Tolerance);
        
%% Displaying
        if display==1            
            timeX = 0:Tobs/(length(x)-1):Tobs;
            Yf1   = spectrogram(x(1:NObs,:),Lwind,fix(Lwind*OverLap),Lf,Fs);
            freq  = (0:Lf/2-1)/(Lf/2)*Fs;
            timeESub = 0:Tobs/(length(ESubNorm)-1):Tobs;
            
            fig = figure(1);
            
            subplot(511);
            plot(timeX,real(x),'k');
            
            subplot(5,1,2);
            imagesc(timeESub,freq,(abs(Yf1(1:Lf/2,:))));
            caxis([0 max(max(abs(Yf1(1:Lf/2,:))))]);
            axis xy;

            subplot(5,1,3);
            imagesc(timeESub,freq,(abs(Yf1(Lf:-1:Lf/2,:))));
            caxis([0 max(max(abs(Yf1(1:Lf/2,:))))])
            
            subplot(514);
            plot(timeESub,smooth(ESubNorm(:,1),10,'mean'),'k',...
                timeESub,ESubNormF(:,1),'m',...
                timeESub,ThresholdHits,'r');
            hold on;
            plot(timeESub(PositionHits),ESubNormF(PositionHits,1),'gs');
            plot(timeESub,(Mpv_pos(1)+mean(ESubNormF(:,1))), 'c-.',...
                timeESub,(Mpv_pos(2)+mean(ESubNormF(:,1))),'c-.');
            plot(timeESub,ThresholdHits(:,2),'r-.');
            plot(timeESub(EmbPos),ESubNormF(EmbPos,1),'r*')
            axis([0 10 0 max(ThresholdHits(:,1))]);
            hold off

            subplot(515);
            plot(timeESub,ESubNorm(:,2),'k',...
                timeESub,ESubNormF(:,2),'m',...
                timeESub,ThresholdArts,'r--');
            hold on;
            plot(timeESub(PositionArts),ESubNormF(PositionArts,2),'g*');
            plot(timeESub,ESubNormF(:,2)'+Mpv_apriori(2),'b',...
                timeESub,2*Med_apriori(2),'g');
            axis([0 10 0 max(ThresholdArts)]);
            hold off
        end
        
        
%         disp(['-------------------------Final RESULTS-------------------------' ])
%         disp(['Total number of 10s-blocks:                       ', int2str(k) ]);
%         disp(['SNR quality of the Doppler signal:         ', int2str(DynamicSNR(k)) ]);
%         disp(['Hits Threshold similarity (%):            ', int2str(SimilarityThresholdRatio(k)) ]);
%         disp(['Number of Detected ARTEFACTS:             ', int2str(length(PositionArts)) ]);%disp(NHits)
%         disp(['Number of Detected EMBOLUS:               ', int2str(NEmb) ]);%disp(NHits)
%         disp(['--------------------------------------------------------------' ])
        
%% Output extraction
        timeX=0:Tobs/(length(ESubNorm)-1):Tobs;

        for kEmb=(1:NEmb)
            EmboleTroncon(k).troncon(kEmb) = k;
            EmboleTroncon(k).Amp(kEmb) = EmbEner(NEmb);
            EmboleTroncon(k).posEch(kEmb) = EmbPos(NEmb);

            EmbtimePos=(k-1)*10+timeX(EmbPos(NEmb));
            EmboleTroncon(k).pos(kEmb) = EmbtimePos;
            EmboleTroncon(k).length(kEmb) = EmbDur(NEmb);
        end
        TotalNemb=TotalNemb+NEmb;
    end

end
%%
kEmb=1;
for k=1:length(EmboleTroncon)
    for k2=1:length(EmboleTroncon(k).troncon)
        Embole.troncon(kEmb) = EmboleTroncon(k).troncon(k2);
        Embole.Amp(kEmb)     = EmboleTroncon(k).Amp(k2);
        Embole.posEch(kEmb)  = EmboleTroncon(k).posEch(k2);
        Embole.pos(kEmb)     = EmboleTroncon(k).pos(k2);
        Embole.length(kEmb)  = EmboleTroncon(k).length(k2);
        kEmb=kEmb+1;
    end
end