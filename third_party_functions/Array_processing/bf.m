function y=bf(eobj,th_look,ph_look,beamType,beamParam,Nbf,nFFT_bf)
%% ======   Function Description   ======
% This function generate and applies beamforming on multichannel data from a microphone array which is contained in an EARO object.
% The beamforming is performed in the frequency and spatial domains and the single array output is transformed back to time domain without any spatiality.
% version 1.0 - created on 20/07/14 by David Alon
%
%assumes that the eobj.data contains measurments in microphone array from a single source
%===== Syntax  ==================
%   y=bf(eobj)
%   y=bf(eobj,th_look,ph_look)
%   y=bf(eobj,th_look,ph_look,beamType)
%
%===== Input variables  =========
%   eobj-  Earo object containing eobj.data gatherd by the micropohne array assuming - 1 source, n_t samples , Q microphones.
%           The dataDomain are assumed to be 'Time' and 'Space'
%   th_look- beamformer look direction elevation [rad] , range [0:pi] , default value={pi/2}
%   ph_look- beamformer look direction azimuth [rad] , range [0:2*pi] , default value={0}
%   beamType- beamformer type indicator: (1) max DI ;(2) max WNG ; (3)  ; (4) . default value={0}
%   [Nbf]
%   [Snn_inv]
%===== Output variables  ========
%   w - beamformer output.
%
%===== Default arguments values =
% y=bf(eobj,{pi/2},{0},{1})
if nargin<2, th_look=pi/2; end
if nargin<3, ph_look=0; end
if nargin<4, beamType=1; end
if (nargin<6), if (~isempty(eobj.orderN)), Nbf=eobj.orderN; else Nbf=sqrt(length(eobj.micGrid.azimuth)/2)-1; end; end;          % Default = Beamforming order is the same as microphone array order
if (nargin<7),nFFT_bf=size(eobj.data,2); end
if nargin<5, beamParam=inf; a_max=inf; end;          % Default = no soft-limiting
% if (nargin<7)&&(beamType==3), Snn=repmat(eye(length(eobj.micGrid.azimuth)),[1,1,size(eobj.data,2)]); end

%% ======   Main  =======================

% nFFT=2^ceil(log2(size(eobj.data,2)));
nFFT=size(eobj.data,2); % signal FFT length (not necessarily the same as the nFFT_bf )
% [mic_sig_f,f_vec]=eobj.toFreq(nFFT);
% c=340.3;
% k_vec=(2*pi*f_vec)/c;

L=length(th_look);
y_out_length=size(eobj.data,2);   % array output is in the same dimensions as the input
w_ir_length=1;                    % assumes that the eobj contains a signal which is sufficiently zero padded
useSOFiA = false;                     % change to true if you wish to generate the radial functions using SOFiA toolbox
% instead of EARO's internal func.




% Transform to frequency domain, if needed
if strcmp(eobj.dataDomain{1},'TIME')
    %temporal zero padding
    eobj.data=cat(2,eobj.data,zeros(size(eobj.data,1),(w_ir_length-1),size(eobj.data,3)));
    [eobj,fVec]=eobj.toFreq(nFFT);
else
    %       nFFT=size(eobj.data,2);
    fVec = (0:(nFFT-1))*(eobj.fs/nFFT);
end


% Compute speed of sound
if ~isempty(eobj.avgAirTemp) && isnumeric(eobj.avgAirTemp)
    c = 331.3 * sqrt(1+(eobj.avgAirTemp/273.15));
else
    c = 343.5; % default value;
end

% Construct k vector
kVec=double(fVec(1:ceil((nFFT+1)/2))*2*pi/c);
    fVec_bf = (0:(nFFT_bf-1))*(eobj.fs/nFFT_bf);
kVec_bf=double(fVec_bf(1:ceil((nFFT_bf+1)/2))*2*pi/c);




% p_t=squeeze(eobj.data).'; %obj.data:[J,nData,Q] p_t:[Q,nData]
% PP=fft(p_t.').';
% nFFT=size(PP,2);
% P=PP(:,1:ceil((nFFT+1)/2)); %take the single side fft - for even or odd nData
%
% P=shiftdim(mic_sig_f.data(1,1:ceil((nFFT+1)/2),:),1).'; %take the single side fft - for even or odd nData
% P(:,1) = 0; % put zero DC (f=0).  P:[Q,~nFFT/2]
%
% f_vec=[0:(nFFT-1)]*eobj.fs/nFFT;
% Transforming data to the SH domain


%% ??????????????????%%%%
% if beamType==2, % calculates Snn for MVDR beamformer
%     if ~exist(Snn,'var')
%         [mic_n_f,f_vec]=mic_n.toFreq(nFFT);
%         Pn=squeeze(mic_n_f.data(1,1:ceil((nFFT+1)/2),:)).'; %take the single side fft - for even or odd nData
%         P(:,1) = 0; % put zero DC (f=0).  P:[Q,~nFFT/2]
%     end
% end

%%

% calculating beamforming weights in the spatial domain
%%
% create tmp earo with metadata
% tmpe = earo;
% tmpe.orderN = eobj.orderN;
switch beamType
    case [{6},{10}]
        % performing bf in spatial domain
        beamType_space=beamType/2;
        a_max=beamParam{1};
        Snn_inv=beamParam{2};%permute((:,:,1:ceil((nFFT+1)/2)),[3,1,2]);
        % Transform to Space domain, if needed
        if strcmp(eobj.dataDomain{2},'SH')
            eobj=eobj.toSpace(Nbf,'MIC');
        end
        P=repmat(squeeze(eobj.data(1,1:ceil((nFFT+1)/2),:)),[1,1,L]);                            %[nFFT/2 x (N+1)^2 x L]
        W=bf_design_spatial(eobj.orderN,eobj.micGrid.r,eobj.micGrid.elevation,eobj.micGrid.azimuth,th_look,ph_look,kVec_bf,beamType_space,eobj.scatterer,Snn_inv); %[L,Q,nFFT/2]
%         W=W.';
%       W=bf_design_spatial(Na,r_array,th_array,ph_array,th_look,ph_look,k_vec,bf_flag,sphere_flag,Snn_inv)
    case 7
        % Transform to Space domain, if needed
        if strcmp(eobj.dataDomain{2},'SH')
            eobj=eobj.toSpace(Nbf,'MIC');
        end
        P=repmat(squeeze(eobj.data(1,1:ceil((nFFT+1)/2),:)),[1,1,L]);                            %[nFFT/2 x (N+1)^2 x L]
        W=bf_design_spatial(eobj.orderN,eobj.micGrid.r,eobj.micGrid.elevation,eobj.micGrid.azimuth,th_look,ph_look,kVec_bf,1,eobj.scatterer); %[L,Q,nFFT/2]
    otherwise
        %    performing bf in SH domain
        
        % Transform to SH domain, if needed
        if strcmp(eobj.dataDomain{2},'SPACE')
            eobj=eobj.toSH(Nbf,'MIC');
        elseif strcmp(eobj.dataDomain{2},'SH')
            if ~(eobj.orderN==Nbf)  % if Nbf and orderN are not the same, then fix it
                eobj=eobj.toSpace('MIC');
                eobj=eobj.toSH(Nbf,'MIC');
            end
        end
        
        %     W=bf_design_spatial(eobj.orderN,eobj.micGrid.r,eobj.micGrid.elevation,eobj.micGrid.azimuth,th_look,ph_look,k_vec(1:size(P,2)),beamType,eobj.scatterer); %[L,Q,nFFT/2]
        switch beamType
            case 4,
                for ang_ind=1:length(th_look),
                    W(:,:,ang_ind)=bf_design_SH(eobj.orderN,eobj.micGrid.r,th_look(ang_ind),ph_look(ang_ind),kVec_bf,beamType,eobj.scatterer,beamParam); %[nFFT/2 x (N+1)^2 x L]
                end
            case 5,
                %% beamParam:[{a_max},{a_m_max},{Nm},{D},{inv_DDH}]
                kVec=double(kVec);
                kr=double(kVec*eobj.micGrid.r(1));
                if length(beamParam)>2,
                    Nm=beamParam{3};
                else
                    Nm=ceil(max(kr));  %assuming there is no temporal aliasing this maximal order should be satisfied
                end                
                %calculates matrix 'D' and 'inv(DD^H)'
                [D3,DDH3]=getD(Nm,Nbf,kr,eobj.scatterer,eobj.micGrid.elevation,eobj.micGrid.azimuth,beamParam);
                beamParam=[beamParam,{D3},{DDH3}];
                fprintf('calculates W for angle ind:               ')
                 for ang_ind=1:length(th_look),
                    W(:,:,ang_ind)=bf_design_SH(eobj.orderN,eobj.micGrid.r,th_look(ang_ind),ph_look(ang_ind),kVec_bf,beamType,eobj.scatterer,beamParam); %[nFFT/2 x (N+1)^2 x L]
                    for print_ind=[0:(floor(log10(length(th_look))+1)+ceil(log10(ang_ind)))], fprintf('\b');  end;     fprintf('%d/%d', ang_ind,length(th_look)); % delete previous counter display and than display new
                end
               fprintf('\n');

            otherwise
                W=bf_design_SH(eobj.orderN,eobj.micGrid.r,th_look,ph_look,kVec_bf,beamType,eobj.scatterer,beamParam); %[nFFT/2 x (N+1)^2 x L]
                
        end
        P=repmat(squeeze(eobj.data(1,1:ceil((nFFT+1)/2),:)),[1,1,L]);                            %[nFFT/2 x (N+1)^2 x L]
        
end

clear D3 DDH3 beamParam
%% Zero DC and Fs/2 frequencies
W(1,:)=0;
W(end,:)=real(W(end,:));


%% checks theoretical and actual beampattern
if 0,
    doaPh1=linspace(0,2*pi,73);
    doaTh1=pi/2*ones(1,73);
    [doaP1(1,:),doaP1(2,:),doaP1(3,:)]=sph2cart(doaPh1,pi/2-doaTh1,1);
    doaP=getRotationMat(pi/2-ph_look(1),pi/2-th_look(1),ph_look(1)-pi/2)*doaP1;
    [doa1(1,:),doa1(2,:),r_tmp]=cart2sph(doaP(1,:),doaP(2,:),doaP(3,:));
    doa1(2,:)=pi/2-doa1(2,:);
    doa1=sortrows(doa1.',1).';
    doaTh=doa1(2,:);     doaPh=doa1(1,:); 
%     figure; plot(doaTh); hold on; plot(doaPh,'r'); legend([{'\theta_{doa}'};{'\phi_{doa}'}]); ylabel('doa [rad]'); xlabel('direction index');
    
    L=length(doaTh);
    kr=double(kVec_bf*eobj.micGrid.r(1));

    Wnm=bf_design_SH(eobj.orderN,eobj.micGrid.r,doaTh,doaPh,kVec_bf,beamType,eobj.scatterer,a_max); %[nFFT/2 x (N+1)^2 x L]
    Wnm3=permute(Wnm,[1,3,2]);

    
    % checks theoretical beampattern
    Ypw=shiftdim(sh2(eobj.orderN,th_look(1),ph_look(1)).',-1);
    Ypw3=repmat(Ypw,[length(kr),L,1]);
    bn = radialMatrix(eobj.orderN,kr,eobj.scatterer,inf,true);
    bn3=repmat(permute(bn,[1,3,2]),[1,L,1]);
    Pnm3=bn3.*conj(Ypw3);
    % preform filtering
    Y=sum(conj(Wnm3).*Pnm3,3);
%     [DoaPh,KR]=meshgrid(doaPh,kr);
%     figure; surf(DoaPh*180/pi,KR,abs(Y),'EdgeColor','none'); view(2); colorbar; axis([min(doaPh)*180/pi,max(doaPh)*180/pi,min(kr),max(kr)]); xlabel('doa [deg]'); ylabel('kr'); title('Theoretical beampattern');    
    Yavg=mean(abs(Y)); Yavg=Yavg/max(abs(Yavg));
    Y=Y./repmat(max(abs(Y),[],2),[1,size(Y,2)]);
    figure; l_ind=601; h_ind=4001; plot(doaPh(:)*180/pi,abs(Y(h_ind,:)));   
    hold on; plot(doaPh(:)*180/pi,abs(Y(l_ind,:)),'r'); plot(doaPh(:)*180/pi,abs(Yavg),'k'); xlim([min(doaPh),max(doaPh)]*180/pi); xlabel('doa [deg]'); title('Theoretical normalized beampattern'); legend([{['f=',num2str(round(fVec(h_ind))),' Hz']};{['f=',num2str(round(fVec(l_ind))),' Hz']};{['average beampattern']}]);
    figure; imagesc(doaPh*180/pi,kr,abs(Y));  xlabel('doa [deg]'); ylabel('kr'); title('Theoretical normalized beampattern');    

    % checks actual beampattern

    
    Pnm=P(:,:,1);
    Pnm3=repmat(permute(Pnm,[1,3,2]),[1,L,1]);
    Y=sum(conj(Wnm3).*Pnm3,3);
    Yavg=mean(abs(Y)); Yavg=Yavg/max(abs(Yavg));
    Y=Y./repmat(max(abs(Y),[],2),[1,size(Y,2)]);
    figure; l_ind=601; h_ind=4001; plot(doaPh(:)*180/pi,abs(Y(h_ind,:)));   
    hold on; plot(doaPh(:)*180/pi,abs(Y(l_ind,:)),'r'); plot(doaPh(:)*180/pi,abs(Yavg),'k'); xlim([min(doaPh),max(doaPh)]*180/pi); xlabel('doa [deg]'); title('Actual normalized beampattern'); legend([{['f=',num2str(round(fVec(h_ind))),' Hz']};{['f=',num2str(round(fVec(l_ind))),' Hz']};{['average beampattern']}]);
    figure; imagesc(doaPh*180/pi,kr(10:end),abs(Y(10:end,:)));  xlabel('doa [deg]'); ylabel('kr'); title('Actual normalized beampattern');    
    
end
%%
% P=shiftdim(mic_sig_f.data(1,1:ceil((nFFT+1)/2),:),1).'; %take the single side fft - for even or odd nData
% P(:,1) = 0; % put zero DC (f=0).  P:[Q,~nFFT/2]

%% Checks the beamformer weights magnitude
if 0,
    figure; plot(fVec(1:size(W,1)).',20*log10(abs(W(:,:,1)))); ylim([-65,15]);
    figure; plot(fVec(1:size(W,1)).',(angle(W(:,1,1)))); ylim([0,2]); xlim([0,4000]);
    
end
%% Checks the impulse response of filter for each microphone
if 1,
    WW=cat(1,W, conj(flipdim(W(2:(nFFT_bf-size(W,1)+1),:,:),1)) );
    W_t=ifft(WW,[],1);
    s_T=1/eobj.fs;
    s_N=nFFT_bf;
    t_vec=[0:(s_N-1)].'*s_T;
    figure; plot(t_vec,20*log10(abs(W_t(:,:,1)))); xlabel(' t [sec]'); ylabel('h_{mic} [dB]');

    pw_win_N=256;
    [pw_ps,pw_f_vec]=pwelch(W_t(:,1,1)/max(abs(W_t(:,1,1))),pw_win_N,pw_win_N/2,pw_win_N,eobj.fs);
    [pw_ps2,pw_f_vec]=pwelch(W_t(:,end,1)/max(abs(W_t(:,end,1))),pw_win_N,pw_win_N/2,pw_win_N,eobj.fs);
     figure; plot(pw_f_vec,10*log10(pw_ps)); hold on; plot(pw_f_vec,10*log10(pw_ps2),'r');  xlabel('f [Hz]'); ylabel('Power/frequency (dB/Hz/sample)'); 
     if beamType<=5,
         legend([{'order (n,m)=(0,0) '};{'order (n,m)=(3,3)'}]);
     else
         legend([{'mic 1'};{'mic32'}]);
     end
     
     clear W_t
     
     %% Implement here filtering in time
    PP=cat(1,P, conj(flipdim(P(2:(nFFT-size(P,1)+1),:,:),1)) );
     Wt_Cnj=ifft(conj(WW),[],1);
     W_t_cs=circshift(Wt_Cnj,round(size(Wt_Cnj,1)/2));
%     W_t_cs=W_t;%circshift(W_t,round(length(W_t)/2));
         figure; plot(t_vec,20*log10(abs(W_t_cs(:,:,1)))); xlabel(' t [sec]'); ylabel('h_{mic} [dB]');
conv_length=size(W_t_cs,1)+size(PP,1)-1;
         
YY=squeeze(sum(fft(W_t_cs,conv_length).*[fft(ifft(PP),conv_length)],2));%single side fft of array output Y:[ nFFT/2 X L]
% YY(ceil(size(YY,1)/2)+1,:)=real(YY(ceil(size(YY,1)/2)+1,:));
% YY=[Y;conj(flipud(Y(2:(nFFT-size(Y,1)+1),:)))]; %YY size [ nFFT x L]
% y=ifft(YY,size(eobj.data,2)); %size [nTime x L]
y=ifft(YY); %size [nTime x L]
% for mic_i=1:32,
% y1(:,mic_i)=conv(W_t(:,mic_i),ifft(PP(:,mic_i)));
% end
% y1=sum(y1,2);
% for mic_i=1:32,
% y1_cs(:,mic_i)=conv(W_t_cs(:,mic_i),ifft(PP(:,mic_i)));
% end
% y1_cs=sum(y1_cs,2);
% soundsc(y1,16000);
% soundsc(y1_cs,16000);
% % figure; plot(real(y));
 y=y(round(size(Wt_Cnj,1)/4)+[1:nFFT],:);
% soundsc(y_cs(8000:40000),16000);
% soundsc(y(16000:end),16000);

end
%% preform spatial filtering
% % appliyng beamfomer weights to calculate array output
% % % Y=squeeze(sum(conj(W).*P,2));%single side fft of array output Y:[ nFFT/2 X L]
% % % YY=[Y;conj(flipud(Y(2:(nFFT-size(Y,1)+1),:)))]; %YY size [ nFFT x L]
% % % y=ifft(YY,size(eobj.data,2)); %size [nTime x L]
% if w_ir_length>1, y=y(1:end-(w_ir_length-1),:); end
% 

%%%% removes zero padding


