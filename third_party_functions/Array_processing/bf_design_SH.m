function W_nm=bf_design_SH(Na,r_array,th_look,ph_look,kVec,beamType,sphereType,beamParam)
%% ======   Function Description   ======
% This function generates a beamfer weights for multicahnnel (Q-channels) data from a microphone array. 
% The beamforming weights are calculated for different frequencies in spatial domain. Different types of beamformers can be design (currently only maxDI).
%
% version 1.0 - created on 23/07/14 by David Alon 
%
%===== Syntax  ==================
%   W=bf_design_SH(Na,r_array,th_array,ph_array,th_look,ph_look,kVec,bf_flag,sphere_flag)
%
%===== Input variables  =========
%   Na- microphone array SH order, range [0:Inf]   
%%%   r_array- microphone array radius [m], size [1], range [0:Inf]  
%   th_look- beamformer look direction elevation [rad] , range [0:pi] , default value={pi/2}
%   ph_look- beamformer look direction azimuth [rad] , range [0:2*pi] , default value={0}
%%%   kVec- wave number vector in which the beamformer is calculated [rad/m]
%   beamType-  beamformer type indicator: (1) max DI ;(2) max WNG ; (3) Peled optimal weights  ; (4) SH domain MVDR. ; (5) aliasing cancellation
%   sphereType- sphere scattering indicator: (1) rigid ;(0) open
%
%              a_max=beamParam{1}; sigmaA=beamParam{2}; sigmaI=beamParam{3};
%   beamParam - struct with the following fields:
%                 a_max - soft limiting amplification limit [dB]
%                 sigmaA - acoustic noise std
%                 sigmaI - sensor noise std
%                 Snn_inv - sensor noise autocorrelation matrix
%                 D - aliasing matrix
%                 invDDH - inverse aliasing matrix 
%                  Nm - sound field model order - for aliaisng cancellation 
% beamParam: case 1:a_max ,  case 3: [sigmaA,sigmaI] ,  case 4: Snn_inv, case 5: [{D},{invDDH}];

%  !!!! note: in case of beamType =3 (Peled weights) sigma_I is assumed to be divided by Q, the numner of microphones. sigma_I_original=sigma_I_here*Q
                
%===== Output variables  ========
%   W - matrix containing the beamformer weights (in spatial domain) over differnet frequencies. size [Q,length(kVec)]
%



  useSOFiA = false;                     % change to true if you wish to generate the radial functions using SOFiA toolbox
                                        % instead of EARO's internal func.
    

  if nargin<8, beamParam=struct('a_max',inf); end;          % Default = no soft-limiting
  if nargin<7, sphereType=1; end;           % Default = Rigid sphere
  if nargin<6, beamType=1; end;             % Default = max DI
  
  % Construct kr vector
  kVec=double(kVec);
  kr=double(kVec*r_array(1));
  
%   bn = radialMatrix(Na,kr,sphereType,inf,true,useSOFiA);
%   bn=fix_bn(bn,1,0.2,kr);

            Nfft=length(kVec);
            L=length(th_look);

  % Generate beamforming coefficients
    switch beamType
        case 1 %max DI
               bn = radialMatrix(Na,kr,sphereType,beamParam.a_max,true,useSOFiA);   % Generate radial functions
%  SL_power=1.5; limitParam=[a_max,SL_power];
% bn = radialMatrix2(Na,kr,sphereType,limitParam,true,useSOFiA);   % Generate radial functions
             % figure;
%              bn = radialMatrix2(Na+1,kr,sphereType,a_max,true,useSOFiA); figure; plot(kr,20*log10(abs(bn)));  ylabel('|b_n(f)| [dB]'); xlabel('kr'); axis([0,5,-10,30]);
            Y_look=shMatrix(Na,th_look,ph_look);
            dnm=permute(repmat(Y_look,[1,1,Nfft]),[3,1,2]);
        case 2 %Delay and Sum
             bn = radialMatrix(Na,kr,sphereType,inf,true,useSOFiA);   % Generate radial functions
             Y_look=shMatrix(Na,th_look,ph_look);
             dnm=permute(repmat(Y_look,[1,1,Nfft]),[3,1,2]).*repmat(abs(bn).^2,[1,1,L]);
        case 3 %optimal weights [Peled]: trade-off between acoustic noise and sensor noise
%              a_max=beamParam{1}; sigmaA=beamParam{2}; sigmaI=beamParam{3};
             bn = radialMatrix(Na,kr,sphereType,beamParam.a_max,true,useSOFiA);   % Generate radial functions
%              bn(1,:)=bn(2,:);
             bnDiag=bn(:,([0:Na]+1).^2);  %???Dimensions
             Y_look=shMatrix(Na,th_look,ph_look);
             b=(1:2:2*Na+1).'/(4*pi);
             %constracting matrix T   (T is multiplied by 'J' comparing to the paper version)
             T=zeros(Nfft,(Na+1),(Na+1));
             [Ind1,Ind2]=meshgrid([1:Nfft],[1:(Na+1)]);
             T(sub2ind(size(T),Ind1,Ind2,Ind2).')=repmat(4*pi*b.',Nfft,1)./abs(bnDiag).^2;
%              T=diag(b./(abs(bnDiag.').^2));  %???Dimensions
             %constracting matrix P   
             P=repmat(shiftdim(diag(b)/(4*pi),-1),[Nfft,1,1]);  % multiplied by 4*pi because of the definition of 'b' ["Objective performance analysis..." paper] 
             %constracting matrix R^-1
             Rinv=(beamParam.sigmaA^2*P+beamParam.sigmaI^2*T);
             Rinv(Rinv~=0)=Rinv(Rinv~=0).^-1;
             Rinv_b=sum(Rinv.*repmat(b.',[Nfft,1,Na+1]),3);
             bH_Rinv_b=sum(Rinv_b.*repmat(b',[Nfft,1]),2);
             % calculating optimal (Na+1) coefficients
             dn=Rinv_b./(repmat(bH_Rinv_b,1,(Na+1)));
             % ca
             dn2=zeros(Nfft,(Na+1)^2);
             indMat=zeros((Na+1)^2,Na+1);
             for ind1=0:Na,
                 indMat(((ind1)^2)+[1:(2*ind1+1)],ind1+1)=1;
             end
             dn2=(indMat*dn.').';
             dnm=permute(repmat(Y_look,[1,1,Nfft]),[3,1,2]).*repmat(dn2,[1,1,L]);  
               figure; plot(kr,20*log10(abs(bn)));  ylabel('|b_n(f)| [dB]'); xlabel('kr'); axis([0,5,-10,30]);  grid on;
               figure; plot(kr,20*log10(abs(dn)));  ylabel('|d_n(f)| [dB]'); xlabel('kr'); axis([0,5,-10,30]);  grid on;
                legend_str=[{'d_0(kr)'};{'d_1(kr)'};{'d_2(kr)'};{'d_3(kr)'};{'d_4(kr)'};{'d_5(kr)'}];
               legend(legend_str(1:(Na+1)));
               figure; plot(kr,20*log10(abs(bn(:,[1,3,7])./dn)));  ylabel('|1/g_n(f)| [dB]'); xlabel('kr'); axis([0,5,-30,30]); grid on;
                legend_str=[{'1/g_0(kr)'};{'1/g_1(kr)'};{'1/g_2(kr)'}];
               legend(legend_str(1:(Na+1)));
        case 4
            
            Snn_inv=permute(beamParam.Snn_inv(:,:,1:Nfft),[3,1,2]);
        
        bn = radialMatrix(Na,kr,sphereType,beamParam.a_max,true,false);   % Generate radial functions
        Br_array3=repmat(permute(bn,[1,3,2]),[1,(Na+1)^2,1]); %size [Nfft x (Na+1)^2 x (N+1)^2]
        Y_look=shMatrix(Na,th_look(1),ph_look(1));
        Y_l_conj3=repmat(permute(conj(Y_look(:,1)),[2,3,1]),[Nfft,(Na+1)^2,1]); %size [Qx(N+1)^2xNfft]
        V_l=sum(Br_array3.*Y_l_conj3,3);
        V_l3=repmat(permute(V_l,[1,3,2]),[1,(Na+1)^2,1]);
        
        SV=sum(Snn_inv.*V_l3,3);  %summing and multiplying over the rows of Snn and rows of V_l
        SVS1=sum(conj(V_l).*SV,2);
        W_nm=repmat(SVS1.^-1,[1,(Na+1)^2]).*SV;
           
%             a_max=beamParam{1};
%             bn = radialMatrix(Na,kr,sphereType,a_max,true,useSOFiA);   % Generate radial functions
%             Y_look=shMatrix(Na,th_look(1),ph_look(1));
%             
%             Y_l_conj2=repmat(Y_look',[Nfft,1]); %size [???]
%             V_l2=bn.*Y_l_conj2;
%             SV2=squeeze(sum(Snn_inv.*repmat(V_l2,[1,1,(Na+1)^2]),2));
%             SVS1=sum(conj(V_l2).*SV2,2);
%             W_nm=squeeze(repmat(SVS1.^-1,[1,(Na+1)^2]).*SV2);
%             
        case 5
            %calculates dnm
%             a_max=beamParam{1};
            Nm=beamParam{3};
            D=beamParam{4};
            invDDH=beamParam{5};
            clear beamParam;
            SL_power=1.5; limitParam=[beamParam.a_max,SL_power];
            bn = radialMatrix2(Na,kr,sphereType,limitParam,true,useSOFiA);   % Generate radial functions
%             bn = radialMatrix(Na,kr,sphereType,a_max,true,useSOFiA);   % Generate radial functions
            for ang_ind=1:length(th_look),
                Y_look=shMatrix(beamParam.Nm,th_look(ang_ind),ph_look(ang_ind));
                YDH=sum(repmat(shiftdim(Y_look,-2),[Nfft,(Na+1)^2,1]).*conj(beamParam.D),3);
                dnm_tmp=4*pi*sum(repmat(YDH,[1,1,(Na+1)^2]).*beamParam.invDDH,2);
                dnm(:,:,ang_ind)=permute(dnm_tmp,[1,3,2]);
            end

            
            
    end
    
    if beamType~=4
        % Generate beamforming weights
        Wnm_conj=dnm./repmat(bn,[1,1,L]);                
        % Output
        W_nm = conj(Wnm_conj);
    end
  

end



