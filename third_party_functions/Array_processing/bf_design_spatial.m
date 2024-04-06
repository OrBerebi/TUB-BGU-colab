function W=bf_design_spatial(Na,r_array,th_array,ph_array,th_look,ph_look,k_vec,bf_flag,sphere_flag,Snn_inv)
%% ======   Function Description   ======
% This function generates a beamformer weights for multicahnnel (Q-channels) data from a microphone array.
% The beamforming weights are calculated for different frequencies over L look directions in spatial domain. 
% Different types of beamformers can be design (currently only maxDI).
%
% version 1.0 - created on 23/07/14 by David Alon
%             - edited on 22/11/14 by Jonathan Sheaffer
%
%===== Syntax  ==================
%   W=bf_design_spatial(Na,r_array,th_array,ph_array,th_look,ph_look,k_vec,bf_flag,sphere_flag)
%
%===== Input variables  =========
%   Na- microphone array SH order, range [0:Inf]
%   r_array- array microphones radii [m], size [1xQ], range [0:Inf]
%   th_array- array microphones elevation [rad], size [1xQ], range [0:pi]
%   ph_array- array microphones azimuth [rad] , size [1xQ], range [0:2*pi]
%   th_look- beamformer look direction elevation [rad], size [1xL], range [0:pi] , default value={pi/2}
%   ph_look- beamformer look direction azimuth [rad], size [1xL], range [0:2*pi] , default value={0}
%   k_vec- wave number vector in which the beamformer is calculated [rad/m]
%   bf_flag - beamformer type indicator: 
%                            1 = Maximum Directivity
%                            2 = Maximum White Noise Gain (WNG)
%                            3 = Maximum Variance Distortionless Response
%                            (MVDR)
%                            4 = MVDR / type II
%
%   sphere_flag- sphere scattering indicator: (1) rigid ;(0) open
%
%===== Output variables  ========
%   W - matrix containing the beamformer weights (in spatial domain) over differnet frequencies. size [LxQxlength(k_vec)]
%

Nfft=length(k_vec);  %number of frequencies
Q=length(th_array);  % number of microphones
L=length(th_look);  % number of microphones

bn_thresh=2*10^-1;
bn_handling=0; %flag inidicates how to overcome the low bessel function problem:
%(1)-clliping low bessel function values; (2)-loading bessel function values; (3)- removing low bessel function orders
%Note: for the case of array with different radii method (3) is not defined
switch bf_flag
    case 1 %Max_DI
        Y_array=sh2(Na,th_array,ph_array).';
        Y_look=sh2(Na,th_look,ph_look);
        Br_array=zeros(length(k_vec),(Na+1)^2);
        if length(unique(r_array))==1,
            disp('Generating beamforming weights - single sphere');
            % Y_array=sh2(mic_s.orderN,mic_s.micGrid{elevation},mic_s.micGrid{azimuth}).';
            Ypsd=(Y_array'*Y_array)\Y_array';
            % Pnm=Ypsd*P;%(:,1);  %Pnm:[(orderN+1)^2,N_fft/2]
            if sphere_flag==0,
                display('warning: the b_n calculation for open sphere may contain errors');
                Br_array(2:end,:)=B2(Na,k_vec(2:end)*r_array(1),k_vec(2:end)*r_array(1),sphere_flag);
            else
%                Br_array(2:end,:)=b_n(Na,k_vec(2:end)*r_array(1),1,1).';
%                  Br_array(2:end,:)=b_n(Na,k_vec(2:end)*r_array(1),-1,-1).';
                %                 Br_array(2:end,:)=b_n(Na,k_vec(2:end)*r_array(1),-1,1).';
      timeDep = true;
      a_max = inf;
                 Br_array(2:end,:) = radialMatrix (Na,k_vec(2:end)*r_array(1),sphere_flag,a_max,timeDep);
                1;
            end
            Br_array(1,:)=Br_array(2,:);
            Br_array=fix_bn(Br_array,bn_handling,bn_thresh,k_vec.'*r_array(1));
            Br_inv=zeros(size(Br_array,1),size(Br_array,2));
            Br_inv((Br_array~=0))=Br_array(Br_array~=0).^(-1);
%             W_hermit=4*pi/(Na+1)^2*(repmat((Y_look.'),length(k_vec),1).*Br_inv)*Ypsd;
            Y_look4=repmat(Y_look.',[1,1,Q,Nfft]);
            Br_inv4(1,:,1,:)=Br_inv.';
            Br_inv4=repmat(Br_inv4(1,:,1,:),[L,1,Q,1]);
            Y_psd4(1,:,:,1)=Ypsd;
            Y_psd4=repmat(Y_psd4,[L,1,1,Nfft]);
            W_hermit3=4*pi/(Na+1)^2*sum(Y_look4.*Br_inv4.*Y_psd4,2);
            W_hermit3=(permute(W_hermit3,[4,3,1,2]));
            W=conj(W_hermit3);
        else
            if bn_handling==3, bn_handling=1; disp('bn_handling=3 method is not defined for diferent radii array'); end
            disp('Generating beamforming weights - varaible radii');
            fprintf('frequency:          ');
            W=zeros(length(th_array),length(k_vec));
            for f_ind=2:length(k_vec)
                Br_array=B2(Na,k_vec(f_ind)*r_array,k_vec(f_ind)*r_array,sphere_flag);
                Br_array=fix_bn(Br_array,bn_handling,bn_thresh,k_vec*r_array(1));
                BY_a=Br_array.*Y_array;
                C=(BY_a'*BY_a)\BY_a';
                W(:,f_ind)=4*pi/(Na+1)^2*C'*conj(Y_look);
                if (mod(f_ind,100)==0), for print_ind=[0:(floor(log10(f_ind-1))+2+floor(log10(length(k_vec))))], fprintf('\b');  end;     fprintf('%d/%d',f_ind,(length(k_vec))); % delete previous counter display and than display new
                end
            end
            fprintf('\n');
        end
    case 3 %MVDR 08/2014
        
        Y_array=shMatrix(Na,th_array,ph_array).';
        Y_look=shMatrix(Na,th_look,ph_look);
        Br_array=zeros(length(k_vec),(Na+1)^2);
        %         tic
        %             if bn_handling==3, bn_handling=1; disp('bn_handling=3 method is not defined for diferent radii array'); end
        %             disp('Generating MVDR beamforming weights');
        %             fprintf('frequency:          ');
        %             W=zeros(Q,Nfft);
        %             for f_ind=2:length(k_vec)
        %                 Br_array=B2(Na,k_vec(f_ind)*r_array,k_vec(f_ind)*r_array,sphere_flag);
        %                 Br_array=fix_bn(Br_array,bn_handling,bn_thresh,k_vec*r_array(1));
        %                 BY_a=Br_array.*Y_array;
        %                 V=BY_a*conj(Y_look);
        %                 W(:,f_ind)=Snn_inv(:,:,f_ind)*V/(V'*Snn_inv(:,:,f_ind)*V);
        %                 if (mod(f_ind,100)==0), for print_ind=[0:(floor(log10(f_ind-1))+2+floor(log10(length(k_vec))))], fprintf('\b');  end;     fprintf('%d/%d',f_ind,(length(k_vec))); % delete previous counter display and than display new
        %                 end
        %             end
        %             fprintf('\n');
        %             toc
        %             W_tmp=W;
        tic
%         Br_array=zeros(length(k_vec),(Na+1)^2);
%         Br_array3=zeros(1,(Na+1)^2,Nfft);
%         Br_array(2:end,:)=B2(Na,k_vec(2:end)*r_array(1),k_vec(2:end)*r_array(1),sphere_flag);
%         Br_array(1,:)=Br_array(2,:);
%         Br_array2=fix_bn(Br_array,bn_handling,bn_thresh,k_vec*r_array(1)).';
%         Br_array3(1,:,:)=Br_array2;
                bn = radialMatrix(Na,k_vec*r_array(1),sphere_flag,8,true,false);   % Generate radial functions
        Br_array3=permute(bn,[3,2,1]);
        Br_array3=repmat(Br_array3,[Q,1,1]); %size [Qx(N+1)^2xNfft]
        Y_array3=repmat(Y_array,[1,1,Nfft]);
        BY3=(Br_array3.*Y_array3);
        Y_l_conj3=repmat(Y_look(:,1)',[Q,1,Nfft]); %size [Qx(N+1)^2xNfft]
        V_l=sum(BY3.*Y_l_conj3,2);
        V_l_trans(1,:,:)=V_l;
        SV2=sum(Snn_inv.*repmat(V_l_trans,[Q,1,1]),2);
        SVS1=sum(conj(V_l).*SV2,1);
        W=(squeeze(repmat(SVS1.^-1,[Q,1,1]).*SV2)).';
        toc
    case 5 %MVDR 09/2014
        Snn_inv=permute(Snn_inv,[3,1,2]);
        
        %          Y_l_conj2=repmat(Y_look',[Nfft,1]); %size [???]
        bn = radialMatrix(Na,k_vec*r_array(1),sphere_flag,8,true,false);   % Generate radial functions
        Br_array3=repmat(permute(bn,[1,3,2]),[1,Q,1]); %size [Nfft x Q x (N+1)^2]
        clear bn;
        Y_array=shMatrix(Na,th_array,ph_array).';
        Y_array3=repmat(permute(Y_array,[3,1,2]),[Nfft,1,1]);
        Y_look=shMatrix(Na,th_look(1),ph_look(1));
        Y_l_conj3=repmat(permute(conj(Y_look(:,1)),[2,3,1]),[Nfft,Q,1]); %size [Qx(N+1)^2xNfft]
        V_l=sum(Br_array3.*Y_array3.*Y_l_conj3,3);
%         V_l_conj3=repmat(permute(conj(V_l),[1,3,2]),[1,Q,1]);
        V_l3=repmat(permute(V_l,[1,3,2]),[1,Q,1]);
        clear Br_array3 Y_array3 Y_l_conj3;
%         V_l2=bn.*Y_l_conj2;        
        SV=sum(Snn_inv.*V_l3,3);  %summing and multiplying over the rows of Snn and rows of V_l
        SVS1=sum(conj(V_l).*SV,2);
        W=repmat(SVS1.^-1,[1,Q]).*SV;
%         W=conj(W_conj);

        
        
%         Y_array=sh2(Na,th_array,ph_array).';
% %             Y_look=shMatrix(Na,th_look(1),ph_look(1));
%         Y_look=sh2(Na,th_look,ph_look);
% %         Br_array=zeros(length(k_vec),(Na+1)^2);
% %         %         tic
%         %             if bn_handling==3, bn_handling=1; disp('bn_handling=3 method is not defined for diferent radii array'); end
%         %             disp('Generating MVDR beamforming weights');
%         %             fprintf('frequency:          ');
%         %             W=zeros(Q,Nfft);
%         %             for f_ind=2:length(k_vec)
%         %                 Br_array=B2(Na,k_vec(f_ind)*r_array,k_vec(f_ind)*r_array,sphere_flag);
%         %                 Br_array=fix_bn(Br_array,bn_handling,bn_thresh,k_vec*r_array(1));
%         %                 BY_a=Br_array.*Y_array;
%         %                 V=BY_a*conj(Y_look);
%         %                 W(:,f_ind)=Snn_inv(:,:,f_ind)*V/(V'*Snn_inv(:,:,f_ind)*V);
%         %                 if (mod(f_ind,100)==0), for print_ind=[0:(floor(log10(f_ind-1))+2+floor(log10(length(k_vec))))], fprintf('\b');  end;     fprintf('%d/%d',f_ind,(length(k_vec))); % delete previous counter display and than display new
%         %                 end
%         %             end
%         %             fprintf('\n');
%         %             toc
%         %             W_tmp=W;
%         tic
% %         Br_array=zeros(length(k_vec),(Na+1)^2);
% %         Br_array3=zeros(1,(Na+1)^2,Nfft);
% %         Br_array(2:end,:)=B2(Na,k_vec(2:end)*r_array(1),k_vec(2:end)*r_array(1),sphere_flag);
% %         Br_array(1,:)=Br_array(2,:);
% %         Br_array2=fix_bn(Br_array,bn_handling,bn_thresh,k_vec*r_array(1)).';
% 
%         Y_array3=repmat(Y_array,[1,1,Nfft]);
%         BY3=(Br_array3.*Y_array3);
%         Y_l_conj3=repmat(Y_look(:,1)',[Q,1,Nfft]); %size [Qx(N+1)^2xNfft]
%         V_l=sum(BY3.*Y_l_conj3,2);
%         V_l_trans(1,:,:)=V_l;
%         SV2=sum(Snn_inv.*repmat(V_l_trans,[Q,1,1]),2);
%         SVS1=sum(conj(V_l).*SV2,1);
%         W=squeeze(repmat(SVS1.^-1,[Q,1,1]).*SV2);
%         toc
        
end

