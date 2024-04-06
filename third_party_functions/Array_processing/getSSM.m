function Snn_inv=getSSM(mic_n_ref_path,Snn_loading,Snn_Nfft,fs_analysis,SnnDomain)
%generates the inverse of the spatial spectral matrix or the cross-spectrum matrix in SH domain
% SnnDomain=1 [default] for Space domain Snn ; 2 for SH domain Snn

% Snn_inv(:,:,f_ind) is a result of E{p*p^H} where p is the zero mean pressure measured by the array at freq f_ind

if nargin<5
    SnnDomain=1;
end
Snn_var_loading=5*logspace(0,-2,ceil((Snn_Nfft+1)/2));

if SnnDomain==1
    Snn_file_path=['MVDR_out\Snn_inv_load_',num2str(Snn_loading),'Nfft_',num2str(Snn_Nfft),'.mat'];
else
    Snn_file_path=['MVDR_SH_out\Snn_inv_load_',num2str(Snn_loading),'Nfft_',num2str(Snn_Nfft),'.mat'];
end
% mic_n_ref_path=[root_dir,exp_path,'\',exp_path,'_',source_label{source_noise_ind},'\mic_s.mat'];
%size(mic_n.data,2)
% if (size(mic_n_ref.data,2))>(Snn_avg_num*Snn_Nfft)*(1-overlap)+1*(overlap) %checks that there is enough data for estimating the autocorrelation matrix
% end

if ~exist(Snn_file_path,'file')
    if 1,%exist(['MVDR_out\Snn_inv_load_0.mat'],'file')
        
        %loading reference noise record
        disp('loading mic_n from file ');
        load(mic_n_ref_path);
        mic_n_ref=mic_s;
        clear mic_s
        if mic_n_ref.fs~=fs_analysis,
        mic_n_ref=mic_n_ref.resampleData(fs_analysis);
        end
        % Transform to Time domain, if needed
        if ~strcmp(mic_n_ref.dataDomain{1},'TIME')
            mic_n_ref=mic_n_ref.toTime;
        end
        
        % Transform to SH domain, if needed
        if strcmp(mic_n_ref.dataDomain{2},'SPACE') && SnnDomain==2,
            mic_n_ref=mic_n_ref.toSH(mic_n_ref.orderN,'MIC');
        end
        
        
        p_t=squeeze(mic_n_ref.data).';
        disp('compute STFT - for p_t');
        switch SnnDomain,
            case 1, fprintf('Mic number:  ');
            case 2, fprintf('SH order number:  ');
        end
        %I = length( P(:,1) );
        
        
        Snn_avg_num=32;
        overlap=0.5;
        window_len=Snn_Nfft;
        overlap_len=overlap*window_len;
        
        clear mic_n_ref
        
        
        for mic_index = 1:size(p_t,1)   % P_stft : (mic x f x n)
            %     [P_stft(mic_index,:,:) num_freq num_frames] = STFT(p_t(mic_index,:).',window,overlap,frame_length,M);
            [P_stft_tmp(mic_index,:,:) stft_freq stft_frames]=spectrogram(p_t(mic_index,:),window_len,overlap_len,Snn_Nfft,fs_analysis);
            %     [P_stft(mic_index,:,:) num_freq num_frames] = spectrogram(p_t(mic_index,:).',window,overlap,frame_length,M);
            %-----STFT parameters-----
            %Hanning window of length T = 256 samples (16 msec), 75% overlap between adjacent frames , and 256 point FFT length
            %--STFT output-----------
            % for a 20 sec signal
            
            for print_ind=[0:floor(log10(mic_index))], fprintf('\b');  end;     fprintf('%d', mic_index); % delete previous counter display and than display new
            %[P_stft{mic_index} num_freq num_frames] = STFT(p_t(mic_index,:).',window,overlap,frame_length,M);
        end
        fprintf('\n');
        switch SnnDomain,
            case 1,         num_freq=length(stft_freq);
            case 2,         num_freq=ceil((length(stft_freq)+1)/2);
        end
        
        num_frames=length(stft_frames);
        disp('compute Snn from STFT ');
        fprintf('freq=         ');
        
        P_stft=P_stft_tmp;
        P_stft_mean=mean(P_stft,3);
        clear P_stft_tmp
        
        Snn_inv=zeros(size(P_stft,1),size(P_stft,1),num_freq);
        for f_ind=2:num_freq
            for Snn_ind=1:num_frames
                Snn_tmp(:,:,Snn_ind)=(P_stft(:,f_ind,Snn_ind)-P_stft_mean(:,f_ind))*(P_stft(:,f_ind,Snn_ind)-P_stft_mean(:,f_ind))';
                %         Snn_tmp(:,:,Snn_ind)=(P_stft(:,f_ind,Snn_ind))*(P_stft(:,f_ind,Snn_ind))';
            end
            Snn_tmp1=mean(Snn_tmp,3);
            [U,S,V] = svd(Snn_tmp1);
            if Snn_loading>=0 %positive value for fixed diagonal loading
                Snn_tmp2=U*(S+eye((size(p_t,1)))*Snn_loading)*V';
            else  %for variable loading
                Snn_tmp2=U*(S+eye((size(p_t,1)))*Snn_var_loading(f_ind))*V';
            end
            Snn_inv(:,:,f_ind)=inv(Snn_tmp2);
            %     Snn_inv(:,:,f_ind)=V*diag(diag((S+eye((size(p_t,1)))*Snn_loading)).^-1)*U';
            Snn_eigen(:,f_ind)=diag(S);
            %     Snn_min(f_ind)=min(min(abs(diag(Snn(:,:,f_ind)))));
            %     Snn_max(f_ind)=max(max(abs(diag(Snn(:,:,f_ind)))));
            Snn_cond(f_ind)=cond(Snn_tmp2);
            for print_ind=[0:(floor(log10(f_ind-1))+floor(log10(num_freq))+2)],     fprintf('\b');  end;     fprintf('%d/%d', f_ind,num_freq); % delete previous counter display and than display new
        end
        Snn_inv(:,:,1)=Snn_inv(:,:,2);
        fprintf('\n');
        clear P_stft;
        clear P_stft_mean;
        % figure;
        % subplot(2,1,1);
        % plot(stft_freq,10*log10(Snn_min).'); hold on; plot(stft_freq,10*log10(Snn_max).','r');  %plot(Snn_max,'g');
        % subplot(2,1,2);
        figure; plot(stft_freq(1:num_freq),10*log10(Snn_cond).'); title('Cond of Snn'); xlabel('f [Hz]');
        % imagesc(abs(Snn(:,:,1000))); colorbar;
        
        save(Snn_file_path,'Snn_inv','Snn_cond','Snn_eigen');
        
        figure;
        plot(10*log10(Snn_eigen(:,1500))); title(['Snn eigen values @ f=',num2str(round(stft_freq(1500))),'Hz']); ylabel('eigen value [dB]'); xlabel('index'); grid on;
        
        figure;
        plot(stft_freq(1:num_freq),10*log10(Snn_eigen(:,:))); ylim([-30,20]); ylabel('eigen value [dB]'); xlabel('freq [Hz]');
        % Snn=getSSM(mic_n_ref,size(mic_n.data,2));
    end
else
    disp('loading Snn from file ');
    load(Snn_file_path);
end
