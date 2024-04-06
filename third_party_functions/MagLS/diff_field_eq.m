function [H_l_nm_low,H_r_nm_low] = diff_field_eq(Y_high,H_high,H_l_nm_low,H_r_nm_low,cutOffFreqInd_high)
    H_l = squeeze(double(H_high(:,:,1))).'; % left HRTF
    H_r = squeeze(double(H_high(:,:,2))).'; % right HRTF
    pY_high = pinv(Y_high); % pseudo-inverse of matrix Y_high
    H_l_nm_high = H_l*pY_high; % SH transform
    H_r_nm_high = H_r*pY_high; % SH transform
    
%    for freqInd = 1 : cutOffFreqInd_high
    for freqInd = 1 : size(H_r_nm_high,1)
        H_SH_high = [H_l_nm_high(freqInd,:); H_r_nm_high(freqInd,:)].';
        H_SH = [H_l_nm_low(freqInd,:); H_r_nm_low(freqInd,:)].';

        % Eq. (4.63)
        X = chol(H_SH_high' * H_SH_high);
        X_est = chol(H_SH' * H_SH);

        [U,~,V] = svd(X_est'*X);

        H_SH_corr = ( H_SH * inv(X_est) * V * U' * X).';

        H_l_nm_low(freqInd,:) = H_SH_corr(1,:);
        H_r_nm_low(freqInd,:) = H_SH_corr(2,:);
    end
end