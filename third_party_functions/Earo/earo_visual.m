function earo_visual(eobj, anm, f, Nvis)
% function earo_visual(eobj, anm, f, [Nvis])
%
% Plots a SH function in spacc, e.g. plane wave decomposition coefficients
%
%    eobj       - Input EARO object
%    anm        - Function to plot, in the SH domain. This can be obtained
%                   e.g. by using the pwd_earo function.
%    f          - Frequency [Hz] to plot
%    [Nvis]     - (optional) ploting resolution in the space domain. if
%                   unspecified then this is computed from the order of anm
%
% August 2014, Jonathan Sheaffer, Ben-Gurion University
% Part of the EARS Beamforming toolbox
%

    if nargin<4, Nvis=floor(sqrt(size(anm,2))-1); end;  % Default is order of anm
    
    % find nearest frequency
    nFFT=size(eobj.data,2);
    fVec = (0:(nFFT-1))*(eobj.fs/nFFT);

    freqIdx=findnearest(f,fVec);
    if (fVec(freqIdx)~=f) && (eobj.shutUp==false)
        fprintf('Chosen frequency %d Hz does not exist in data.\nUsing nearest frequency, %d Hz, instead.\n',f,fVec(freqIdx));
    end
    
    % Generate the visualization grid
    [tmp,th,ph]=equiangle_sampling(Nvis);
    
    % Interpolate spatial data
    Y=shMatrix(sqrt(size(anm,2))-1,th,ph);
    ak=anm(freqIdx,:)*Y;
    phi=ph(1:2*(Nvis+1));
    the=th(1:2*(Nvis+1):end);

    % Plot spatial data
    w2=(reshape(ak,2*(Nvis+1),2*(Nvis+1)).');
    [c,h]=contourf(phi*180/pi,the*180/pi,abs(w2));
    grid on;

    set(gca,'FontSize',14);
    colorbar;
    set(colorbar,'FontSize',14);
    %colormap(flipud(gray));
    xlabel('\phi','FontSize',16);
    ylabel('\theta','FontSize',16);

end


%% Internal Functions
function [r,c,V] = findnearest(srchvalue,srcharray,bias)

% Usage:
% Find the nearest numerical value in an array to a search value
% All occurances are returned as array subscripts
%
% Output:
%
% For 2D matrix subscripts (r,c) use:
%
%       [r,c] = findnearest(srchvalue,srcharray,gt_or_lt)
%
%
% To also output the found value (V) use:
%
%       [r,c,V] = findnearest(srchvalue,srcharray,gt_or_lt)
%
%
% For single subscript (i) use:
%
%         i   = findnearest(srchvalue,srcharray,gt_or_lt)
% 
%
% Inputs:
%
%    srchvalue = a numerical search value
%    srcharray = the array to be searched
%    bias      = 0 (default) for no bias
%                -1 to bias the output to lower values
%                 1 to bias the search to higher values
%                (in the latter cases if no values are found
%                 an empty array is ouput)
%
%
% By Tom Benson (2002)
% University College London
% t.benson@ucl.ac.uk

    if nargin<2
        error('Need two inputs: Search value and search array')
    elseif nargin<3
        bias = 0;
    end

    % find the differences
    srcharray = srcharray-srchvalue;

    if bias == -1   % only choose values <= to the search value

        srcharray(srcharray>0) =inf;

    elseif bias == 1  % only choose values >= to the search value

        srcharray(srcharray<0) =inf;

    end

    % give the correct output
    if nargout==1 | nargout==0

        if all(isinf(srcharray(:)))
            r = [];
        else
            r = find(abs(srcharray)==min(abs(srcharray(:))));
        end 

    elseif nargout>1
        if all(isinf(srcharray(:)))
            r = [];c=[];
        else
            [r,c] = find(abs(srcharray)==min(abs(srcharray(:))));
        end

        if nargout==3
            V = srcharray(r,c)+srchvalue;
        end
    end
end

function [a,th,ph]=equiangle_sampling(N);

% [a,th,ph]=equiangle_sampling(N);
%
% a - weights
% th - theta angles for all points
% ph - phi angles for all points
% Total no. of points is length(a)=length(th)=length(ph)
%
% Equiangle sampling for order N
% From Driscoll and Healy 1994
%
% Boaz Rafaely 21 July 2003

    L=2*(N+1);       % no. samples over the angles

    %theta_eq=0:pi/Leq:(Leq-1)*pi/Leq;               % start at 0
    theta = (0:pi/L:(L-1)*pi/L)+(pi/(2*L)); % symmetric
    phi   = 0:2*pi/L:(L-1)*2*pi/L;

    count=0;
    for j=1:L,
        for k=1:L,
            count=count+1;
            th(count)=theta(j);
            ph(count)=phi(k);
        end;
    end;

    S=0;
    for k=0:N, 
        S=S+(1/(2*k+1))*sin((2*k+1)*th); 
    end;

    a = ((8*pi)/(L^2)) * S .*sin(th);
end