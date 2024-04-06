function [df] = earo_df (eobj, N, kr, theta0, phi0, a_max)
%calculates directivity factor

    az = eobj.micGrid.azimuth;
    el = eobj.micGrid.elevation;
    r = ones(1,length(az))*eobj.micGrid.r;
    k = kr/eobj.micGrid.r;

    % Generate the B-Matrix
    B=zeros(length(az),(N+1)^2);  % B is a Q x (N+1)^2 x length(k)
    Y = shMatrix(N,el,az).';
    
    if eobj.scatterer, sphereType = 1; else sphereType = 0; end;
    
    bn = radialMatrix(N, kr, sphereType, a_max, true, false);
    B = Y*diag(bn);
    Bdag = pinv(B);

    fprintf('.');    
    
    
    % Generate plane wave 
    p0 = exp(-1i.*kr.*angleTwoPoints(el,az,theta0,phi0));

    
    % Generate integration quadrature
    [a,theta_l,phi_l]=equiangle_sampling(N);
    Yl=shMatrix(N,theta_l,phi_l);
    
    for ll=1:length(phi_l)      % iterate through look directions
        tmp=Bdag.'*Yl(:,ll);
        wconj=(4*pi/((N+1)^2))*tmp;  
        wconjArr(ll,:)=wconj;

        % Beamformer output
        y(ll)=wconj'*p0.';
    end
    
    % Find closest index in quadrature
    [Xg,Yg,Zg]   = sph2cart(phi_l,theta_l-pi/2,1);
    [Xp,Yp,Zp]   = sph2cart(phi0,theta0-pi/2,1);
    distances    = sqrt((Xg-Xp).^2+(Yg-Yp).^2+(Zg-Zp).^2);
    id           = max(find(distances == min(distances)));

    % Compute DF
    df = abs(y(id))^2 / ((1/(4*pi)) * sum(a.*abs(y).^2));    


end

function out = angleTwoPoints(theta1,phi1,theta2,phi2)
% function out = angleTwoPoints(theta1,phi1,theta2,phi2)
%
% calculate angle between two points on a sphere
%

  out = sin(theta1).*cos(phi1).*sin(theta2).*cos(phi2) + sin(theta1).*sin(phi1).*sin(theta2).*sin(phi2) + cos(theta1).*cos(theta2);
end

function [a,th,ph]=equiangle_sampling(N)

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
