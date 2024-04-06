%FIG_SH_BALLOON generates Figures 1.5-1.8, which show blloon plots of 
% spherical harmonics functions for various orders and degrees and from 
% various view points. 
%
% Fundmentals of Spherical Array Processing
% Boaz Rafaely, 2017.

close all;
clear all;

path(path,'../../math');
path(path,'../../plot');

S=100;
N=4;
[x,y,z] = sphere(S);
[th,ph,r]=c2s(x,y,z);
th1=reshape(th,(S+1)^2,1);
ph1=reshape(ph,(S+1)^2,1);
SH=sh2(N,th1,ph1);

for fig=1:4,
    
figure(fig);
for n=0:4,
    for m=-n:n,
        subplot(5,9,9*n+5+m),
        if m<0, c1=imag(SH); else c1=real(SH); end;
        c=reshape(c1(n^2+n+m+1,:),S+1,S+1);
        [X,Y,Z]=s2c(th,ph,abs(c));
        
        h3=surf(X,Y,Z,sign(c)); 
        axis([-1 1 -1 1 -1 1]);
        axis off

        colormap(flipud(cool));
        brighten(0.9);
        shading faceted;
        set(h3,'EdgeAlpha',0.01);
        caxis([-1 1]);
        
        lighting phong
        if fig==1, view(3);     camlight('right'); camzoom(3.5); end; %standard Matlab view
        if fig==2, view(0,90);  camlight('right'); camzoom(3);   end; % view from the +z axis (y vs. x)
        if fig==3, view(90,0);  camlight('right'); camzoom(2.3); end; % view from the +x axis (z vs. y)
        if fig==4, view(180,0); camlight('right'); camzoom(2.3); end; % view from the +y axis (z vs. x)

   end;
end;

end;

% figure(1); print -dpng -r600 ../../figures/chapter01/fig_sh_balloon_angle.png % angular view
% figure(2); print -dpng -r600 ../../figures/chapter01/fig_sh_balloon_z.png     % +z axis view
% figure(3); print -dpng -r600 ../../figures/chapter01/fig_sh_balloon_x.png     % +x axis view
% figure(4); print -dpng -r600 ../../figures/chapter01/fig_sh_balloon_y.png     % +y axis view

