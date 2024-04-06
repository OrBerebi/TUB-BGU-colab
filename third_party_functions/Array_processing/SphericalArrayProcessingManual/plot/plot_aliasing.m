function plot_aliasing(E)
%PLOT_ALIASING plots an aliasing matrix E.
% plot_aliasing(E).
%
% Fundmentals of Spherical Array Processing
% Boaz Rafaely, 2017.

AxisFontSize=14;

N=floor(sqrt(size(E,1)-1));
NN=floor(sqrt(size(E,2)-1));
    
E=abs(E);

E(E<1e-10)=0;
E=E/max(E(:)); % normalize
normFactor=min(E(E~=0));
if(isempty(normFactor)), normFactor=1; end;
E=E+normFactor/100; % Avoid zeros...
normFactorDb=floor(20*log10(normFactor)/10)*10;
clims = [ normFactorDb-10 0 ];
imagesc(0:size(E,2)-1,0:size(E,1)-1,20*log10(abs(E)),clims );
axis('image');
set(gca,'YDir','normal','FontSize',AxisFontSize)
colormap(gray);
h=colorbar;
set(gca,'FontSize',AxisFontSize);

% colorbar title
colorbarTitleHandler=get(h,'title');
set(colorbarTitleHandler,'String','(dB)','FontSize',AxisFontSize);

xlabel('$n''^2+n''+m''$','Interp','Latex','FontSize',AxisFontSize);
ylabel('$n^2+n+m$','Interp','Latex','FontSize',AxisFontSize);

hold on;
picLimX=[0 size(E,2)-1];
picLimY=[0 size(E,1)-1];

ns=0:N;
qLimY=ns.^2-0.5;
for curQ=qLimY
    plot(picLimX,[curQ curQ],'c-');
end

ns=0:NN;
qLimX=ns.^2-0.5;    
for curQ=qLimX
    plot([curQ curQ],picLimY,'c-');
end
end