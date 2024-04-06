function hobj = sofaToEaro(SOFAfile,outFile)
%
% SOFAfile = 'FABIAN_multiHATO_2degGCD.sofa';
% outFile = 'earoFABIAN.mat';
noSave = 0;
if nargin<2
    noSave = 1;
end

mobj=SOFAload(SOFAfile,'nochecks');
eobj=earo;

%% Cart2SPH
if strcmp(mobj.SourcePosition_Type,'cartesian')
    [r_grid_sofa,el_grid_sofa,az_grid_sofa]=c2s(mobj.SourcePosition(:,1),mobj.SourcePosition(:,2),mobj.SourcePosition(:,3));
elseif strcmp(mobj.SourcePosition_Type,'spherical')
    if strcmp(mobj.SourcePosition_Units,'degree, degree, metre')
        r_grid_sofa = mobj.SourcePosition(:,3);
        el_grid_sofa = mobj.SourcePosition(:,2)*pi/180;
        el_grid_sofa = pi/2 - el_grid_sofa;
        az_grid_sofa = mobj.SourcePosition(:,1)*pi/180;
        az_grid_sofa(az_grid_sofa<0) = az_grid_sofa(az_grid_sofa<0) + 2*pi;
    elseif strcmp(mobj.SourcePosition_Units,'radian, radian, metre')
        r_grid_sofa = mobj.SourcePosition(:,3);
        el_grid_sofa = mobj.SourcePosition(:,2);
        el_grid_sofa = pi/2 - el_grid_sofa;
        az_grid_sofa = mobj.SourcePosition(:,1);
    else
        error('Unknown sourcePosition units')
    end
else
    error('Unknown sourcePosition type')
end

eobj.name=mobj.GLOBAL_ListenerShortName;
eobj.context=mobj.GLOBAL_DatabaseName;
eobj.date   = mobj.GLOBAL_DateCreated;
eobj.engineer=mobj.GLOBAL_AuthorContact;
if isfield(mobj,'GLOBAL_ReceiverDescription')
    hobj.microphone     = mobj.GLOBAL_ReceiverDescription;
else
    hobj.microphone     = [];
end
if isfield(mobj,'GLOBAL_EmitterDescription')
    eobj.source     = mobj.GLOBAL_EmitterDescription;
else
    eobj.source     = [];
end
if isfield(mobj,'GLOBAL_Comment')
    eobj.comments=mobj.GLOBAL_Comment;
else
    eobj.comments= [];
end
eobj.isRobot=false;
eobj.type=mobj.GLOBAL_DataType;
eobj.fs=mobj.Data.SamplingRate;
eobj.taps=size(mobj.Data.IR,3);
eobj.nData=size(mobj.Data.IR,1);

eobj.scatterer=true;

eobj.dataDesc='HRIR Data sources x nData x receivers (1=left ear)';
eobj.dataDomain{1} = 'TIME';
eobj.dataDomain{2} = 'SPACE';

eobj.sourceGrid         = struct('azimuth',az_grid_sofa,'elevation',el_grid_sofa,'r',r_grid_sofa,'quadType',mobj.SourcePosition_Type ,'quadWeight',0);
eobj.micGrid            = struct('azimuth',[],'elevation',[],'r',[],'quadType',[],'quadWeight',[]);

eobj.angles             = 'RAD';

eobj.data = permute(mobj.Data.IR,[1 3 2]);

hobj=eobj;

if ~noSave
    save(outFile,'hobj');
end