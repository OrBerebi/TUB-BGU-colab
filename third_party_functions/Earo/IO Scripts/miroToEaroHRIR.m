clear
close
clc

load ../Raw/HRIR_L2702.mat;
mobj=HRIR_L2702; clear HRIR_L2702;

eobj=earo;

eobj.name=mobj.name;
eobj.context=mobj.context;
eobj.location=mobj.location;
eobj.engineer=mobj.engineer;
eobj.contact=mobj.contact;
eobj.comments=mobj.comments;
eobj.isRobot=false;
eobj.type=mobj.type;
eobj.fs=mobj.fs;
eobj.taps=mobj.taps;
eobj.nData=mobj.nIr;
eobj.excitationSignal=mobj.excitationSignal;
eobj.gapTime=mobj.gapTime;
eobj.microphone=mobj.microphone;
eobj.source=mobj.source;
eobj.audioInterface=mobj.audioInterface;
eobj.micPreamp=mobj.micPreamp;
eobj.capturingSystem=mobj.capturingSystem;
eobj.systemLoopLatency=mobj.systemLoopLatency;
eobj.latencyCompensated=mobj.latencyCompensated;
eobj.headCut=mobj.headCut;
eobj.avgAirTemp=mobj.avgAirTemp;
eobj.avgRelHumidity=mobj.avgRelHumidity;
eobj.positionReference=mobj.positionReference;
eobj.postProcessing=mobj.postProcessing;
eobj.scatterer=mobj.scatterer;

eobj.dataDesc='HRIR Data sources x nData x receivers (1=left ear)';
eobj.dataDomain{1} = 'TIME';
eobj.dataDomain{2} = 'SPACE';

eobj.sourceGrid.azimuth = mobj.azimuth;
eobj.sourceGrid.elevation = mobj.elevation;
eobj.sourceGrid.quadType = mobj.quadGrid;
eobj.sourceGrid.quadWeight = mobj.quadWeight;
eobj.sourceGrid.r = mobj.sourceDistance;

eobj.micGrid.r = mobj.radius;

eobj.data(:,:,1) = mobj.irChOne.';
eobj.data(:,:,2) = mobj.irChTwo.';

hobj=eobj;

save('../Raw/earoHRIR.mat','hobj');
