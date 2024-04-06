clear
close
clc

load ../Raw/SBS_VSA_1202RS_PAL.mat;
mobj=SBS_VSA_1202RS_PAL; clear SBS_VSA_1202RS_PAL;

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

eobj.dataDesc='VSA array IR data, single source';
eobj.dataDomain{1} = 'TIME';
eobj.dataDomain{2} = 'SPACE';

eobj.micGrid.azimuth = mobj.azimuth;
eobj.micGrid.elevation = mobj.elevation;
eobj.micGrid.quadType = mobj.quadGrid;
eobj.micGrid.r = mobj.radius;

eobj.sourceGrid.r = mobj.sourceDistance;

eobj.data(1,:,:) = mobj.irChOne;

save('../Raw/earoSBS.mat','eobj');
