clear
close
clc

[data,fs] = wavread('../Raw/Casa_em32.wav');
load ../Raw/EM_parameters.mat
newfs=48e3;

data=resample(data,newfs,fs);

eobj=earo;

eobj.name='Casa De La Musica';
eobj.context='SRIR';
eobj.location='Italy';
eobj.engineer='Angelo Farina';
eobj.isRobot=false;
eobj.fs=newfs;
eobj.microphone='Eigenmike';
eobj.scatterer=true;

eobj.dataDesc='EM32 array IR data, single source';
eobj.dataDomain{1} = 'TIME';
eobj.dataDomain{2} = 'SPACE';

eobj.micGrid.azimuth = ph_a;
eobj.micGrid.elevation = th_a;
eobj.micGrid.quadType = 'Eigenmike';
eobj.micGrid.r = 0.042;

eobj.data(1,:,:) = data;

save('../Raw/earoCasa.mat','eobj');
