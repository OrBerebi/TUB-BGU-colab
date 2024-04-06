% analyze and create EARO from FABIAN_HRTFs_DAGA15

clear
close
clc

load FABIAN_HRTFs_DAGA15;

fs_orig=44.1e3;
fs=48e3;
hobj_meas = earo;
hobj_bem = earo;

% MEASURED HRTF
hobj_meas.name='FABIAN HRTFs Measured. DAGA15';
hobj_meas.context='HRTF Set';
hobj_meas.engineer='Fabian Brinkman';
hobj_meas.comments='Transformed to EARO by Jonathan Sheaffer';
hobj_meas.isRobot=false;
hobj_meas.type='HRIR';
hobj_meas.fs=fs_orig;
hobj_meas.taps=size(lMeasured,1);
hobj_meas.nData=size(lMeasured,2);

hobj_meas.scatterer=true;

hobj_meas.dataDesc='HRIR Data sources x nData x receivers (1=left ear)';
hobj_meas.dataDomain{1} = 'TIME';
hobj_meas.dataDomain{2} = 'SPACE';

hobj_meas.sourceGrid.azimuth = degtorad(gridMeasured(:,1));
hobj_meas.sourceGrid.elevation = pi/2-degtorad(gridMeasured(:,2));
hobj_meas.sourceGrid.quadType = 'Unknown';
hobj_meas.sourceGrid.quadWeight = 0;
hobj_meas.sourceGrid.r = 1.7;


hobj_meas.data(:,:,1) = lMeasured.';
hobj_meas.data(:,:,2) = rMeasured.';

% BEM HRTF
hobj_bem.name='FABIAN HRTFs BEM. DAGA15';
hobj_bem.context='HRTF Set';
hobj_bem.engineer='Fabian Brinkman';
hobj_bem.comments='Transformed to EARO by Jonathan Sheaffer';
hobj_bem.isRobot=false;
hobj_bem.type='HRIR';
hobj_bem.fs=fs_orig;
hobj_bem.taps=size(lModeled,1);
hobj_bem.nData=size(lModeled,2);

hobj_bem.scatterer=true;

hobj_bem.dataDesc='HRIR Data sources x nData x receivers (1=left ear)';
hobj_bem.dataDomain{1} = 'TIME';
hobj_bem.dataDomain{2} = 'SPACE';

hobj_bem.sourceGrid.azimuth = degtorad(gridModeled(:,1));
hobj_bem.sourceGrid.elevation = pi/2-degtorad(gridModeled(:,2));
hobj_bem.sourceGrid.quadType = 'Unknown';
hobj_bem.sourceGrid.quadWeight = 0;
hobj_bem.sourceGrid.r = 1.7;


hobj_bem.data(:,:,1) = lModeled.';
hobj_bem.data(:,:,2) = rModeled.';

% Resample all data
disp('Resample Measured HRTFs...');
hobj_meas = hobj_meas.resampleData(fs);
disp('Resample BEM HRTFs...');
hobj_bem = hobj_bem.resampleData(fs);

% Save files
disp('Saving files...');
hobj=hobj_meas;
save('FABIAN_HRTFs_DAGA_MEAS','hobj');
hobj=hobj_bem;
save('FABIAN_HRTFs_DAGA_BEM','hobj');

