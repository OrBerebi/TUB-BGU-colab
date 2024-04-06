function mic_s = EM32toEARO(mic_filename,save_flag,exp_location,start_nData_sec,max_nData_sec)
%% Funcion to transform wav file of measurements from EigenMic to earo object
%
% INPUTS:
%       mic_filename    : string, file with EM recordings
%       save_flag       : if true, save earo object to file
%       exp_location    : string for earo object 'location' field
%       start_nData_sec : start point (in seconds) to trim the recording
%       max_nData_sec   : maximum length of recording to use (in seconds)
% 
% OUTPUTS:
%       mic_s           : earo object contains EM recordings
%
% Updated by Zamir Ben-Hur
% 25.7.18

if nargin<2
    mic_filename='Casa_em32';
    %     source_label='S';
    save_flag=0;
end
if nargin<3
    exp_location='exp_location';
end
if nargin<5
    max_nData_sec=10; %maximal length of data [sec]
    start_nData_sec=0; %staring point of data [sec]
end

mic_s=earo;

load('EM32_parameters.mat')

[filepath,exp_file_name,~] = fileparts(mic_filename);

exp_file= fullfile(filepath,exp_file_name);

%%
mic_s.name               = [exp_file_name];
mic_s.context            = [];
mic_s.location           = exp_location;
mic_s.date               = 2018;
mic_s.engineer           = 'researcher';
mic_s.contact            = '###@post.bgu.ac.il';
mic_s.comments           = [];
mic_s.naoPhase           = [];
mic_s.isRobot            = false;
mic_s.earoVersion        = 1.0;
mic_s.type               = 'MICARRAY';
mic_s.fs                 = 48000;
mic_s.taps               = [];
mic_s.nData              = [];
mic_s.excitationSignal   = [];
mic_s.gapTime            = [];
mic_s.microphone         = [];
mic_s.source             = [];
mic_s.audioInterface     = [];
mic_s.micPreamp          = [];
mic_s.capturingSystem    = 'EM32';
mic_s.systemLoopLatency  = [];
mic_s.latencyCompensated = [];
mic_s.headCut            = [];
mic_s.avgAirTemp         = [];
mic_s.avgRelHumidity     = [];
mic_s.positionReference  = [];
mic_s.postProcessing     = [];
mic_s.scatterer          = true;   %rigid sphere
mic_s.sourceGrid         = struct('azimuth',[],'elevation',[],'r',[],'quadType',[]);
mic_s.micGrid            = struct('azimuth',ph_mic,'elevation',th_mic,'r',r_mic,'quadType',[]);
mic_s.orderN             = 4;
mic_s.dataDesc           = [];
mic_s.dataDomain         = [{'TIME'},{'SPACE'}];
mic_s.data               = [];
mic_s.shutUp             = false;
mic_s.angles             = 'RAD';


%%

disp('Generating EARO data structure');

fprintf('loading file number:   ');
[temp,fs] = audioread([exp_file,'.wav']);
if fs~=48e3
     [p,q]=rat(48000/fs);
     temp=resample(temp,p,q);
end
if size(temp,1)>(max_nData_sec+start_nData_sec)*mic_s.fs
    temp=temp(round(start_nData_sec*mic_s.fs)+[1:round(max_nData_sec*mic_s.fs)],:);
end
mic_s.data(1,:,:)=temp;

if save_flag
    eobj = mic_s;
    save(['./eobj_',exp_file_name,'.mat'],'eobj');
end
1;