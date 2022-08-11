clear all;
%% file load
% Loading time should be quite long. So, if you want to use existing data
% in workspace, set load_newfile = 'no'

load_newfile = 'no'; %% 'yes', 'no'

if (~exist('PA_RcvData', 'var')) || strcmp(load_newfile, 'yes')
    clear;
    close all;

    [LoadFilename,LoadFilepath] = uigetfile('..\220729'); %% Choose RcvData
    dir_list = dir(LoadFilepath);

    LoadFilename_part = LoadFilename(1:14);
    
    %load every data of the experiment (US_ImgData, PA_ImgData, PA_RcvData, Parameters)
    for i = 1:size(dir_list,1)
        if (contains(dir_list(i).name, LoadFilename_part))
            load(fullfile(LoadFilepath, dir_list(i).name));
        end
    end
end

%% setting menu
num_frame = 114; %% set frame number of the volume data

%% fixed Parameters
Parameters.NumDepthSample = double(Parameters.NumDepthSample);
Parameters.ImgDepthSample = size(PA_ImgData, 1) * 2;
Parameters.Fs = single(Receive(1).decimSampleRate*1e6); %Hz
Parameters.SoS = single(Resource.Parameters.speedOfSound); % 1450 for US phantom, 1500 for custom phantom
Parameters.factor_interp = 4;
Parameters.Ts = 1/Parameters.Fs;
Parameters.Pitch = single(Trans.spacingMm* 1e-3); % m
Parameters.AxialStep = Parameters.Ts*Parameters.SoS;
Parameters.AxialStep_interp = Parameters.AxialStep / Parameters.factor_interp;
EndDepthSample = Parameters.ImgDepthSample*Parameters.factor_interp;

%% RawData reorder

RawData_temp = PA_RcvData(:,:,num_frame);
RawData = zeros(Parameters.NumDepthSample, size(PA_RcvData,2)*size(RawData_temp, 1)/Parameters.NumDepthSample);
for i = double(1:size(RawData_temp, 1)/Parameters.NumDepthSample)
    RawData(:,(1:size(PA_RcvData,2)) + (i - 1)*size(PA_RcvData,2)) = ...
        RawData_temp((1:Parameters.NumDepthSample) + (i - 1)*Parameters.NumDepthSample, :);
end
% RawData = RawData_temp(1:Parameters.NumDepthSample, :);


%% Parameters
time_series_data = RawData;
sound_speed = Parameters.SoS;        % 1500[m/s]
N_elements = 384;


%acquisition_dict
acquisition_dict.ad_sampling_rate = Parameters.Fs;
acquisition_dict.sizes = size(time_series_data);
acquisition_dict.speed_of_sound = sound_speed;

%Coordinate
device_dict.general.num_detectors = N_elements;
device_dict.general.num_illuminators = 0;
fov = [ 0; (N_elements)*Parameters.Pitch; 0; 0; Parameters.AxialStep_interp; EndDepthSample*Parameters.AxialStep_interp];
device_dict.general.field_of_view = fov;

center_freq = 9.5*1e6;
bw = 95;  % bandwidth

frequency_response=[center_freq, bw];
for det = 1:N_elements
    index = strcat("deleteme", sprintf( '%010d', (det-1) ));
    device_dict.detectors.(index).detector_position = [Parameters.Pitch*(det-0.5) Parameters.AxialStep_interp 0]; 
    device_dict.detectors.(index).detector_geometry_type = "CUBOID";
    device_dict.detectors.(index).detector_geometry = [0 0 0];
    device_dict.detectors.(index).frequency_response = frequency_response;
end


%% Save to hdf5 form

disp("Exporting to the IPASC data format...")
data = pa_data(time_series_data, acquisition_dict, device_dict);
pacfish.write_data(convertCharsToStrings(LoadFilename)+"_"+int2str(num_frame)+".hdf5", data, 1);
