function viewphotometry_fox_odor(root, animals)
%VIEWPHOTOMETRY_FOX_ODOR Colour map and average PSTH of simplified fiber photometry data.
%   VIEWPHOTOMETRY(ROOT, ANIMALS) creates colour map and average PSTH of fiber
%   photometry data recorded from single sessions of from path ROOT with 
%   folders named ANIMALS. PETHs are aligned to stimON trigger saved in the
%   processed photometry data. 

%   See also PERIEVENT.M

%   Bálint Király, Panna Hegedus, Balazs Hangya
%   Institute of Experimental Medicine
%   hangya.balazs@koki.mta.hu

% Load data and FiberEvents
datadir = root; %location of fiber photometry data
resdir = fullfile(root, filesep, 'Plots', filesep, 'NoHighpassFilter');

if ~isfolder(resdir)
    mkdir(resdir)
end

if strcmp(animals, 'SOM')
    animalID = {'SOM_HDB_3' 'SOM_HDB_4' 'SOM_HDB_5' 'SOM_HDB_6'};
    channels = {'dff_D' 'dff_A' 'dff_D' 'dff_A'};
    sessionID = '230812';
end

if strcmp(animals, 'PV')
    animalID = {'FS1' 'FS2' 'FS4' 'FS5'};
    channels = {'dff_D' 'dff_D' 'dff_D' 'dff_D'};
    sessionID = '230812';
end

raw_data = cell(1,length(animals));

for i = 1:length(animalID)
    animal = animalID{i};
    path = [datadir filesep animal filesep sessionID filesep];
    raw_data{i} = load([path filesep 'proF.mat']);
    
end

sr = raw_data{1,1}.sr ;
window = [-120 120]; % in seconds
data_points_needed = window(2)*sr;
dt = 1/sr;

% Creating time vector
time_series = window(1):(1/sr):window(2);
numTrials = length(animalID);

matrix = nan(numTrials,length(time_series));
s_matrix = nan(numTrials,length(time_series));

for m=1:numTrials
    stimon_inx = raw_data{1,m}.stimON;
    if strcmp(channels{m}, 'dff_D')
        matrix(m,:)= raw_data{1,m}.dff_D((stimon_inx(end)-data_points_needed):(stimon_inx(end)+data_points_needed));
    elseif strcmp(channels{m}, 'dff_A')
        matrix(m,:)= raw_data{1,m}.dff_A((stimon_inx(end)-data_points_needed):(stimon_inx(end)+data_points_needed));
    end
    % Smoothing PSTH
    sigma = 1;
    [spsth, conv_psth] = smoothed_psth(matrix(m,:),dt,sigma);
    s_matrix(m,:)=spsth;
end

% mean PETH
mean_psth = nanmean(s_matrix);
std_psth = nanstd(s_matrix,1)/sqrt(size(s_matrix,1));

% plot PETH
figure
errorshade(time_series(1:100:end), mean_psth(1:100:end), std_psth(1:100:end),...
    'LineColor',[1 0.5 1],'ShadeColor',[1 0.5 1])