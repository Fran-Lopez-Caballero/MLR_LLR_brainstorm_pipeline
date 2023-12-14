%% Brainstorm EEG + MEG analysis pipeline

% INFORMATION: Script to perform all analyses of Project project data (human side)
    % CNRL, Fran Lopez Caballero 05/03/2020
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set Brainstorm nogui on completely headless mode for LLR 
% Set a window in server first, then run this section and run scripts normally 
% (it won't break despite bad remote connection this way -16/03/2020-)

% Set up the Brainstorm files
clear
addpath('~/matlab/brainstorm3'); 
BrainstormDbDir = '~/brainstorm_db';
% Start Brainstorm
if ~brainstorm('status')
    brainstorm server
end
bst_set('BrainstormDbDir',BrainstormDbDir)
% Select the correct protocol
ProtocolName = 'Project_LLRraw'; % Enter the name of your protocol
sProtocol.Comment = ProtocolName;
% sProtocol.SUBJECTS = [home 'anat'];
sProtocol.SUBJECTS = '~/brainstorm_db/Project_LLRraw/anat';
sProtocol.STUDIES = '~/brainstorm_db/Project_LLRraw/data';
db_edit_protocol('load',sProtocol);
% Get the protocol index
iProtocol = bst_get('Protocol', ProtocolName);
if isempty(iProtocol)
    error(['Unknown protocol: ' ProtocolName]);
end
% Select the current procotol
gui_brainstorm('SetCurrentProtocol', iProtocol);

%% Define variables for LLR

% Depending on whether first step was already launched in a screen
% clear
% load '/private/path/project/Brainstorm_pipelines/Brainstorm_Vars_LLRraw'
% parloop()
root_dir = '/private/path/project';
root_dir_bs = '~/brainstorm_db/Project_LLRraw'; 

% get protocol name, in case we don't run it with server mode
pos_last = find(root_dir_bs == '/', 1, 'last');
ProtocolName = root_dir_bs(pos_last+1:end); clear pos_last

participant = {
    'xxxx'
};

participant_general_list = participant; % so that it has a different name
% Load info of sessions and blocks per subject
load([root_dir '/Brainstorm_pipelines/session_block_array_original.mat'])
% MAY WANNA USE session_block_array_original_LLR since 2320F has no triggers for LLR only
% load([root_dir '/Brainstorm_pipelines/session_block_array_sources.mat'])
load([root_dir '/Brainstorm_pipelines/additional_bad_chans.mat'])
% stimulus_artefact = {'2193', '2262', '2274'; {'A', 'B', 'D'}, {'A', 'B', 'C', 'E', 'F'}, {'A', 'C', 'D', 'G'}};
% stimulus artefact correction not needed for LLR
session = {'A' 'B' 'C' 'D' 'E' 'F' 'G' 'H' 'I' 'J'}; % default one, will be redefined
condition = {'11' '12' '13' '31' '32' '33' '51' '52' '53' '71' '72' '73' '91' '92' '93'};
% If condition changes (for instance, running shortest ISI alone, careful with c=4:length...
shortest_ISI = {'11' '12' '13'};
loudest_condition = {'13' '33' '53' '73' '93'}; % only for stimulus artefact correction
Exp_cond = {'Quietest', 'Medium_dB', 'Loudest', 'Fastest', 'Fast',...
                'Medium_ISI', 'Slow', 'Slowest'};
block = {'b1' 'b2' 'b3'}; % The default one, will actually be adapted with the block_and_runs variable
runs = {'1' '2' '3' '4' '5' '6' '7' '8' '9'}; % same than before
modality_data = {'EEG','MEG','both_mod'};
wave = {'LLR'}; % inherited from previous structure, there is loop with only one value
epoch_wave = {[-0.15, 0.4]}; % LLR only
epoch_wave_bas = {[-0.15, 0.6]}; % LLR, only to create templates for correction of shortest ISI
epoch_baseline = {[-0.15, 0]}; % LLR baseline correction
epoch_covar = {[-3.5, 0.015]}; % from which we calculate noise cov matrices
time_noise_cov = [-3.4, -0.01]; % time window from which to obtain noise covariance values
reject_EEG = [0, 100]; % peak to peak threshold (change name maybe)
reject_EEG_absolute = [0, 50]; % absolute threshold
reject_MEG_GRAD_absolute = [0, 2500];
reject_MEG_MAG_absolute = [0, 2500];
reject_MEG_grad = [0, 5000];
reject_MEG_mag = [0, 5000];
MLR_highpass = 0; % should be 10, but  already high passsed at preprocessing
MLR_lowpass = 200;
LLR_highpass = 0; % None (because it's already filtered at  0.5 from before)
LLR_lowpass = 20; % Tobias wanted to try something like 80 Hz. Maybe 55 for  time-frequency analysis?
tranband_MLR = 3.7; % Transition bandwith, useful to avoid filter artifacts and reduce filter order
tranband_LLR = 0; % Leave at 0 for now
resample_LLR = 500; % LLR resample from 1500 to 500Hz
wave_name = {'low'}; % since both MLR and LLR will only be low_pass in Brainstorm
reref_option = 1; % 0 = NO 1 = YES Yes or no rereference
ref_EEG = 'M1, M2'; % in case we use it, but we won't for now: Alternative: 'AVERAGE'
detrend_option = 0; % 0 = NO detrend; 1 = Brainstorm all epoch; 2 = BVA from first to last 50ms
detrend_option_templates = 0; % 0 = NO detrend; 1 = Brainstorm all epoch; 2 = BVA from first to last 50ms
subtraction_approach =  1; % 0 = No subtraction produced;   1 = YES.
covar_epochs = 1; % 0 = NO covar epochs; 1 = YES, make them
sensor_analysis = 4; % 1 = EEG only; 2 = MEG only; 3 = Combined EEG and MEG only; 4 = ALL
source_analysis = 4; % 1 = EEG only; 2 = MEG only; 3 = Combined EEG and MEG only. 4 = ALL
sources_single_session = 0; % 0 = NO; 1 = YES % Average sources so that we have each 11, 12, etc per session (unnecessary)
source_noise_options = {'reg'}; % For 'NoiseMethod'; 'median'; 'diag'
source_noise_tags = {'Source_Regul'}; % 'Source_Eigen' 'Source_diag'
source_inverse_measure = 'dspm2018'; % 'amplitude' or 'dspm2018'
reset_channel_files = 0; % 0 = leave current channel files, 1 = reset them to backup channel files
delete_previous_file = 1; % 1 = delete previous head models and source kernels if reran these steps
copy_all_transf = 1; % 0 = copy only last manual corregistration. 1 = copy all
string_to_look = 'Apr-2021';
ERROR_AVERAGE = {}; % variable to store errors when no average is possible
group_default_cortex = 'tess_cortex_pial_02.mat'; % cortical surface where all cortical sources of invididual subjects will be projected
backup_corregistration_folder = [root_dir '/Channel_backup/NEW_APRIL_2021_EEG_projected'];

% Variables for subtraction approach
hann_wind_subtract = 0; % 0 = NO; 1 = YES Apply hanning window to vector to subtract
hann_wind_template = 1; % 0 = NO; 1 = YES Apply half hanning to end of template
hann_start = 0.4; % For now apply from 400 to the end (600 ms). NOT CONSIDERING BASELINE!
all_val = 2; % 1 = all files 2 = all files that are going to be analized
% Note: if you ever change the length of the baseline portion for MLR or
% LLR, change it in the script to correct with subtraction
template = 2; % 1 = block average; 2 = session average; 3 = total average 
% PROBABLY BEST OPTION FOR TEMPLATE IS 2 (Tradeoff between cleanliness of
% templates and affecting MEG data)
delimiter = ',';
startRow = 12;
formatSpec = '%*q%q%q%*s%*s%*s%[^\n\r]';
Non_averaged_across_blocks = {}; naab = 1;

initialVars = who; % variables up until here, which won't be deleted afterwards
initialVars = who; % twice so that InitialVars itself is not deleted

%% Run all sections for 1 subject at a time (except if they belong to diff protocols)

% If willing to do so, put the sections in this loop
for par = 1:length(participant_general_list) 
    participant = {}; % rebuild the participant variable everytime
    participant{1} = participant_general_list{par};
    %%%%%%%%%%% CODE SECTIONS HERE %%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%% Import raw data LLR (EEG/MEG) (import Anatomy with HCP is to be used in local with CNRL_import_anatomy)
% Its going to import files liisted in session_block_array. Double check
% that they are the same 
tic
disp(' ');      
disp('-------------------------');  
disp('IMPORTING EEG/MEG DATA FOR LLR');
disp(datetime)
disp('-------------------------');     
disp(' '); 

Import_error_log_LLR = {};
error_import_LLR = 1;

for p = 1:length(participant)
    
    % Define sessions within this participant
    pos_par = find(strcmp(session_block_array(:,1), participant{p}(1:4)));
    session = session_block_array{pos_par,2}(1,:);
    
    for s = 1:length(session)
        
        folders = dir([root_dir '/analysis/' participant{p}(1:4) session{s} '/']);
        % Define blocks within this session
        pos_ses = find(strcmp(session_block_array{pos_par,2}(1,:), session{s}));
        block = session_block_array{pos_par,2}{2,pos_ses}(1,:);
        for b = 1:length(block)
            % Define runs within this block
            pos_block = find(strcmp(session_block_array{pos_par,2}{2,pos_ses}(1,:), block{b}));
            runs_in_block = session_block_array{pos_par,2}{2,pos_ses}{2,pos_block}(1,:);
            
            for r = 1:length(runs_in_block)    
                results =  endsWith({folders.name},'tsss_AMICA.fif') & (...
                    contains({folders.name},[participant{p}(1:4) session{s} '_' runs_in_block{r}]));
                infolder = find(results);
                if isempty(infolder)
                    Import_error_log_LLR{error_import_LLR} = ([participant{p}(1:4) session{s} '_' runs_in_block{r} ': tsss AMICA not found']);
                    error_import_LLR = error_import_LLR +1;
                    continue
                end
                file_name = [root_dir '/analysis/' participant{p}(1:4) session{s} '/' folders(infolder).name]; 
                
                disp(' ');      
                disp('-------------------------');  
                disp(['Importing data EEG/MEG LLR data for ' participant{p} '_' session{s} '_' block{b} '_' runs_in_block{r}]);
                disp(datetime)
                disp(' '); 

                sFiles = [];
                % Process: Create link to raw file
                sFiles = bst_process('CallProcess', 'process_import_data_raw', sFiles, [], ...
                    'subjectname',    participant{p}, ...
                    'datafile',       {file_name, 'FIF'}, ...
                    'channelreplace', 0, ...
                    'channelalign',   0, ...
                    'evtmode',        'value');          
            end
        end        
    end
end
save([root_dir '/Events/Import_error_log_LLR.mat'],'Import_error_log_LLR')
clearvars('-except', initialVars{:});

disp 'DONE WITH IMPORTING DATA EEG/MEG for LLR!!!'
disp(datetime)
toc

%% Remove ANY previous active projectors for LLR before analyzing the data
tic
disp(' ');      
disp('-------------------------');  
disp('REMOVING ANY ACTIVE PROJECTORS FOR LLR BEFORE ANYTHING');  
disp(datetime)
disp('-------------------------');     
disp(' '); 

for p = 1:length(participant)
    
    % Define sessions within this participant
    pos_par = find(strcmp(session_block_array(:,1), participant{p}(1:4)));
    session = session_block_array{pos_par,2}(1,:);
    
    for s = 1:length(session)
    % Find the raw folders for each participant and session
        folders = dir([root_dir_bs '/data/' participant{p} '/']);
        results =  contains({folders.name},'raw') & contains({folders.name},[participant{p}(1:4) session{s}]) & endsWith({folders.name},'tsss_AMICA');
        infolder = find(results);
            if isempty(infolder)
                continue
            end
        for l= 1:length(infolder)
            line = infolder(l);
            file_name = [root_dir_bs '/data/' participant{p} '/' folders(line).name '/channel_vectorview306_acc1.mat'];
            eval(['load ' file_name])
            
            disp(' ');      
            disp('-------------------------');  
            disp(['Removing active projectors for LLR before anything: ' participant{p} '_' session{s}]);
            disp(datetime)
            disp(' '); 
            
            if ~isempty(Projector)
                num_proj = size(Projector,2);
                for npj = 1:num_proj
                    Projector(npj).Status = 0;
                end
            end
            save(file_name, 'Channel', 'Comment', 'HeadPoints', 'History', 'IntraElectrodes', 'MegRefCoef', 'Projector', 'SCS', 'TransfEeg', 'TransfEegLabels', 'TransfMeg', 'TransfMegLabels')
        end
    end
end


clearvars('-except', initialVars{:});
disp 'DONE WITH REMOVING ANY ACTIVE PROJECTORS FOR LLR BEFORE ANYTHING!!!'
disp(datetime)
toc

%% Import events to raw data LLR
tic
disp(' ');      
disp('-------------------------');  
disp('IMPORTING EVENTS TO RAW DATA FOR LLR');  
disp(datetime)
disp('-------------------------');     
disp(' '); 

for p = 1:length(participant)
    
    disp(' ');      
    disp('-------------------------');
    disp(['importing events for LLR participant ' participant{p}]);
    disp(datetime)
    disp(' '); 
    
    prot_subs = bst_get('ProtocolSubjects');
    current_sub = find(strcmp({prot_subs.Subject.Name}, participant{p}));
    db_reload_conditions(current_sub);
    
    % Define sessions within this participant
    pos_par = find(strcmp(session_block_array(:,1), participant{p}(1:4)));
    session = session_block_array{pos_par,2}(1,:);
    
    for s = 1:length(session)
    % Find the raw folders for each participant and session
        folders = dir([root_dir_bs '/data/' participant{p} '/']);
        results =  contains({folders.name},'raw') & contains({folders.name},[participant{p}(1:4) session{s}]) & endsWith({folders.name},'tsss_AMICA');
        infolder = find(results);
            if isempty(infolder)
                continue
            end
        file_names = {};
        for l= 1:length(infolder)
            line = infolder(l);
            file_names{l} = [participant{p} '/' folders(line).name '/data_0raw_' folders(line).name(5:end) '.mat'];
        end
        
        disp(' ');      
        disp('-------------------------');  
        disp(['Importing events LLR for ' participant{p} '_' session{s}]);
        disp(datetime)
        disp(' '); 
            
        % Input files (contains all runs found for that participant and session)
        sFiles = file_names;

        if isempty(sFiles)
            error('for whatever reason, sFiles is empty, my friend...')
        end
        
        % Process: CNRL rename EEG Channels
        sFiles = bst_process('CallProcess', 'process_cnrl_rename_eeg', sFiles, [], ...
            'action', 2);  % EasyCap
        
        % Process: Set channels type
        sFiles = bst_process('CallProcess', 'process_channel_settype', sFiles, [], ...
        'sensortypes', 'HEOG,VEOG', ...
        'newtype',     'EOG');

        % Process: Set channels type
        sFiles = bst_process('CallProcess', 'process_channel_settype', sFiles, [], ...
        'sensortypes', 'ECG', ...
        'newtype',     'ECG');
        
        if reref_option == 1     
        % Process: Re-reference EEG
        sFiles = bst_process('CallProcess', 'process_eegref', sFiles, [], ...
            'eegref',      ref_EEG, ...
            'sensortypes', 'EEG');     
        end
        
        % Process: Set channels type: addition on December 2020
        sFiles = bst_process('CallProcess', 'process_channel_settype', sFiles, [], ...
            'sensortypes', 'M1,M2', ...
            'newtype',     'EEG REF');
        
        % Process: Run Matlab command
        sFiles = bst_process('CallProcess', 'process_matlab_eval', sFiles, [], ...
            'matlab',      ['% Available variables: Data, TimeVector' 10 '' 10 'TimeVector = TimeVector - TimeVector(1);' 10 ''], ...
            'sensortypes', '');

        % Process: Import events file
        sFiles = bst_process('CallProcess', 'process_CNRL_evt_import', sFiles, [], ...
            'evtname', '11, 12, 13, 31, 32, 33, 51, 52, 53, 71, 72, 73, 91, 92, 93, boundary');

        % Process: Convert to simple event
        sFiles = bst_process('CallProcess', 'process_evt_simple', sFiles, [], ...
            'eventname', 'boundary', ...
            'method',    1);  % Keep the start of the events

        % Process: Rename event
        sFiles = bst_process('CallProcess', 'process_evt_rename', sFiles, [], ...
            'src',  'boundary', ...
            'dest', 'boundaryStart');

        % Process: Import events file
        sFiles = bst_process('CallProcess', 'process_CNRL_evt_import', sFiles, [], ...
            'evtname', 'boundary');

        % Process: Convert to simple event
        sFiles = bst_process('CallProcess', 'process_evt_simple', sFiles, [], ...
            'eventname', 'boundary', ...
            'method',    3);  % Keep the end of the events

        % Process: Rename event
        sFiles = bst_process('CallProcess', 'process_evt_rename', sFiles, [], ...
            'src',  'boundary', ...
            'dest', 'boundaryEnd');

        for c = 1:length(condition)
        % Process: Remove simultaneous End
        sFiles = bst_process('CallProcess', 'process_evt_remove_simult', sFiles, [], ...
            'remove', condition{c}, ...
            'target', 'boundaryEnd', ...
            'dt',     0.15, ...
            'rename', 0);
        
        % Process: Remove simultaneous Start
        sFiles = bst_process('CallProcess', 'process_evt_remove_simult', sFiles, [], ...
            'remove', condition{c}, ...
            'target', 'boundaryStart', ...
            'dt',     0.4, ...
            'rename', 0);
        end
        
    end
end

clearvars('-except', initialVars{:});
disp 'DONE WITH IMPORTING EVENTS TO RAW DATA FOR LLR!!!'
disp(datetime)
toc

%% Filtering
tic

disp(' ');      
disp('-------------------------');  
disp('FILTERING FOR LLR');  
disp(datetime)
disp('-------------------------');     
disp(' '); 

% Import event data this time (matlab folders)
for p = 1:length(participant)
    
    disp(' ');      
    disp('-------------------------');
    disp(['loading participant ' participant{p}]);
    disp(datetime)
    disp(' '); 
    
    prot_subs = bst_get('ProtocolSubjects');
    current_sub = find(strcmp({prot_subs.Subject.Name}, participant{p}));
    db_reload_conditions(current_sub);
    
    % Define sessions within this participant
    pos_par = find(strcmp(session_block_array(:,1), participant{p}(1:4)));
    session = session_block_array{pos_par,2}(1,:);
    
    for s = 1:length(session)
        % Find the event folders for each participant and session
        folders = dir([root_dir_bs '/data/' participant{p} '/']);
        results = contains({folders.name},'raw') & contains({folders.name},[participant{p}(1:4) session{s}]) & endsWith({folders.name},'matlab');
        infolder = find(results);
        
        if isempty(infolder)
            continue
        end
        file_names = {};
        for l= 1:length(infolder)
            line = infolder(l);
            file_names{l} = [participant{p} '/' folders(line).name '/data_0raw_' folders(line).name(5:end) '.mat'];
        end

        % Input files
        sFiles = file_names;
        
        if isempty(sFiles)
            error('for whatever reason, sFiles is empty, my friend...')
        end
        
        disp(' ');      
        disp('-------------------------');  
        disp(['Filtering for LLR ' participant{p} '_' session{s}]);
        disp(datetime)
        disp(' ');  

        % Process: LLR filter (low pass)
        LLR_Files = bst_process('CallProcess', 'process_bandpass', sFiles, [], ...
            'sensortypes', '', ...
            'highpass',    LLR_highpass, ... 
            'lowpass',     LLR_lowpass, ... 
            'tranband',    tranband_LLR, ...
            'attenuation', 'strict', ...  % 60dB
            'ver',         '2019', ...  % 2019
            'mirror',      0, ...
            'read_all',    0);

        % Process: Add tag: LLR
        LLR_Files = bst_process('CallProcess', 'process_add_tag', LLR_Files, [], ...
            'tag',           'LLR', ...
            'output',        1);  %#ok<*NASGU> % Add to comment
        
        %%% Apparently, after filtering projectors are no longer active, so
        %%% we have to activate the last one applied before epoching the data
        if reref_option == 1 
        folders = dir([root_dir_bs '/data/' participant{p} '/']);
        results =  contains({folders.name},'raw') & contains({folders.name},[participant{p}(1:4) session{s}]) & endsWith({folders.name},'low');
        infolder = find(results);
        for l = 1:length(infolder) 
            % We loop again through each folder ending in 'low', only this time is the whole path
            line = infolder(l);
            file_name = ['~/brainstorm_db/' ProtocolName '/data/' participant{p} '/' folders(line).name '/channel_vectorview306_acc1.mat'];
            eval(['load ' file_name])
            if ~isempty(Projector)
                npj = size(Projector,2); % last projector applied
                Projector(npj).Status = 1;
            end
            save(file_name, 'Channel', 'Comment', 'HeadPoints', 'History', 'IntraElectrodes', 'MegRefCoef', 'Projector', 'SCS', 'TransfEeg', 'TransfEegLabels', 'TransfMeg', 'TransfMegLabels')
            % Not necessary to reload afterwards, as in next steps the
            % whole subject folder is reloaded
        end
        end
    end
end

clearvars('-except', initialVars{:});
disp 'DONE WITH FILTERING DATA FOR MLR!!!'
disp(datetime)
toc

%% Epoching to create TEMPLATES (Part 1)

if subtraction_approach ==  1
tic
disp(' ');      
disp('-------------------------');  
disp('OBTAINING ELONGATED LLR EPOCHS FOR BASELINE CORRECTION (PART 1)'); 
disp(datetime)
disp('-------------------------');     
disp(' ');

% Import filtered data this time (matlab_high and matlab_low folders)
for p = 1:length(participant)
    
    disp(' ');      
    disp('-------------------------');
    disp(['loading participant ' participant{p}]);
    disp(datetime)
    disp(' ');     
    
    prot_subs = bst_get('ProtocolSubjects');
    current_sub = find(strcmp({prot_subs.Subject.Name}, participant{p}));
    db_reload_conditions(current_sub);
    
    % Define sessions within this participant
    pos_par = find(strcmp(session_block_array(:,1), participant{p}(1:4)));
    session = session_block_array{pos_par,2}(1,:);
    
    
    for s = 1:length(session)
        % Define blocks within this session
        pos_ses = find(strcmp(session_block_array{pos_par,2}(1,:), session{s}));
        block = session_block_array{pos_par,2}{2,pos_ses}(1,:);
        
        for b = 1:length(block)
            % Define runs within this block
            pos_block = find(strcmp(session_block_array{pos_par,2}{2,pos_ses}(1,:), block{b}));
            runs_in_block = session_block_array{pos_par,2}{2,pos_ses}{2,pos_block}(1,:);
            
            
            for w = 1:length(wave)
            % Inherited from previious structure and left because it makes things less complicated: 
            % only that w will always be 1
                
                
            % Find the filtered folders for each participant and session
            folders = dir([root_dir_bs '/data/' participant{p} '/']);
            
            file_names = {};
            for r = 1:length(runs_in_block)
                results =  contains({folders.name},'raw') & endsWith({folders.name},wave_name{w}) & (...
                    contains({folders.name},[participant{p}(1:4) session{s} '_' runs_in_block{r}]));
                infolder = find(results);

                if isempty(infolder) % If there are no coincidences, continue to next block/session
                    continue
                end

                file_names{r} = [participant{p} '/' folders(infolder).name '/data_0raw_' folders(infolder).name(5:end) '.mat'];
            end
            % Input files original from filtered files
            file_names = file_names(~cellfun('isempty', file_names')); % to avoid empty cells
            sFiles = file_names;
            
            if isempty(sFiles)
            error('for whatever reason, sFiles is empty, my friend...')
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Renaming events for templates %%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            disp(' ');      
            disp('-------------------------');
            disp(['Renaming LLR events for templates ' participant{p} '_' wave{w} '_' session{s} '_' block{b}]);
            disp(datetime)
            disp(' ');  

            % Process: Rename events
            for c =1:length(condition)
            sFiles = bst_process('CallProcess', 'process_evt_rename', sFiles, [], ...
                'src',  condition{c}, ...
                'dest', [wave{w} '_' condition{c} session{s} '_' block{b} '_temp_baseline']);
            end            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End renaming events %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Making epochs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            disp(' ');      
            disp('-------------------------');  
            disp(['Making LLR template epochs for ' participant{p} '_' wave{w} '_' session{s} '_' block{b}]);
            disp(datetime)
            disp(' ');  

            % Process: epoch data for baseline subtraction purposes
            sFiles = bst_process('CallProcess', 'process_import_data_event', sFiles, [], ...
                'subjectname',  participant{p}, ...
                'condition',    '', ...
            ...%    'datafile',     RawFiles, ...
                'eventname',    [wave{w} '_31' session{s} '_' block{b} '_temp_baseline,'...
                 wave{w} '_32' session{s} '_' block{b} '_temp_baseline,' wave{w} '_33' session{s} '_' block{b} '_temp_baseline,' wave{w} '_51' session{s} '_' block{b} '_temp_baseline,' ...
                 wave{w} '_52' session{s} '_' block{b} '_temp_baseline,' wave{w} '_53' session{s} '_' block{b} '_temp_baseline,' wave{w} '_71' session{s} '_' block{b} '_temp_baseline,' ...
                 wave{w} '_72' session{s} '_' block{b} '_temp_baseline,' wave{w} '_73' session{s} '_' block{b} '_temp_baseline,' wave{w} '_91' session{s} '_' block{b} '_temp_baseline,'...
                wave{w} '_92' session{s} '_' block{b} '_temp_baseline,' wave{w} '_93' session{s} '_' block{b} '_temp_baseline'], ...
                'timewindow',   [], ...
                'epochtime',    epoch_wave_bas{w}, ...
                'createcond',   1, ...
                'ignoreshort',  1, ...
                'channelalign', 0, ...
                'usectfcomp',   0, ...
                'usessp',       1, ...
                'freq',         resample_LLR, ...
                'baseline',     []); 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End making template epochs %%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Detrend %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if detrend_option_templates ==  1
            disp(' ');      
            disp('-------------------------');  
            disp(['Detrend LLR template epochs for ' participant{p} '_' wave{w} '_' session{s} '_' block{b}]);
            disp(datetime)
            disp(' ');  
            % Process: Remove linear trend: All file
            sFiles = bst_process('CallProcess', 'process_detrend', sFiles, [], ...
                'timewindow',  [], ...
                'sensortypes', '', ...
                'overwrite',   1);
            elseif detrend_option_templates ==  2 %   BVA option (to be completed)
            % Process: Remove linear trend: BVA option
            sFiles = bst_process('CallProcess', 'CNRL_BVA_process_detrend_LLR', sFiles, [], ...
                'timewindow',  epoch_wave{w}, ...
                'sensortypes', '', ...
                'overwrite',   1);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Baseline correction templates %%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            disp(' ');      
            disp('-------------------------');  
            disp(['Baseline correction for LLR templates ' participant{p} '_' wave{w} '_' session{s} '_' block{b}]);
            disp(datetime)
            disp(' ');

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Baseline correction first %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
            % Process: DC offset correction: [-50ms,1ms]
            sFiles = bst_process('CallProcess', 'process_baseline_norm', sFiles, [], ...
                'baseline',    epoch_baseline{w}, ...
                'sensortypes', '', ...
                'method',      'bl', ...  % DC offset correction:    x_std = x - &mu;
                'overwrite',   1);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
        end
    end
end

clearvars('-except', initialVars{:});
disp 'DONE WITH OBTAINING LLR ELONGATED EPOCHS FOR BASELINE CORRECTION (PART 1)!!!'
disp(datetime)
toc
end

%% Rename events in TEMPLATES back to their original names

if subtraction_approach ==  1
% From previous step they were renamed to create template folders,
% now we need to rename them back to create normal epoch folders
tic
disp(' ');      
disp('-------------------------');  
disp('RENAMING FILTERED FILES BACK TO ORIGINAL NAMES (LLR TEMPLATE)');  
disp(datetime)
disp('-------------------------');     
disp(' '); 

for p = 1:length(participant)
    
    disp(' ');      
    disp('-------------------------');
    disp(['loading participant ' participant{p}]);
    disp(datetime)
    disp(' '); 
    
    prot_subs = bst_get('ProtocolSubjects');
    current_sub = find(strcmp({prot_subs.Subject.Name}, participant{p}));
    db_reload_conditions(current_sub);

    % Define sessions within this participant
    pos_par = find(strcmp(session_block_array(:,1), participant{p}(1:4)));
    session = session_block_array{pos_par,2}(1,:);
    
    for s = 1:length(session)
        
        % Define blocks within this session
        pos_ses = find(strcmp(session_block_array{pos_par,2}(1,:), session{s}));
        block = session_block_array{pos_par,2}{2,pos_ses}(1,:);
        
        for b = 1:length(block)
            
            % Define runs within this block
            pos_block = find(strcmp(session_block_array{pos_par,2}{2,pos_ses}(1,:), block{b}));
            runs_in_block = session_block_array{pos_par,2}{2,pos_ses}{2,pos_block}(1,:);

            for w = 1:length(wave)
                
            % Find the filtered folders for each participant and session
            folders = dir([root_dir_bs '/data/' participant{p} '/']);
            
            file_names = {};
            for r = 1:length(runs_in_block)
                results =  contains({folders.name},'raw') & endsWith({folders.name},wave_name{w}) & (...
                    contains({folders.name},[participant{p}(1:4) session{s} '_' runs_in_block{r}]));
                infolder = find(results);

                if isempty(infolder) % If there are no coincidences, continue to next block/session
                    continue
                end

                file_names{r} = [participant{p} '/' folders(infolder).name '/data_0raw_' folders(infolder).name(5:end) '.mat'];
            end

            % Input files original from filtered files
            file_names = file_names(~cellfun('isempty', file_names')); % to avoid empty cells
            sFiles = file_names;
            
            if isempty(sFiles)
            error('for whatever reason, sFiles is empty, my friend...')
            end
            
            % Process: Rename events to its original form
            for c =1:length(condition)
                disp(' ');      
                disp('-------------------------');  
                disp(['Renaming template LLR back ' participant{p} '_' wave{w} '_' condition{c} session{s} '_' block{b}]);
                disp(datetime)
                disp(' ');  
                
            sFiles = bst_process('CallProcess', 'process_evt_rename', sFiles, [], ...
                'src',  [wave{w} '_' condition{c} session{s} '_' block{b} '_temp_baseline'], ...
                'dest', condition{c});
            end
            end
        end
    end
end

clearvars('-except', initialVars{:});
disp 'DONE WITH RENAMING BACK EVENTS FROM FILTERED DATAFILES (LLR TEMPLATE)!!!'
disp(datetime)
toc
end

%% Epoching NORMALLY (PART 1)
tic
disp(' ');      
disp('-------------------------');  
disp('EPOCHING NORMALLY (PART 1)');  
disp(datetime)
disp('-------------------------');     
disp(' ');

% Import filtered data this time (matlab_high and matlab_low folders)
for p = 1:length(participant)
    
    disp(' ');      
    disp('-------------------------');
    disp(['loading participant ' participant{p}]);
    disp(datetime)
    disp(' '); 
    
    prot_subs = bst_get('ProtocolSubjects');
    current_sub = find(strcmp({prot_subs.Subject.Name}, participant{p}));
    db_reload_conditions(current_sub);
    
    % Define sessions within this participant
    pos_par = find(strcmp(session_block_array(:,1), participant{p}(1:4)));
    session = session_block_array{pos_par,2}(1,:);
    
    for s = 1:length(session)
        % Define blocks within this session
        pos_ses = find(strcmp(session_block_array{pos_par,2}(1,:), session{s}));
        block = session_block_array{pos_par,2}{2,pos_ses}(1,:);
        
        for b = 1:length(block)
            
            % Define runs within this block
            pos_block = find(strcmp(session_block_array{pos_par,2}{2,pos_ses}(1,:), block{b}));
            runs_in_block = session_block_array{pos_par,2}{2,pos_ses}{2,pos_block}(1,:);
            
            for w = 1:length(wave) 
            % Inherited from previious structure and left because it makes things less complicated: 
            % only that w will always be 1
            % Find the filtered folders for each participant and session
            folders = dir([root_dir_bs '/data/' participant{p} '/']);
            
            file_names = {};
            for r = 1:length(runs_in_block)
                results =  contains({folders.name},'raw') & endsWith({folders.name},wave_name{w}) & (...
                    contains({folders.name},[participant{p}(1:4) session{s} '_' runs_in_block{r}]));
                infolder = find(results);

                if isempty(infolder) % If there are no coincidences, continue to next block/session
                    continue
                end

                file_names{r} = [participant{p} '/' folders(infolder).name '/data_0raw_' folders(infolder).name(5:end) '.mat'];
            end
            % Input files original from filtered files
            file_names = file_names(~cellfun('isempty', file_names')); % to avoid empty cells
            sFiles = file_names;

            if isempty(sFiles)
            error('for whatever reason, sFiles is empty, my friend...')
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Renaming normal epochs %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            disp(' ');      
            disp('-------------------------');
            disp(['Renaming events for normal epochs LLR ' participant{p} '_' wave{w} '_' session{s} '_' block{b}]);
            disp(datetime)
            disp(' ');  

            % Process: Rename events
            for c =1:length(condition)
            sFiles = bst_process('CallProcess', 'process_evt_rename', sFiles, [], ...
                'src',  condition{c}, ...
                'dest', [wave{w} '_' condition{c} session{s} '_' block{b} '_normal']);
            end            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End renaming events %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Making normal epochs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            disp(' ');      
            disp('-------------------------');  
            disp(['Making normal LLR epochs (part 1) for ' participant{p} '_' wave{w} '_' session{s} '_' block{b}]);
            disp(datetime)
            disp(' ');  
            
            % Process: epoch data for normal epochs
            sFiles = bst_process('CallProcess', 'process_import_data_event', sFiles, [], ...
                'subjectname',  participant{p}, ...
                'condition',    '', ...
            ...%    'datafile',     RawFiles, ...
                'eventname',    [wave{w} '_11' session{s} '_' block{b} '_normal,'...
                wave{w} '_12' session{s} '_' block{b} '_normal,' wave{w} '_13' session{s} '_' block{b} '_normal,' wave{w} '_31' session{s} '_' block{b} '_normal,' ...
                 wave{w} '_32' session{s} '_' block{b} '_normal,' wave{w} '_33' session{s} '_' block{b} '_normal,' wave{w} '_51' session{s} '_' block{b} '_normal,' ...
                 wave{w} '_52' session{s} '_' block{b} '_normal,' wave{w} '_53' session{s} '_' block{b} '_normal,' wave{w} '_71' session{s} '_' block{b} '_normal,' ...
                 wave{w} '_72' session{s} '_' block{b} '_normal,' wave{w} '_73' session{s} '_' block{b} '_normal,' wave{w} '_91' session{s} '_' block{b} '_normal,'...
                wave{w} '_92' session{s} '_' block{b} '_normal,' wave{w} '_93' session{s} '_' block{b} '_normal'], ...
                'timewindow',   [], ...
                'epochtime',    epoch_wave{w}, ...
                'createcond',   1, ...
                'ignoreshort',  1, ...
                'channelalign', 0, ...
                'usectfcomp',   0, ...
                'usessp',       1, ...
                'freq',         resample_LLR, ...
                'baseline',     []); 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End making normal epochs %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Detrend %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if detrend_option ==  1     
            disp(' ');      
            disp('-------------------------');  
            disp(['Detrend normal LLR epochs (part 1) for ' participant{p} '_' wave{w} '_' session{s} '_' block{b}]);
            disp(datetime)
            disp(' ');  
                
            % Process: Remove linear trend: All file
            sFiles = bst_process('CallProcess', 'process_detrend', sFiles, [], ...
                'timewindow',  [], ...
                'sensortypes', '', ...
                'overwrite',   1);
            elseif detrend_option ==  2 %   BVA option
            % Process: Remove linear trend: BVA option
            sFiles = bst_process('CallProcess', 'CNRL_BVA_process_detrend_LLR', sFiles, [], ...
                'timewindow',  epoch_wave{w}, ...
                'sensortypes', '', ...
                'overwrite',   1);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Baseline correction first %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
            disp(' ');      
            disp('-------------------------');  
            disp(['Baseline correcting normal LLR epochs (part 1) for ' participant{p} '_' wave{w} '_' session{s} '_' block{b}]);
            disp(datetime)
            disp(' ');  
            
            % Process: DC offset correction: [-50ms,1ms]
            sFiles = bst_process('CallProcess', 'process_baseline_norm', sFiles, [], ...
                'baseline',    epoch_baseline{w}, ...
                'sensortypes', '', ...
                'method',      'bl', ...  % DC offset correction:    x_std = x - &mu;
                'overwrite',   1);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
        end
    end
end

clearvars('-except', initialVars{:});
disp 'DONE WITH NORMAL EPOCHING (PART 1)!!!'
disp(datetime)
toc

%% Rename NORMAL epochs names back to original
tic
disp(' ');      
disp('-------------------------');  
disp('RENAMING FILTERED FILES BACK TO ORIGINAL NAMES (NORMAL EPOCHS)');  
disp(datetime)
disp('-------------------------');     
disp(' '); 

for p = 1:length(participant)
    
    disp(' ');      
    disp('-------------------------');
    disp(['loading participant ' participant{p}]);
    disp(datetime)
    disp(' '); 
    
    prot_subs = bst_get('ProtocolSubjects');
    current_sub = find(strcmp({prot_subs.Subject.Name}, participant{p}));
    db_reload_conditions(current_sub);

    % Define sessions within this participant
    pos_par = find(strcmp(session_block_array(:,1), participant{p}(1:4)));
    session = session_block_array{pos_par,2}(1,:);
    
    for s = 1:length(session)
        
        % Define blocks within this session
        pos_ses = find(strcmp(session_block_array{pos_par,2}(1,:), session{s}));
        block = session_block_array{pos_par,2}{2,pos_ses}(1,:);
        
        for b = 1:length(block)
            
            % Define runs within this block
            pos_block = find(strcmp(session_block_array{pos_par,2}{2,pos_ses}(1,:), block{b}));
            runs_in_block = session_block_array{pos_par,2}{2,pos_ses}{2,pos_block}(1,:);
            
            for w = 1:length(wave)

            % Find the filtered folders for each participant and session
            folders = dir([root_dir_bs '/data/' participant{p} '/']);
            
            file_names = {};
            for r = 1:length(runs_in_block)
                results =  contains({folders.name},'raw') & endsWith({folders.name},wave_name{w}) & (...
                    contains({folders.name},[participant{p}(1:4) session{s} '_' runs_in_block{r}]));
                infolder = find(results);

                if isempty(infolder) % If there are no coincidences, continue to next block/session
                    continue
                end

                file_names{r} = [participant{p} '/' folders(infolder).name '/data_0raw_' folders(infolder).name(5:end) '.mat'];
            end

            % Input files original from filtered files
            file_names = file_names(~cellfun('isempty', file_names')); % to avoid empty cells
            sFiles = file_names;

            if isempty(sFiles)
            error('for whatever reason, sFiles is empty, my friend...')
            end
            
            % Process: Rename events to its original form
            for c =1:length(condition)
                disp(' ');      
                disp('-------------------------');  
                disp(['Renaming back ' participant{p} '_' wave{w} '_' condition{c} session{s} '_' block{b}]);
                disp(datetime)
                disp(' ');  
                
            sFiles = bst_process('CallProcess', 'process_evt_rename', sFiles, [], ...
                'src',  [wave{w} '_' condition{c} session{s} '_' block{b} '_normal'], ...
                'dest', condition{c});
            end
            end
        end
    end
end

clearvars('-except', initialVars{:});
disp 'DONE WITH RENAMING BACK EVENTS FROM FILTERED DATAFILES (NORMAL EPOCHS)!!!'
disp(datetime)
toc

%% Make epochs for noise covariance analyses (Part 1)
if covar_epochs == 1

tic
disp(' ');      
disp('-------------------------');  
disp('OBTAINING EPOCHS FOR COVARIANCE ANALYSES (Part 1)');  
disp(datetime)
disp('-------------------------');     
disp(' ');

% Import filtered data this time (matlab_high and matlab_low folders)
for p = 1:length(participant)
    
    disp(' ');      
    disp('-------------------------');
    disp(['loading participant ' participant{p}]);
    disp(datetime)
    disp(' '); 
    
    prot_subs = bst_get('ProtocolSubjects');
    current_sub = find(strcmp({prot_subs.Subject.Name}, participant{p}));
    db_reload_conditions(current_sub);
    
    % Define sessions within this participant
    pos_par = find(strcmp(session_block_array(:,1), participant{p}(1:4)));
    session = session_block_array{pos_par,2}(1,:);
    
    for s = 1:length(session)
        
        % Define blocks within this session
        pos_ses = find(strcmp(session_block_array{pos_par,2}(1,:), session{s}));
        block = session_block_array{pos_par,2}{2,pos_ses}(1,:);
        
        for b = 1:length(block)
            
            % Define runs within this block
            pos_block = find(strcmp(session_block_array{pos_par,2}{2,pos_ses}(1,:), block{b}));
            runs_in_block = session_block_array{pos_par,2}{2,pos_ses}{2,pos_block}(1,:);
            
            for w = 1:length(wave)  
            % Find the filtered folders for each participant and session
            folders = dir([root_dir_bs '/data/' participant{p} '/']);
            
            file_names = {};
            for r = 1:length(runs_in_block)
                results =  contains({folders.name},'raw') & endsWith({folders.name},wave_name{w}) & (...
                    contains({folders.name},[participant{p}(1:4) session{s} '_' runs_in_block{r}]));
                infolder = find(results);

                if isempty(infolder) % If there are no coincidences, continue to next block/session
                    continue
                end

                file_names{r} = [participant{p} '/' folders(infolder).name '/data_0raw_' folders(infolder).name(5:end) '.mat'];
            end
            % Input files original from filtered files
            file_names = file_names(~cellfun('isempty', file_names')); % to avoid empty cells
            sFiles = file_names;

            if isempty(sFiles)
            error('for whatever reason, sFiles is empty, my friend...')
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%% Deleting events too close to boundaries %%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            disp(' ');      
            disp('-------------------------');
            disp(['Deleting events too close to boundaryEnd (cov) for ' participant{p} '_' wave{w} '_' session{s} '_' block{b}]);
            disp(datetime)
            disp(' ');  
            
            for c = 13:length(condition) % JUST THE 91, 92 AND 93
            sFiles = bst_process('CallProcess', 'process_evt_remove_simult', sFiles, [], ...
                'remove', condition{c}, ...
                'target', 'boundaryEnd', ...
                'dt',     3.5, ...
                'rename', 0);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Renaming events for covariance %%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            disp(' ');      
            disp('-------------------------');
            disp(['Renaming events for covariance ' participant{p} '_' wave{w} '_' session{s} '_' block{b}]);
            disp(datetime)
            disp(' ');  

            % Process: Rename events
            for c = 13:length(condition) % JUST THE 91, 92 AND 93
            sFiles = bst_process('CallProcess', 'process_evt_rename', sFiles, [], ...
                'src',  condition{c}, ...
                'dest', [wave{w} '_' condition{c} session{s} '_' block{b} '_covariance']);
            end            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End renaming events %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Making epochs (covariance) %%%%%%%%%%%%%%%%%%%%%%%%   
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            disp(' ');      
            disp('-------------------------');  
            disp(['Making covariance epochs for ' participant{p} '_' wave{w} '_' session{s} '_' block{b}]);
            disp(datetime)
            disp(' ');  

            % Process: epoch data for covariance purposes
            sFiles = bst_process('CallProcess', 'process_import_data_event', sFiles, [], ...
                'subjectname',  participant{p}, ...
                'condition',    '', ...
            ...%    'datafile',     RawFiles, ...
                'eventname',    [wave{w} '_91' session{s} '_' block{b} '_covariance,'...
                wave{w} '_92' session{s} '_' block{b} '_covariance,' wave{w} '_93' session{s} '_' block{b} '_covariance'], ...
                'timewindow',   [], ...
                'epochtime',    epoch_covar{w}, ...
                'createcond',   1, ...
                'ignoreshort',  1, ...
                'channelalign', 0, ...
                'usectfcomp',   0, ...
                'usessp',       1, ...
                'freq',         resample_LLR, ... % CHECK IF WE HAVE ENOUGH DATAPOINTS WITH THIS...
                'baseline',     []); 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End making covariance epochs %%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Baseline correction of covariance epochs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
            disp(' ');      
            disp('-------------------------');  
            disp(['Baseline correcting covariance epochs for ' participant{p} '_' wave{w} '_' session{s} '_' block{b}]);
            disp(datetime)
            disp(' ');  
            
            % Process: DC offset correction: [-150ms,1ms]
            sFiles = bst_process('CallProcess', 'process_baseline_norm', sFiles, [], ...
                'baseline',    epoch_baseline{w}, ...
                'sensortypes', '', ...
                'method',      'bl', ...  % DC offset correction:    x_std = x - &mu;
                'overwrite',   1);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            end
        end
    end
end
            
clearvars('-except', initialVars{:});
disp 'DONE WITH OBTAINING EPOCHS FOR COVARIANCE ANALYSES (Part 1)!!!'
disp(datetime)
toc

end      
       
%% Rename covariance epochs back to original

if covar_epochs == 1

tic
disp(' ');      
disp('-------------------------');  
disp('RENAMING FILTERED FILES BACK TO ORIGINAL NAMES (COVARIANCE EPOCHS)');  
disp(datetime)
disp('-------------------------');     
disp(' '); 

for p = 1:length(participant)
    
    disp(' ');      
    disp('-------------------------');
    disp(['loading participant ' participant{p}]);
    disp(datetime)
    disp(' '); 
    
    prot_subs = bst_get('ProtocolSubjects');
    current_sub = find(strcmp({prot_subs.Subject.Name}, participant{p}));
    db_reload_conditions(current_sub);
    
    % Define sessions within this participant
    pos_par = find(strcmp(session_block_array(:,1), participant{p}(1:4)));
    session = session_block_array{pos_par,2}(1,:);
    
    for s = 1:length(session)
        
        % Define blocks within this session
        pos_ses = find(strcmp(session_block_array{pos_par,2}(1,:), session{s}));
        block = session_block_array{pos_par,2}{2,pos_ses}(1,:);
        
        for b = 1:length(block)
            
            % Define runs within this block
            pos_block = find(strcmp(session_block_array{pos_par,2}{2,pos_ses}(1,:), block{b}));
            runs_in_block = session_block_array{pos_par,2}{2,pos_ses}{2,pos_block}(1,:);
            
            for w = 1:length(wave)  
            % Find the filtered folders for each participant and session
            folders = dir([root_dir_bs '/data/' participant{p} '/']);
            
            file_names = {};
            for r = 1:length(runs_in_block)
                results =  contains({folders.name},'raw') & endsWith({folders.name},wave_name{w}) & (...
                    contains({folders.name},[participant{p}(1:4) session{s} '_' runs_in_block{r}]));
                infolder = find(results);

                if isempty(infolder) % If there are no coincidences, continue to next block/session
                    continue
                end

                file_names{r} = [participant{p} '/' folders(infolder).name '/data_0raw_' folders(infolder).name(5:end) '.mat'];
            end
            % Input files original from filtered files
            file_names = file_names(~cellfun('isempty', file_names')); % to avoid empty cells
            sFiles = file_names;
    
            if isempty(sFiles)
            error('for whatever reason, sFiles is empty, my friend...')
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%% Renaming events back from covariance %%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            disp(' ');      
            disp('-------------------------');
            disp(['Renaming events back from covariance ' participant{p} '_' wave{w} '_' session{s} '_' block{b}]);
            disp(datetime)
            disp(' ');  

            % Process: Rename events
            for c = 13:length(condition)
            sFiles = bst_process('CallProcess', 'process_evt_rename', sFiles, [], ...
                'src',  [wave{w} '_' condition{c} session{s} '_' block{b} '_covariance'], ...
                'dest', condition{c});
            end            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%% End renaming events back %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
        end
    end
end

clearvars('-except', initialVars{:});
disp 'DONE RENAMING FILTERED FILES BACK TO ORIGINAL NAMES (COVARIANCE EPOCHS)!!!'
disp(datetime)
toc

end

%% FROM HERE ON, COPY/MOVE FILES AND CHANGE PROTOCOL NAME

tic
disp(' ');      
disp('-------------------------');  
disp('MOVING LLR EPOCHED FILES TO SEGMENTED PROTOCOL');  
disp(datetime)
disp('-------------------------');     
disp(' '); 

new_root_dir_bs = '~/brainstorm_db/Project_LLRseg'; % depends on whether is LLR or MLR

for par = 1:length(participant_general_list)
    
    disp(' ');      
    disp('-------------------------');
    disp(['Moving LLR epoched files to segmented protocol for participant ' participant_general_list{par}]);
    disp(datetime)
    disp(' ');  
    
    % COPY default study info
    if isfile([root_dir_bs '/data/' participant_general_list{par} '/@default_study/brainstormstudy.mat'])
        copyfile([root_dir_bs '/data/' participant_general_list{par} '/@default_study/brainstormstudy.mat'], [new_root_dir_bs '/data/' participant_general_list{par} '/@default_study/'])
    end
    % COPY intra study info (at this point there should be nothing else there)
    if isfile([root_dir_bs '/data/' participant_general_list{par} '/@intra/brainstormstudy.mat'])
        copyfile([root_dir_bs '/data/' participant_general_list{par} '/@intra/brainstormstudy.mat'], [new_root_dir_bs '/data/' participant_general_list{par} '/@intra/'])
    end
    
    % MOVE epoched files
    folders = dir([root_dir_bs '/data/' participant_general_list{par} '/']);
    results = startsWith({folders.name},'LLR'); % ALL EPOCHED START WITH THAT
    infolder = find(results);
    for l = 1:length(infolder) % for each coincidence
        line = infolder(l);
        movefile([root_dir_bs '/data/' participant_general_list{par} '/' folders(line).name],[new_root_dir_bs '/data/' participant_general_list{par} '/']);
    end
end

% CHANGE BRAINSTORM PROTOCOL PARAMETERS AND LOAD NEW PROTOCOL
bst_set('BrainstormDbDir',BrainstormDbDir)
% Select the correct protocol
ProtocolName = 'Project_LLRseg'; % Enter the name of your protocol
sProtocol.Comment = ProtocolName;
% sProtocol.SUBJECTS = [home 'anat'];
sProtocol.SUBJECTS = '~/brainstorm_db/Project_LLRseg/anat';
sProtocol.STUDIES = '~/brainstorm_db/Project_LLRseg/data';
db_edit_protocol('load',sProtocol);
% Get the protocol index
iProtocol = bst_get('Protocol', ProtocolName);
if isempty(iProtocol)
    error(['Unknown protocol: ' ProtocolName]);
end
% Select the current procotol
gui_brainstorm('SetCurrentProtocol', iProtocol);

% Change protocol name from here on
root_dir_bs = new_root_dir_bs; % FROM NOW ON, IT IS SEGMENTED DATA

clearvars('-except', initialVars{:}); % FROM HERE ON THOUGH, root_dir_bs will be changed to the segmented one
disp 'DONE MOVING LLR EPOCHED FILES TO SEGMENTED PROTOCOL!!!'
disp(datetime)
toc

%% Run all sections for 1 subject at a time (except if they belong to diff protocols)

% If willing to do so, put the sections in this loop
for par = 1:length(participant_general_list) 
    participant = {}; % rebuild the participant variable everytime
    participant{1} = participant_general_list{par};
    %%%%%%%%%%% CODE SECTIONS HERE %%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%% Epoching to create TEMPLATES (Part 2)

if subtraction_approach ==  1
tic
disp(' ');      
disp('-------------------------');  
disp('OBTAINING ELONGATED LLR EPOCHS FOR BASELINE CORRECTION (PART 2)'); 
disp(datetime)
disp('-------------------------');     
disp(' ');

for p = 1:length(participant)
    
    disp(' ');      
    disp('-------------------------');
    disp(['loading participant ' participant{p}]);
    disp(datetime)
    disp(' '); 
    
    prot_subs = bst_get('ProtocolSubjects');
    current_sub = find(strcmp({prot_subs.Subject.Name}, participant{p}));
    db_reload_conditions(current_sub);
    
    % Define sessions within this participant
    pos_par = find(strcmp(session_block_array(:,1), participant{p}(1:4)));
    session = session_block_array{pos_par,2}(1,:);
    
    for s = 1:length(session)
        
        % Define blocks within this session
        pos_ses = find(strcmp(session_block_array{pos_par,2}(1,:), session{s}));
        block = session_block_array{pos_par,2}{2,pos_ses}(1,:);
        
        for b = 1:length(block)
            for c = 4:length(condition) % ONLY 31 ONWARDS
                for w = 1:length(wave)
                    
                    folders = dir([root_dir_bs '/data/' participant{p} '/']);
                    results = contains({folders.name},[wave{w} '_' condition{c} session{s} '_' block{b}]) & contains({folders.name},'_temp_baseline');
                    infolder = find(results);
                    
                    if isempty(infolder)
                        continue
                    end
                    if size(infolder,2)> 1 % in case of more than one coincidence, error
                        error('More than one coincidence');
                    end
                        
                    file_names = {};
                    dir_sweeps = dir([root_dir_bs '/data/' participant{p} '/' folders(infolder).name]);
                    sweeps_norm = contains({dir_sweeps.name},'_trial');
                    sweeps_list = find(sweeps_norm);
                    for j= 1:length(sweeps_list)
                        position = sweeps_list(j);
                        file_names{j} = [participant{p} '/' folders(infolder).name '/' dir_sweeps(position).name];
                    end
                sFiles = file_names; % all sweeps from that block and condition contained here now
                
                if isempty(sFiles)
                    % It can be that sometimes certain sweeps types (2235
                    % 52 C b4 do not exist while others do, because for
                    % whatever reason there were no 52 in a block). Since
                    % it has happened, better to continue here
                	continue
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Amplitude threshold normal epochs %%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if (sensor_analysis == 1) || (sensor_analysis == 4) % EEG only
                
                disp(' ');      
                disp('-------------------------');  
                disp(['Cleaning LLR template epochs for ' participant{p} '_' wave{w} '_' condition{c} session{s} '_' block{b} '(EEG)']);
                disp(datetime)
                disp(' '); 
                
                % Process: Detect bad trials: Absolute threshold EEG
                sFiles_EEG = bst_process('CallProcess', 'process_CNRL_detectbad', sFiles, [], ...
                    'timewindow', [], ...
                    'meggrad',    [0, 0], ...
                    'megmag',     [0, 0], ...
                    'eeg',        reject_EEG_absolute, ...
                    'ieeg',       [0, 0], ...
                    'eog',        [0, 0], ...
                    'ecg',        [0, 0], ...
                    'rejectmode', 2);  % Reject the entire trial
                
                % save and reset bad trials
                trials_file = [root_dir_bs '/data/' participant{p} '/' folders(infolder).name '/brainstormstudy.mat'];
                eval(['load ' trials_file]); % load them
                % save current ones in events folder
                if ~exist([root_dir '/Events/BadTrials_Template/' ProtocolName(7:9) '/EEG/' participant{p} '/' folders(infolder).name], 'dir')
                    mkdir([root_dir '/Events/BadTrials_Template/'], [ProtocolName(7:9) '/EEG/' participant{p} '/' folders(infolder).name]);
                end
                save([root_dir '/Events/BadTrials_Template/' ProtocolName(7:9) '/EEG/' participant{p} '/' folders(infolder).name '/BadTrials.mat'],'BadTrials'); % save bad trials
                % reset them and save back to original file
                BadTrials = cell([], 1); % exact structure it has when empty
                save(trials_file,'BadTrials', 'DateOfStudy', 'Name');
                end
                
                % BadTrials, DateOfStudy, Name
                
                if (sensor_analysis == 2) || (sensor_analysis == 4) % MEG only
                    
                disp(' ');      
                disp('-------------------------');  
                disp(['Cleaning LLR template epochs for ' participant{p} '_' wave{w} '_' condition{c} session{s} '_' block{b} '(MEG)']);
                disp(datetime)
                disp(' '); 
                
                % Process: Detect bad trials: Peak-to-peak  MEG ONLY
                sFiles_MEG = bst_process('CallProcess', 'process_CNRL_detectbad',sFiles, [], ...
                    'timewindow', [], ...
                    'meggrad',    reject_MEG_GRAD_absolute, ...
                    'megmag',     reject_MEG_MAG_absolute, ...
                    'eeg',        [0, 0], ...
                    'ieeg',       [0, 0], ...
                    'eog',        [0, 0], ...
                    'ecg',        [0, 0], ...
                    'rejectmode', 2);  % Reject the entire trial
                
                % save and reset bad trials
                trials_file = [root_dir_bs '/data/' participant{p} '/' folders(infolder).name '/brainstormstudy.mat'];
                eval(['load ' trials_file]); % load them
                % save current ones in events folder
                if ~exist([root_dir '/Events/BadTrials_Template/' ProtocolName(7:9) '/MEG/' participant{p} '/' folders(infolder).name], 'dir')
                    mkdir([root_dir '/Events/BadTrials_Template/'], [ProtocolName(7:9) '/MEG/' participant{p} '/' folders(infolder).name]);
                end
                save([root_dir '/Events/BadTrials_Template/' ProtocolName(7:9) '/MEG/' participant{p} '/' folders(infolder).name '/BadTrials.mat'],'BadTrials'); % save bad trials
                % reset them and save back to original file
                BadTrials = cell([], 1); % exact structure it has when empty
                save(trials_file,'BadTrials', 'DateOfStudy', 'Name');
                end
                
                if (sensor_analysis == 3) || (sensor_analysis == 4) % Combined EEG and MEG
                
                disp(' ');      
                disp('-------------------------');  
                disp(['Cleaning LLR template epochs for ' participant{p} '_' wave{w} '_' condition{c} session{s} '_' block{b} '(both_mod)']);
                disp(datetime)
                disp(' '); 
                
                % Process: Detect bad trials: Peak-to-peak  BOTH EEG AND MEG
                sFiles_both_mod = bst_process('CallProcess', 'process_CNRL_detectbad',sFiles, [], ...
                    'timewindow', [], ...
                    'meggrad',    reject_MEG_GRAD_absolute, ...
                    'megmag',     reject_MEG_MAG_absolute, ...
                    'eeg',        reject_EEG_absolute, ...
                    'ieeg',       [0, 0], ...
                    'eog',        [0, 0], ...
                    'ecg',        [0, 0], ...
                    'rejectmode', 2);  % Reject the entire trial
                
                % save and reset bad trials
                trials_file = [root_dir_bs '/data/' participant{p} '/' folders(infolder).name '/brainstormstudy.mat'];
                eval(['load ' trials_file]); % load them
                % save current ones in events folder
                if ~exist([root_dir '/Events/BadTrials_Template/' ProtocolName(7:9) '/both_mod/' participant{p} '/' folders(infolder).name], 'dir')
                    mkdir([root_dir '/Events/BadTrials_Template/'], [ProtocolName(7:9) '/both_mod/' participant{p} '/' folders(infolder).name]);
                end
                save([root_dir '/Events/BadTrials_Template/' ProtocolName(7:9) '/both_mod/' participant{p} '/' folders(infolder).name '/BadTrials.mat'],'BadTrials'); % save bad trials
                % reset them and save back to original file
                BadTrials = cell([], 1); % exact structure it has when empty
                save(trials_file,'BadTrials', 'DateOfStudy', 'Name');
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Average of sweeps normal epochs %%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if (sensor_analysis == 1) || (sensor_analysis == 4) % EEG only
                    
                % SENSOR AVERAGE EEG (WITH INTERPOLATION)   
                disp(' ');      
                disp('-------------------------');  
                disp(['Averaging LLR template epochs for ' participant{p} '_' wave{w} '_' condition{c} session{s} '_' block{b} '(EEG sensor)']);
                disp(datetime)
                disp(' '); 
                
                sFiles_EEG_averaged_sensor = bst_process('CallProcess', 'process_average', sFiles_EEG, [], ...
                    'avgtype',         5, ...  % By trial group (folder average)
                    'avg_func',        1, ...  % Arithmetic average:  mean(x)
                    'weighted',        0, ...
                    'keepevents', 0);
                
                % Interpolate channels
                sFiles_EEG_averaged_sensor = bst_process('CallProcess', 'process_eeg_interpbad', sFiles_EEG_averaged_sensor, [], ...
                'maxdist',     5, ...
                'sensortypes', 'EEG', ...
                'overwrite',   1);

                % Process: Add tag
                sFiles_EEG_averaged_sensor = bst_process('CallProcess', 'process_add_tag', sFiles_EEG_averaged_sensor, [], ...
                    'tag',           'EEG_average_sensor', ...
                    'output',        2);  % Add to file name (1 to add a tag)

                % Process: Set name
                sFiles_EEG_averaged_sensor = bst_process('CallProcess', 'process_set_comment', sFiles_EEG_averaged_sensor, [], ...
                    'tag',           'EEG_average_sensor', ...
                    'isindex',       1);
                
                end
                
                if (sensor_analysis == 2) || (sensor_analysis == 4) % MEG only
                    
                % SOURCE AVERAGE MEG (NO INTERPOLATION)                  
                disp(' ');      
                disp('-------------------------');  
                disp(['Averaging LLR template epochs for ' participant{p} '_' wave{w} '_' condition{c} session{s} '_' block{b} '(MEG sensor)']);
                disp(datetime)
                disp(' '); 
                
                sFiles_MEG_averaged_sensor = bst_process('CallProcess', 'process_average', sFiles_MEG, [], ...
                    'avgtype',         5, ...  % By trial group (folder average)
                    'avg_func',        1, ...  % Arithmetic average:  mean(x)
                    'weighted',        0, ...
                    'keepevents', 0);
                
                % FIND ONE THAT WORKS FOR MEG
%                 % Interpolate channels
%                 sFiles_MEG_averaged_sensor = bst_process('CallProcess', 'process_eeg_interpbad', sFiles_MEG_averaged_sensor, [], ...
%                 'maxdist',     5, ...
%                 'sensortypes', 'MEG', ...
%                 'overwrite',   1);

                % Process: Add tag
                sFiles_MEG_averaged_sensor = bst_process('CallProcess', 'process_add_tag', sFiles_MEG_averaged_sensor, [], ...
                    'tag',           'MEG_average_sensor', ...
                    'output',        2);  % Add to file name (1 to add a tag)

                % Process: Set name
                sFiles_MEG_averaged_sensor = bst_process('CallProcess', 'process_set_comment', sFiles_MEG_averaged_sensor, [], ...
                    'tag',           'MEG_average_sensor', ...
                    'isindex',       1);
                end
                
                if (sensor_analysis == 3) || (sensor_analysis == 4) % Combined EEG and MEG
                    
                % SOURCE AVERAGE BOTH MOD (NO INTERPOLATION)                  
                disp(' ');      
                disp('-------------------------');  
                disp(['Averaging LLR template epochs for ' participant{p} '_' wave{w} '_' condition{c} session{s} '_' block{b} '(both mod sensor)']);
                disp(datetime)
                disp(' '); 
                
                sFiles_both_mod_averaged_sensor = bst_process('CallProcess', 'process_average', sFiles_both_mod, [], ...
                    'avgtype',         5, ...  % By trial group (folder average)
                    'avg_func',        1, ...  % Arithmetic average:  mean(x)
                    'weighted',        0, ...
                    'keepevents', 0);
                
                % FIND ONE THAT WORKS FOR MEG and EEG combined
%                 % Interpolate channels
%                 sFiles_MEG_averaged_sensor = bst_process('CallProcess', 'process_eeg_interpbad', sFiles_MEG_averaged_sensor, [], ...
%                 'maxdist',     5, ...
%                 'sensortypes', 'MEG', ...
%                 'overwrite',   1);

                % Process: Add tag
                sFiles_both_mod_averaged_sensor = bst_process('CallProcess', 'process_add_tag', sFiles_both_mod_averaged_sensor, [], ...
                    'tag',           'both_mod_average_sensor', ...
                    'output',        2);  % Add to file name (1 to add a tag)

                % Process: Set name
                sFiles_both_mod_averaged_sensor = bst_process('CallProcess', 'process_set_comment', sFiles_both_mod_averaged_sensor, [], ...
                    'tag',           'both_mod_average_sensor', ...
                    'isindex',       1);
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End Average of sweeps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                end
            end
        end
    end
end

clearvars('-except', initialVars{:});
disp 'DONE WITH OBTAINING LLR ELONGATED EPOCHS FOR BASELINE CORRECTION (PART 2)!!!'
disp(datetime)
toc

end

%% Average TEMPLATES across blocks

if subtraction_approach ==  1
% In case data from every block is insufficient for substraction approach
tic
disp(' ');      
disp('-------------------------');  
disp('AVERAGING LLR TEMPLATES ACROSS BLOCKS (IN CASE BLOCK DATA IS INSUFFICIENT FOR SUBSTRACTION APPROACH)'); 
disp(datetime)
disp('-------------------------');
disp(' ');

if sensor_analysis == 1
    modality_data = {'EEG'};
elseif sensor_analysis == 2
    modality_data = {'MEG'};
elseif sensor_analysis == 3
    modality_data = {'both_mod'};
end

for p = 1:length(participant)
    
    disp(' ');      
    disp('-------------------------');
    disp(['loading participant ' participant{p}]);
    disp(datetime)
    disp(' '); 
    
    prot_subs = bst_get('ProtocolSubjects');
    current_sub = find(strcmp({prot_subs.Subject.Name}, participant{p}));
    db_reload_conditions(current_sub);
    
    % Define sessions within this participant
    pos_par = find(strcmp(session_block_array(:,1), participant{p}(1:4)));
    session = session_block_array{pos_par,2}(1,:);
    
    for s = 1:length(session)
        
        % Define blocks within this session
        pos_ses = find(strcmp(session_block_array{pos_par,2}(1,:), session{s}));
        block = session_block_array{pos_par,2}{2,pos_ses}(1,:);
        
        for w = 1:length(wave)
            for c = 4:length(condition) % ONLY 31 ONWARDS
                folders = dir([root_dir_bs '/data/' participant{p} '/']);
                results =  contains({folders.name},[wave{w} '_' condition{c} session{s}]) & contains({folders.name},'_temp_baseline');
                infolder = find(results); % So get all blocks you find for that session
                if isempty(infolder)
                    continue
                    % error(['no ' wave{w} '_' condition{c} session{s} '_normal blocks for ' participant{p}]);
                end
                for mode = 1:length(modality_data) % EEG, MEG and both mod
                    file_names = {};
                    for l= 1:length(infolder) % for each block found
                        line = infolder(l);
                        % Check if that block is in the block list, if not, go to next block
                        init_block_name = regexp(folders(line).name, '_b', 'once');
                        fol_nam = contains(block, folders(line).name(init_block_name+1:init_block_name+2));
                        in_fol_nam = find(fol_nam); % index of position in logical
                        if isempty(in_fol_nam); continue; end % the block in the folder name is not present in the list, so go to next
                        % End checking if block exist in list
                        sub_folders = dir([root_dir_bs '/data/' participant{p} '/' folders(line).name '/']);
                        sub_results =  contains({sub_folders.name},[modality_data{mode} '_average_sensor']); % Only those averages that are interpolated
                        sub_infolder = find(sub_results);
                        if isempty(sub_infolder) % if average is not present despite the folder existing
                            continue
                        end
                        file_names{l} = [participant{p} '/' folders(line).name '/' sub_folders(sub_infolder).name]; % sub_results (using logicals)
                    end
                    file_names = file_names(~cellfun('isempty', file_names')); % to avoid empty cells
                    sFiles = file_names;

                    if isempty(sFiles)
                        % It could be that all blocks from a particular
                        % modality have no average because of rejected
                        % trials in all blocks of the session, in which
                        % casse better to continue
                        continue
                    end
                    
                    disp(' ');      
                    disp('-------------------------');  
                    disp(['Averaging Template epochs for ' participant{p} '_' modality_data{mode} '_' wave{w} '_' condition{c} session{s}]);
                    disp(datetime)
                    disp(' ');

                    % If stated, find and delete previous template averages across
                    if delete_previous_file == 1
                        % check if there is already a template average for
                        % each session in intra folder and delete it
                        folders_delete = dir([root_dir_bs '/data/' participant{p} '/@intra/']);
                        % data_LLR_average_200727_1504_Template_both_mod_LLR_51
                        results_delete =  endsWith({folders_delete.name}, ['Template_' modality_data{mode} '_' wave{w} '_' condition{c} session{s} '.mat']); 
                        infolder_delete = find(results_delete);
                        if ~isempty(infolder_delete) % file exists, therefore delete it
                           delete([root_dir_bs '/data/' participant{p} '/@intra/' folders_delete(infolder_delete).name]);
                        end
                    end
                    
                    % Process: Weighted Average: Everything
                    sFiles = bst_process('CallProcess', 'process_average', sFiles, [], ...
                        'avgtype',       1, ...  % Everything
                        'avg_func',      1, ...  % Arithmetic average:  mean(x)
                        'weighted',      1, ...
                        'keepevents',    0);

                    % Process: Add tag
                    sFiles = bst_process('CallProcess', 'process_add_tag', sFiles, [], ...
                        'tag',           ['Template_' modality_data{mode} '_' wave{w} '_' condition{c} session{s}], ...
                        'output',        2);  % Add to file name (1 to add a tag)

                    % Process: Set name
                    sFiles = bst_process('CallProcess', 'process_set_comment', sFiles, [], ...
                        'tag',           ['Template_' modality_data{mode} '_' wave{w} '_' condition{c} session{s}], ...
                        'isindex',       1);
                end
            end
        end
    end
end

% reset modality_data variable to its original
modality_data = {'EEG','MEG','both_mod'};

clearvars('-except', initialVars{:});
disp 'DONE WITH AVERAGING LLR TEMPLATES ACROSS BLOCKS!!!'
disp(datetime)
toc
end

%% Average TEMPLATES across blocks from all sessions

if subtraction_approach ==  1
    
tic
disp(' ');      
disp('-------------------------');  
disp('AVERAGING LLR WAVEFORMS ACROSS BLOCKS FROM ALL SESSIONS (TEMPLATE EPOCHS)'); 
disp(datetime)
disp('-------------------------');
disp(' ');

if sensor_analysis == 1
    modality_data = {'EEG'};
elseif sensor_analysis == 2
    modality_data = {'MEG'};
elseif sensor_analysis == 3
    modality_data = {'both_mod'};
end

for p = 1:length(participant)
    
    disp(' ');      
    disp('-------------------------');
    disp(['loading participant ' participant{p}]);
    disp(datetime)
    disp(' '); 

    % reload subject normally
    prot_subs = bst_get('ProtocolSubjects');
    current_sub = find(strcmp({prot_subs.Subject.Name}, participant{p}));
    db_reload_conditions(current_sub);
    
    % Define sessions within this participant
    pos_par = find(strcmp(session_block_array(:,1), participant{p}(1:4)));
    session = session_block_array{pos_par,2}(1,:);
    
    for w = 1:length(wave)
        for c = 4:length(condition) % ONLY 31 ONWARDS
            for mode =1:length(modality_data)
            file_names = {}; % create one for each p, c and w, so that it includes all sessions and blocks
            list_order = 1;
            for s = 1:length(session)
                
                % Define blocks within this session
                pos_ses = find(strcmp(session_block_array{pos_par,2}(1,:), session{s}));
                block = session_block_array{pos_par,2}{2,pos_ses}(1,:);
                
                folders = dir([root_dir_bs '/data/' participant{p} '/']);
                results =  contains({folders.name},[wave{w} '_' condition{c} session{s}]) & contains({folders.name},'_temp_baseline');
                infolder = find(results); % So get all blocks you find for that session
                if isempty(infolder)
                    continue
                    % error(['no ' wave{w} '_' condition{c} session{s} '_template blocks for ' participant{p}]);
                end
                for l= 1:length(infolder) % for each block found for that session
                    line = infolder(l);
                    % Check if that block is in the block list, if not, go to next block
                    init_block_name = regexp(folders(line).name, '_b', 'once');
                    fol_nam = contains(block, folders(line).name(init_block_name+1:init_block_name+2));
                    in_fol_nam = find(fol_nam); % index of position in logical
                    if isempty(in_fol_nam); continue; end % the block in the folder name is not present in the list, so go to next
                    % End checking if block exist in list
                    sub_folders = dir([root_dir_bs '/data/' participant{p} '/' folders(line).name '/']);
                    sub_results =  contains({sub_folders.name},[modality_data{mode} '_average_sensor']);
                    sub_infolder = find(sub_results);
                    if isempty(sub_infolder) % if average is not present despite the folder existing
                        continue
                    end
                    file_names{list_order} = [participant{p} '/' folders(line).name '/' sub_folders(sub_infolder).name]; % sub_results (using logicals)
                    list_order = list_order + 1;
                end
            end
            file_names = file_names(~cellfun('isempty', file_names')); % to avoid empty cells
            sFiles = file_names;

            if isempty(sFiles)
            error('for whatever reason, sFiles is empty, my friend...')
            end
            
            disp(' ');      
            disp('-------------------------');  
            disp(['Averaging LLR blocks from all sessions (template epochs) for ' participant{p} '_' modality_data{mode} '_' wave{w} '_' condition{c}]);
            disp(datetime)
            disp(' ');

            % If stated, find and delete previous template averages across
            if delete_previous_file == 1
                % check if there is already a template average for
                % each session in intra folder and delete it
                folders_delete = dir([root_dir_bs '/data/' participant{p} '/@intra/']);
                % data_LLR_average_200727_1504_Template_both_mod_LLR_51
                results_delete =  endsWith({folders_delete.name}, ['Template_' modality_data{mode} '_' wave{w} '_' condition{c} '.mat']); 
                infolder_delete = find(results_delete);
                if ~isempty(infolder_delete) % file exists, therefore delete it
                   delete([root_dir_bs '/data/' participant{p} '/@intra/' folders_delete(infolder_delete).name]);
                end
            end

            % Process: Weighted Average: Everything
            sFiles = bst_process('CallProcess', 'process_average', sFiles, [], ...
                'avgtype',       1, ...  % Everything
                'avg_func',      1, ...  % Arithmetic average:  mean(x)
                'weighted',      1, ...
                'keepevents',    0);


            % Process: Add tag
            sFiles = bst_process('CallProcess', 'process_add_tag', sFiles, [], ...
                'tag',           ['Template_' modality_data{mode} '_' wave{w} '_' condition{c}], ...
                'output',        2);  % Add to file name (1 to add a tag)

            % Process: Set name
            sFiles = bst_process('CallProcess', 'process_set_comment', sFiles, [], ...
                'tag',           ['Template_' modality_data{mode} '_' wave{w} '_' condition{c}], ...
                'isindex',       1);   
            end
        end     
    end
end

% reset modality_data variable to its original
modality_data = {'EEG','MEG','both_mod'};

clearvars('-except', initialVars{:});
disp 'DONE WITH AVERAGING LLR WAVEFORMS ACROSS BLOCKS FROM ALL SESSIONS (TEMPLATE EPOCHS)!!!'
disp(datetime)
toc

end

%% Get TEMPLATE data out of Brainstorm and (if chosen) apply half hanning

if subtraction_approach ==  1

% Chose whether you want the sensor data extracted to be the best of each
% (MEG and EEG amplitude thresholds separated), or the combined survivor
% sweeps (by default is the best of each)
extracted_sensor_data = 1; % 1 = separated amplitude threshold 2 = combined
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sensor_analysis == 1
    modality_data = {'EEG'};
    size_mode_iteration = 1;
elseif sensor_analysis == 2
    modality_data = {'MEG'};
    size_mode_iteration = 1;
elseif sensor_analysis == 3
    modality_data = {'both_mod'};
    size_mode_iteration = 1;
elseif sensor_analysis == 4 % only then you can chosee between getting amplitudes from combined or separated
    if extracted_sensor_data == 1
        size_mode_iteration = [1:2]; % thus covering only EEG and MEG
    elseif extracted_sensor_data == 2
        size_mode_iteration = 3;
    end
end
    
    
tic
disp(' ');      
disp('-------------------------');  
disp('GETTING TEMPLATES OUT OF BRAINSTORM (TO PLOT TEMPLATES)');
disp(datetime)
disp('-------------------------');
disp(' ');

% Averages of individual blocks
for p = 1:length(participant)
    
    % reload subject normally
    prot_subs = bst_get('ProtocolSubjects');
    current_sub = find(strcmp({prot_subs.Subject.Name}, participant{p}));
    db_reload_conditions(current_sub);
    
    % Define sessions within this participant
    pos_par = find(strcmp(session_block_array(:,1), participant{p}(1:4)));
    session = session_block_array{pos_par,2}(1,:);
    
    eval(['no_template_averages_' participant{p} ' = {};'])
    for s = 1:length(session)     
        % Define blocks within this session
        pos_ses = find(strcmp(session_block_array{pos_par,2}(1,:), session{s}));
        block = session_block_array{pos_par,2}{2,pos_ses}(1,:);
        
        for b = 1:length(block)
            for mode = size_mode_iteration
                for w = 1:length(wave)
                    for c = 4:length(condition) % ONLY 31 ONWARDS
                        folders = dir([root_dir_bs '/data/' participant{p} '/']);
                        results =  contains({folders.name},[wave{w} '_' condition{c} session{s} '_' block{b} '_temp_baseline']);
                        infolder = find(results);
                        if isempty(infolder) % If that session/condition does not exist, continue
                            continue
                        end
                        sub_folders = dir([root_dir_bs '/data/' participant{p} '/' folders(infolder).name '/']);
                        sub_results =  contains({sub_folders.name},[modality_data{mode} '_average_sensor']); % Only those averages that are interpolated
                        sub_infolder = find(sub_results);
                        if isempty(sub_infolder)
                            eval(['emp_pos = size(cellfun(@isempty,no_template_averages_' participant{p} '),1) + 1;']) % next empty position in cell array
                            eval(['no_template_averages_' participant{p} '{emp_pos,1} = [wave{w} ''_'' condition{c} session{s} ''_'' block{b}];'])
                            continue
                        end
                        file_name = [root_dir_bs '/data/' participant{p} '/' folders(infolder).name '/' sub_folders(sub_infolder).name];
                        % This way so that we can preserve the structure for the plotting script
                        eval(['load ' file_name])
                        if hann_wind_template == 1
                            % Build half hanning window and apply to template
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            if strcmp(wave{w}, 'MLR') % JUST LEAVE IT HERE FOR NOW TO USE IT IN THE MLR SCRIPT
                                hann_start_sampl = (hann_start + 0.05)*1500;
                                pair = 0;
                            elseif strcmp(wave{w}, 'LLR')
                                hann_start_sampl = ((hann_start + 0.150)*resample_LLR);
                                pair = 1;
                            end
                            size_edge = size(F,2)-hann_start_sampl; 
                            edge = linspace(1,0,size_edge);
                            init_window = ones(1,pair + round(hann_start_sampl));
                            han_window = horzcat(init_window,edge);
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            F = F.*han_window(1:size(F,2)); % because sometimes there is one sample more (like 376 vs 377)
                        end
                        
                        if strcmp(modality_data{mode}, 'EEG')
                            F_EEG = F(315:376,:);
                            clear F
                            F = F_EEG;
                            if ~exist([root_dir '/Averages/' participant{p} '/template/' session{s} '/' block{b}], 'dir')
                               mkdir([root_dir '/Averages/'], [participant{p} '/template/' session{s} '/' block{b}]);
                            end
                            save([root_dir '/Averages/' participant{p} '/template/' session{s} '/' block{b} '/EEG_' wave{w} '_' condition{c} '.mat'], 'F');
                            clear F
                        elseif strcmp(modality_data{mode}, 'MEG')
                            F_MEG = F(1:306,:);
                            clear F
                            if ~exist([root_dir '/Averages/' participant{p} '/template/' session{s} '/' block{b}], 'dir')
                               mkdir([root_dir '/Averages/'], [participant{p} '/template/' session{s} '/' block{b}]);
                            end
                            F = F_MEG;
                            save([root_dir '/Averages/' participant{p} '/template/' session{s} '/' block{b} '/MEG_' wave{w} '_' condition{c} '.mat'], 'F');
                            clear F
                        elseif strcmp(modality_data{mode}, 'both_mod')
                            F_EEG = F(315:376,:);
                            F_MEG = F(1:306,:);
                            clear F
                            F = F_EEG;
                            if ~exist([root_dir '/Averages/' participant{p} '/template/' session{s} '/' block{b}], 'dir')
                               mkdir([root_dir '/Averages/'], [participant{p} '/template/' session{s} '/' block{b}]);
                            end
                            save([root_dir '/Averages/' participant{p} '/template/' session{s} '/' block{b} '/EEG_' wave{w} '_' condition{c} '.mat'], 'F');
                            clear F
                            F = F_MEG;
                            save([root_dir '/Averages/' participant{p} '/template/' session{s} '/' block{b} '/MEG_' wave{w} '_' condition{c} '.mat'], 'F');
                        end
                    end
                end
            end
        end
    end
    save([root_dir '/Events/no_template_averages_' participant{p} '.mat'],['no_template_averages_' participant{p}]);
end

% Averages in individual sessions (to compare quality for substraction using individual sessions)
for p = 1:length(participant)
    
    % reload subject normally
    prot_subs = bst_get('ProtocolSubjects');
    current_sub = find(strcmp({prot_subs.Subject.Name}, participant{p}));
    db_reload_conditions(current_sub);
    
    % Define sessions within this participant
    pos_par = find(strcmp(session_block_array(:,1), participant{p}(1:4)));
    session = session_block_array{pos_par,2}(1,:);
    
    for s = 1:length(session)
        for mode = size_mode_iteration
            for w = 1:length(wave)
                for c = 4:length(condition) % ONLY 31 ONWARDS
                    folders = dir([root_dir_bs '/data/' participant{p} '/@intra/']);
                    results =  contains({folders.name},['Template_' modality_data{mode} '_' wave{w} '_' condition{c} session{s} '.mat']);
                    infolder = find(results);
                    if isempty(infolder) % If that session/condition does not exist, continue
                        continue
                    end
                    file_name = [root_dir_bs '/data/' participant{p} '/@intra/' folders(infolder).name];            
                    % This way so that we can preserve the structure for the plotting script
                    eval(['load ' file_name])
                    if hann_wind_template == 1
                        % Build half hanning window and apply to template
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        if strcmp(wave{w}, 'MLR')
                            hann_start_sampl = (hann_start + 0.05)*1500;
                            pair = 0;
                        elseif strcmp(wave{w}, 'LLR')
                            hann_start_sampl = ((hann_start + 0.150)*resample_LLR);
                            pair = 1;
                        end
                        size_edge = size(F,2)-hann_start_sampl; 
                        edge = linspace(1,0,size_edge);
                        init_window = ones(1,pair + round(hann_start_sampl));
                        han_window = horzcat(init_window,edge);
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        F = F.*han_window(1:size(F,2)); % because sometimes there is one sample more (like 376 vs 377)
                    end
                    
                    if strcmp(modality_data{mode}, 'EEG')
                        F_EEG = F(315:376,:);
                        clear F
                        F = F_EEG;
                        if ~exist([root_dir '/Averages/' participant{p} '/template/' session{s}], 'dir')
                           mkdir([root_dir '/Averages/'], [participant{p} '/template/' session{s}]);
                        end
                        save([root_dir '/Averages/' participant{p} '/template/' session{s} '/EEG_' wave{w} '_' condition{c} '.mat'], 'F');
                        clear F
                    elseif strcmp(modality_data{mode}, 'MEG')
                        F_MEG = F(1:306,:);
                        clear F
                        if ~exist([root_dir '/Averages/' participant{p} '/template/' session{s}], 'dir')
                           mkdir([root_dir '/Averages/'], [participant{p} '/template/' session{s}]);
                        end
                        F = F_MEG;
                        save([root_dir '/Averages/' participant{p} '/template/' session{s} '/MEG_' wave{w} '_' condition{c} '.mat'], 'F');
                        clear F
                    elseif strcmp(modality_data{mode}, 'both_mod')
                        F_EEG = F(315:376,:);
                        F_MEG = F(1:306,:);
                        clear F
                        F = F_EEG;
                        if ~exist([root_dir '/Averages/' participant{p} '/template/' session{s}], 'dir')
                           mkdir([root_dir '/Averages/'], [participant{p} '/template/' session{s}]);
                        end
                        save([root_dir '/Averages/' participant{p} '/template/' session{s} '/EEG_' wave{w} '_' condition{c} '.mat'], 'F');
                        clear F
                        F = F_MEG;
                        save([root_dir '/Averages/' participant{p} '/template/' session{s} '/MEG_' wave{w} '_' condition{c} '.mat'], 'F');
                    end
                end
            end
        end
    end
end

% Template averaged across sessions
for p = 1:length(participant)
    if ~exist([root_dir '/Averages/' participant{p} '/template'], 'dir')
       mkdir([root_dir '/Averages/'], [participant{p} '/template']);
    end
    for mode = size_mode_iteration
        for w = 1:length(wave)
            for c = 4:length(condition)  
            folders = dir([root_dir_bs '/data/' participant{p} '/@intra/']);
            results =  contains({folders.name},['Template_' modality_data{mode} '_' wave{w} '_' condition{c} '.mat']);
            infolder = find(results);
            if isempty(infolder) % If that session/condition does not exist, continue
                continue
            end
            file_name = [root_dir_bs '/data/' participant{p} '/@intra/' folders(infolder).name];            
            % This way so that we can preserve the structure for the plotting script
            eval(['load ' file_name])
            if hann_wind_template == 1
                % Build half hanning window and apply to template
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if strcmp(wave{w}, 'MLR')
                    hann_start_sampl = (hann_start + 0.05)*1500;
                    pair = 0;
                elseif strcmp(wave{w}, 'LLR')
                    hann_start_sampl = ((hann_start + 0.150)*resample_LLR);
                    pair = 1;
                end
                size_edge = size(F,2)-hann_start_sampl; 
                edge = linspace(1,0,size_edge);
                init_window = ones(1,pair + round(hann_start_sampl));
                han_window = horzcat(init_window,edge);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                F = F.*han_window(1:size(F,2)); % because sometimes there is one sample more (like 376 vs 377)
            end
                if strcmp(modality_data{mode}, 'EEG')
                    F_EEG = F(315:376,:);
                    clear F
                    F = F_EEG;
                    save([root_dir '/Averages/' participant{p} '/template/EEG_' wave{w} '_' condition{c} '.mat'], 'F');
                    clear F
                elseif strcmp(modality_data{mode}, 'MEG')
                    F_MEG = F(1:306,:);
                    clear F
                    F = F_MEG;
                    save([root_dir '/Averages/' participant{p} '/template/MEG_' wave{w} '_' condition{c} '.mat'], 'F');
                    clear F
                elseif strcmp(modality_data{mode}, 'both_mod')
                    F_EEG = F(315:376,:);
                    F_MEG = F(1:306,:);
                    clear F
                    F = F_EEG;
                    save([root_dir '/Averages/' participant{p} '/template/EEG_' wave{w} '_' condition{c} '.mat'], 'F');
                    clear F
                    F = F_MEG;
                    save([root_dir '/Averages/' participant{p} '/template/MEG_' wave{w} '_' condition{c} '.mat'], 'F');
                end
            end
        end
    end
end


% reset modality_data variable to its original
modality_data = {'EEG','MEG','both_mod'};

clearvars('-except', initialVars{:});
disp 'DONE WITH GETTING LLR TEMPLATES OUT OF BRAINSTORM (TO PLOT TEMPLATES)!!!'
disp(datetime)
toc
end

%% From all 11, 12 and 13, obtain exact time, latency from previous trial and category of previous trial

if subtraction_approach ==  1
% Preparation for subtraction approach
tic
disp(' ');      
disp('-------------------------');  
disp('OBTAINING EVENT TIMES & PREVIOUS ISI AND CATEGORIES (SUBTRACTION APPROACH)');  
disp(datetime)
disp('-------------------------');     
disp(' ');

for p = 1:length(participant)
    
    % Define sessions within this participant
    pos_par = find(strcmp(session_block_array(:,1), participant{p}(1:4)));
    session = session_block_array{pos_par,2}(1,:);
    
    for s = 1:length(session)
        
        % Define blocks within this session
        pos_ses = find(strcmp(session_block_array{pos_par,2}(1,:), session{s}));
        block = session_block_array{pos_par,2}{2,pos_ses}(1,:);
        
        for b = 1:length(block)
            
            % Define runs within this block
            pos_block = find(strcmp(session_block_array{pos_par,2}{2,pos_ses}(1,:), block{b}));
            runs_in_block = session_block_array{pos_par,2}{2,pos_ses}{2,pos_block}(1,:);
            
            for r = 1:length(runs_in_block)
                
                folders = dir([root_dir '/analysis/' participant{p}(1:4) session{s} '/']);
                results =  contains({folders.name},runs_in_block{r}) & endsWith({folders.name},'vmrk');
                infolder = find(results);
                    if isempty(infolder)
                        continue
                    end
                    if size(infolder,2)> 1 % in case of more than one coincidence, vmrk is going to be the same for all, so pick first
                        infolder = infolder(1);
                    end
                filename = folders(infolder).name;

                fullfilename = [root_dir '/analysis/' participant{p}(1:4) session{s} '/' filename];
                if ~exist(fullfilename, 'file')
                    disp (['filename ' fullfilename ' does not exist'])
                    continue
                end
                fileID = fopen(fullfilename,'r');
                textscan(fileID, '%[^\n\r]', startRow-1, 'WhiteSpace', '', 'ReturnOnError', false, 'EndOfLine', '\r\n');
                dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'ReturnOnError', false);
                fclose(fileID);
                current_run = runs_in_block{r}(10:end-1); % extract the number from the "Click_run9_", even with two numbers (10)
                eval(['eve_' participant{p}(1:4) session{s} '_r' current_run '= [dataArray{1:end-1}];'])
                eval(['events = eve_' participant{p}(1:4) session{s} '_r' current_run ';'])
                latencies_11 = []; count_11 = 1;
                latencies_12 = []; count_12 = 1;
                latencies_13 = []; count_13 = 1;
                for e = 1:length(events)
                    if events(e,1)== '11'
                        if events(e-1,1) == 'boundary'
                            continue
                        elseif events(e-1,1) == ''
                            continue
                        end
                        % Absolute time (first column, in seconds)
                        latencies_11(count_11,1) = str2double(events(e,2))/1500; % since this is vmrk archive, no resampling values are necessary
                        % Latency from previous trial (second column, in seconds)
                        latencies_11(count_11,2) = (str2double(events(e,2)) - str2double(events(e-1,2)))/1500; % since this is vmrk archive, no resampling values are necessary
                        % Category of previous trial (third column)
                        latencies_11(count_11,3) = str2double(events(e-1,1));
                        count_11 = count_11 + 1;
                    elseif events(e,1)== '12'
                        if events(e-1,1) == 'boundary'
                            continue
                        elseif events(e-1,1) == ''
                            continue
                        end
                        % Absolute time (first column, in seconds)
                        latencies_12(count_12,1) = str2double(events(e,2))/1500; % since this is vmrk archive, no resampling values are necessary
                        % Latency from previous trial (second column, in seconds)
                        latencies_12(count_12,2) = (str2double(events(e,2)) - str2double(events(e-1,2)))/1500; % since this is vmrk archive, no resampling values are necessary
                        % Category of previous trial (third column)
                        latencies_12(count_12,3) = str2double(events(e-1,1));
                        count_12 = count_12 + 1;
                    elseif events(e,1)== '13' %#ok<*BDSCA>
                        if events(e-1,1) == 'boundary'
                            continue
                        elseif events(e-1,1) == ''
                            continue
                        end
                        % Absolute time (first column, in seconds)
                        latencies_13(count_13,1) = str2double(events(e,2))/1500; % since this is vmrk archive, no resampling values are necessary
                        % Latency from previous trial (second column, in seconds)
                        latencies_13(count_13,2) = (str2double(events(e,2)) - str2double(events(e-1,2)))/1500; % since this is vmrk archive, no resampling values are necessary
                        % Category of previous trial (third column)
                        latencies_13(count_13,3) = str2double(events(e-1,1));
                        count_13= count_13 + 1;
                    end
                end
                eval(['Log_11_' participant{p} session{s} '_r' current_run '= latencies_11;'])
                eval(['Log_12_' participant{p} session{s} '_r' current_run '= latencies_12;'])
                eval(['Log_13_' participant{p} session{s} '_r' current_run '= latencies_13;'])
                clearvars filename fileID dataArray ans;
            end
        end
    end
    save([root_dir '/Events/Latencies_Fastest_condition_' participant{p} '.mat'],'-regexp','^Log_') 
end

clearvars('-except', initialVars{:});
disp 'DONE WITH OBTAINING EVENT TIMES & PREVIOUS ISI AND CATEGORIES (SUBTRACTION APPROACH)!!!'
disp(datetime)
toc
end

%% SUBTRACT from every sweep vectors obtained from templates (LLR)

if subtraction_approach ==  1
% NEVER PUT A SUBJECT NAME (E.G. 2259_M1_20Hz) THAT CONTAINS THE WORD 'RUN' ON IT

tic
disp(' ');      
disp('-------------------------');  
disp('LLR SUBTRACTION OF ARTEFACT FROM SHORTEST ISI (EFFECTS FROM PREVIOUS TRIAL)');  
disp(datetime)
disp('-------------------------');     
disp(' ');

for p = 1:length(participant)
    eval(['load ''' root_dir '/Events/Latencies_Fastest_condition_' participant{p} '.mat'''])
    
    % Define sessions within this participant
    pos_par = find(strcmp(session_block_array(:,1), participant{p}(1:4)));
    session = session_block_array{pos_par,2}(1,:);
    
    for s = 1:length(session)
        
        % Define blocks within this session
        pos_ses = find(strcmp(session_block_array{pos_par,2}(1,:), session{s}));
        block = session_block_array{pos_par,2}{2,pos_ses}(1,:);
        
        for b = 1:length(block)
            
            % Define runs within this block
            pos_block = find(strcmp(session_block_array{pos_par,2}{2,pos_ses}(1,:), block{b}));
            runs_in_block = session_block_array{pos_par,2}{2,pos_ses}{2,pos_block}(1,:);
            
            for w = 1:length(wave) % WILL BE LLR ONLY IN THIS CASE
                for si=1:length(shortest_ISI)
                eval(['ncs_' participant{p} '_' shortest_ISI{si} session{s} '_' block{b} '_' wave{w} '=0;'])    
                folders = dir([root_dir_bs '/data/' participant{p} '/']);
                results = contains({folders.name},[wave{w} '_' shortest_ISI{si} session{s} '_' block{b} '_normal']);
                infolder = find(results);
                    if isempty(infolder)
                        continue
                    end
                    dir_trial = dir([root_dir_bs '/data/' participant{p} '/' folders(infolder).name]);
                    trial_cov = contains({dir_trial.name},'_trial');
                    trial_count = find(trial_cov);
                    
                    disp(' ');      
                    disp('-------------------------');  
                    disp(['LLR Subtraction from shortest ISI for ' participant{p} '_' session{s} '_' block{b} '_' wave{w} '_' shortest_ISI{si}]);
                    disp(datetime)
                    disp(' ');  
                    
                    for j= 1:length(trial_count) % for every trial
                        position = trial_count(j);
                        trial_name = [root_dir_bs '/data/' participant{p} '/' folders(infolder).name '/' dir_trial(position).name];
                        eval(['load ' trial_name]) % trial loaded 
                        % contains(ChannelFlag ColormapType Comment DataType Device DisplayUnits Events F History Leff Std Time nAvg)
                        % Retrieve value and run number of each trial
                        
                        [x,y]= find(contains(History,'import_time'));
                        if isempty(x); continue; end % Just in case
                        
                        brack = regexp(History{x,y+1},'[');com = regexp(History{x,y+1},',');      
                        % latency absolute value (for match with log data)
                        trial_time = str2double(History{x,y+1}(brack +1: com -1));
                        if strcmp(wave{w}, 'MLR')
                            trial_time = trial_time + 0.05; % adjust epoch time to stimuly time
                        elseif strcmp(wave{w}, 'LLR') % WILL ALWAYS BE LLR IN THIS CASE CASE
                            trial_time = trial_time + 0.150;
                        end
                        % run number
                        cell_import = find(contains(History,'Import from:'));
                        init_run_name = regexp(History{cell_import}, 'run', 'once');
                        end_run_name = regexp(History{cell_import}, '_mc_tsss', 'once');
                        trial_run = History{cell_import}(init_run_name+3:end_run_name-1); % gives the two digit number of the run
                        % Check if any Log file is empty before looking into it
                        eval(['empty_logs = isempty(Log_' shortest_ISI{si} '_' participant{p} session{s} '_r' trial_run ');'])
                        if empty_logs == 1 % if the log is empty
                            continue
                        end     
                        % find position in Log
                        eval(['ident_pos = find(abs(trial_time  - Log_' shortest_ISI{si} '_' participant{p} session{s} '_r' trial_run '(:,1))<10^-2);']) % Was three before, maybe that's why we lost some trials
                        if isempty(ident_pos)
                            eval(['ncs_' participant{p} '_' shortest_ISI{si} session{s} '_' block{b} '_' wave{w} ' = ncs_' participant{p} '_' shortest_ISI{si} session{s} '_' block{b} '_' wave{w} ' +1;'])
                            continue
                        end
                        % Previous category info
                        eval(['prev_cat = Log_' shortest_ISI{si} '_' participant{p} session{s} '_r' trial_run '(ident_pos,3);'])
                        % avoid rare random categories that I don't expect (like '8' or '' or 'boundary'), I rather do this
                        if (prev_cat ~= 31) && (prev_cat ~= 32) && (prev_cat ~= 33) &&...
                                (prev_cat ~= 51) && (prev_cat ~= 52) && (prev_cat ~= 53) &&...
                            (prev_cat ~= 71) && (prev_cat ~= 72) && (prev_cat ~= 73) &&...
                            (prev_cat ~= 91) && (prev_cat ~= 92) && (prev_cat ~= 93)
                            continue
                        end
                        eval(['prev_lat = Log_' shortest_ISI{si} '_' participant{p} session{s} '_r' trial_run '(ident_pos,2);'])
                        if (0.249 > prev_lat) || (prev_lat > 0.501)
                            % To fix rare values like Log 2494A_r1_row 4
                            continue
                        end
                        % With this info, retrieve values from templates
                        if template == 1 % block
                            templ_directory = [root_dir '/Averages/' participant{p} '/template/' session{s} '/' block{b}];
                        elseif template == 2 % session  
                            templ_directory = [root_dir '/Averages/' participant{p} '/template/' session{s}];
                        elseif template == 3 % total
                            templ_directory = [root_dir '/Averages/' participant{p} '/template'];
                        end
                        % save trial amplitudes before overwriting F variable
                        trial_amplitudes = F; % always gonna be 392 channels
                        length_trial = size(trial_amplitudes,2);
                        % define init points in vectors to subtract
                        init_point = prev_lat*resample_LLR; % IN THE CASE OF LLR ONLY!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        % IMPORTANT!! was (0.05*1500)+(prev_lat*1500)
                        % before, but I noticed then we do not account
                        % for the baseline effects at all, only for the
                        % portion after the stimulus. So it should be 
                        % (0.05*1500)+(prev_lat*1500) - baseline again
                        % = (prev_lat*1500) alone. Therefore there is
                        % only one init point regardless of the wave
                        % Load templates and subtract with different channels
                        
                        % EEG channels
                        if ~exist([templ_directory '/EEG_' wave{w} '_' num2str(prev_cat) '.mat'], 'file')
                            % It could be that we don't have the template
                            % for the block because no sweeps survided to
                            % create the template average for the block, or
                            % even for a whole session. So, if the template
                            % chosen was the block and that one does not
                            % exist, we may try the session. If only one
                            % block survived and we don't even have the
                            % average of the session, we may use the
                            % subject average. And if the subject average
                            % is empty, we are fucked, so better continue
                            % with next iteration
                            if template == 1 % I chose block and it's empty, let's try with session
                                % Change template directory to session
                                templ_directory = [root_dir '/Averages/' participant{p} '/template/' session{s}];
                                % If session template does not exist either
                                if ~exist([templ_directory '/EEG_' wave{w} '_' num2str(prev_cat) '.mat'], 'file')
                                   % Change template directory to subject
                                    templ_directory = [root_dir '/Averages/' participant{p} '/template'];
                                   % If subject template does not exist, the only option is to continue
                                   if ~exist([templ_directory '/EEG_' wave{w} '_' num2str(prev_cat) '.mat'], 'file')
                                       continue
                                   end
                                end
                            elseif template == 2 % I chose session and it's empty 
                                % change template directory to subject
                                templ_directory = [root_dir '/Averages/' participant{p} '/template'];
                                % If subject template does not exist, the only option is to continue
                                if ~exist([templ_directory '/EEG_' wave{w} '_' num2str(prev_cat) '.mat'], 'file')
                                    continue
                                end
                            elseif template == 3 % I chose subject/total and it's empty
                                % If that template does not exist, continue
                                % with next trial (for which the previous
                                % category may not be the one for which we
                                % don't have the template)
                                % (If EEG template does not exist, MEG one won't either) 
                                continue
                            end
                        end
                        eval(['load ''' templ_directory '/EEG_' wave{w} '_' num2str(prev_cat) '.mat'''])
                        length_templ_vector = size(F(1:62,fix(init_point):end),2);  
                        if length_templ_vector > length_trial
                        subtract_vector = F(1:62,fix(init_point):fix(init_point)+fix(length_trial-1));
                            if hann_wind_subtract == 1 % Apply hanning window to vector to subtract
                                % Build hanning window based on length of vector
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                size_edges = round(0.05*(size(subtract_vector,2))); % 5% of vector's legth
                                edge_1 = linspace(0,1,size_edges); edge_2 = linspace(1,0,size_edges);
                                middle_window = ones(1,size(subtract_vector,2) - size_edges*2);
                                han_window = horzcat(edge_1,middle_window,edge_2);
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                subtract_vector = subtract_vector.*han_window;
                            end  
                            trial_amplitudes(315:376,:) = trial_amplitudes(315:376,:) - subtract_vector;
                        else
                            subtract_vector = F(1:62,fix(init_point):end);
                            if hann_wind_subtract == 1 % Apply hanning window to vector to subtract
                                % Build hanning window based on length of vector
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                size_edges = round(0.05*(size(subtract_vector,2))); % 5% of vector's legth
                                edge_1 = linspace(0,1,size_edges); edge_2 = linspace(1,0,size_edges);
                                middle_window = ones(1,size(subtract_vector,2) - size_edges*2);
                                han_window = horzcat(edge_1,middle_window,edge_2);
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                subtract_vector = subtract_vector.*han_window;
                            end  
                            trial_amplitudes(315:376,1:fix(length_templ_vector)) = trial_amplitudes(315:376,1:fix(length_templ_vector)) - subtract_vector;
                        end
                        % MEG channels
                        if ~exist([templ_directory '/MEG_' wave{w} '_' num2str(prev_cat) '.mat'], 'file')
                            % It could be that we don't have the template
                            % for the block because no sweeps survided to
                            % create the template average for the block, or
                            % even for a whole session. So, if the template
                            % chosen was the block and that one does not
                            % exist, we may try the session. If only one
                            % block survived and we don't even have the
                            % average of the session, we may use the
                            % subject average. And if the subject average
                            % is empty, we are fucked, so better continue
                            % with next iteration
                            if template == 1 % I chose block and it's empty, let's try with session
                                % Change template directory to session
                                templ_directory = [root_dir '/Averages/' participant{p} '/template/' session{s}];
                                % If session template does not exist either
                                if ~exist([templ_directory '/MEG_' wave{w} '_' num2str(prev_cat) '.mat'], 'file')
                                   % Change template directory to subject
                                    templ_directory = [root_dir '/Averages/' participant{p} '/template'];
                                   % If subject template does not exist, the only option is to continue
                                   if ~exist([templ_directory '/MEG_' wave{w} '_' num2str(prev_cat) '.mat'], 'file')
                                       continue
                                   end
                                end
                            elseif template == 2 % I chose session and it's empty 
                                % change template directory to subject
                                templ_directory = [root_dir '/Averages/' participant{p} '/template'];
                                % If subject template does not exist, the only option is to continue
                                if ~exist([templ_directory '/MEG_' wave{w} '_' num2str(prev_cat) '.mat'], 'file')
                                    continue
                                end
                            elseif template == 3 % I chose subject/total and it's empty
                                % If that template does not exist, continue
                                % with next trial (for which the previous
                                % category may not be the one for which we
                                % don't have the template)
                                % (If EEG template does not exist, MEG one won't either) 
                                continue
                            end
                        end
                        eval(['load ''' templ_directory '/MEG_' wave{w} '_' num2str(prev_cat) '.mat'''])
                        length_templ_vector = size(F(1:306,fix(init_point):end),2);
                        if length_templ_vector > length_trial
                            subtract_vector = F(1:306,fix(init_point):fix(init_point)+fix(length_trial-1));
                            if hann_wind_subtract == 1 % Apply hanning window to vector to subtract
                                % Build hanning window based on length of vector
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                size_edges = round(0.05*(size(subtract_vector,2))); % 5% of vector's legth
                                edge_1 = linspace(0,1,size_edges); edge_2 = linspace(1,0,size_edges);
                                middle_window = ones(1,size(subtract_vector,2) - size_edges*2);
                                han_window = horzcat(edge_1,middle_window,edge_2);
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                subtract_vector = subtract_vector.*han_window;
                            end  
                            trial_amplitudes(1:306,:) = trial_amplitudes(1:306,:) - subtract_vector;
                        else
                            subtract_vector = F(1:306,fix(init_point):end);
                            if hann_wind_subtract == 1 % Apply hanning window to vector to subtract
                                % Build hanning window based on length of vector
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                size_edges = round(0.05*(size(subtract_vector,2))); % 5% of vector's legth
                                edge_1 = linspace(0,1,size_edges); edge_2 = linspace(1,0,size_edges);
                                middle_window = ones(1,size(subtract_vector,2) - size_edges*2);
                                han_window = horzcat(edge_1,middle_window,edge_2);
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                subtract_vector = subtract_vector.*han_window;
                            end  
                            trial_amplitudes(1:306,1:fix(length_templ_vector)) = trial_amplitudes(1:306,1:fix(length_templ_vector)) - subtract_vector;
                        end
                        clear F
                        F = trial_amplitudes; % new trial amplitudes to be saved
                        % rebuild trigger file with new F values
                        save(trial_name, 'ChannelFlag', 'ColormapType', 'Comment', 'DataType', 'Device', 'DisplayUnits', 'Events', 'F', 'History', 'Leff', 'Std', 'Time', 'nAvg')
                        clear ChannelFlag ColormapType Comment DataType Device DisplayUnits Events F History Leff Std Time nAvg
                    end
                end
            end
        end
    end
    save([root_dir '/Events/Non_coincident_sweeps_' participant{p} '.mat'],'-regexp','^ncs_');
end

clearvars('-except', initialVars{:});
disp 'DONE WITH LLR SUBTRACTION OF ARTEFACT FROM SHORTEST ISI (EFFECTS FROM PREVIOUS TRIAL)!!!'
disp(datetime)
toc
end

%% Epoching NORMALLY (PART 2)

tic
disp(' ');      
disp('-------------------------');  
disp('EPOCHING NORMALLY (PART 2: after subtraction)');  
disp(datetime)
disp('-------------------------');     
disp(' ');

% Import filtered data this time (matlab_low folders only for LLR)
for p = 1:length(participant)
    
    disp(' ');      
    disp('-------------------------');
    disp(['loading participant ' participant{p}]);
    disp(datetime)
    disp(' '); 
    
    prot_subs = bst_get('ProtocolSubjects');
    current_sub = find(strcmp({prot_subs.Subject.Name}, participant{p}));
    db_reload_conditions(current_sub);
    
    % Define sessions within this participant
    pos_par = find(strcmp(session_block_array(:,1), participant{p}(1:4)));
    session = session_block_array{pos_par,2}(1,:);
    
    for s = 1:length(session)
        
        % Define blocks within this session
        pos_ses = find(strcmp(session_block_array{pos_par,2}(1,:), session{s}));
        block = session_block_array{pos_par,2}{2,pos_ses}(1,:);
        
        for b = 1:length(block)
            for c =1:length(condition)
                for w = 1:length(wave)
                    
                    folders = dir([root_dir_bs '/data/' participant{p} '/']);
                    results = contains({folders.name},[wave{w} '_' condition{c} session{s} '_' block{b} '_normal']);
                    infolder = find(results);
                    
                    if isempty(infolder)
                        continue
                    end
                    if size(infolder,2)> 1 % in case of more than one coincidence, error
                        error('More than one coincidence');
                    end
                        
                    file_names = {};
                    dir_sweeps = dir([root_dir_bs '/data/' participant{p} '/' folders(infolder).name]);
                    sweeps_norm = contains({dir_sweeps.name},'_trial');
                    sweeps_list = find(sweeps_norm);
                    for j= 1:length(sweeps_list)
                        position = sweeps_list(j);
                        file_names{j} = [participant{p} '/' folders(infolder).name '/' dir_sweeps(position).name];
                    end
                sFiles = file_names; % all sweeps from that block and condition contained here now
                
                if isempty(sFiles)
                    % It can be that sometimes certain sweeps types (2235
                    % 52 C b4 do not exist while others do, because for
                    % whatever reason there were no 52 in a block). Since
                    % it has happened, better to continue here
                	continue
                end
                
                if subtraction_approach == 1 % so if subtraction has ocurred, baseline  correct again    
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Baseline correction normal epochs %%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                disp(' ');      
                disp('-------------------------');  
                disp(['Baseline correction for normal epochs ' participant{p} '_' wave{w} '_' condition{c} session{s} '_' block{b}]);
                disp(datetime)
                disp(' ');

                % Process: DC offset correction for normal epochs
                sFiles = bst_process('CallProcess', 'process_baseline_norm', sFiles, [], ...
                    'baseline',    epoch_baseline{w}, ... % because it´s downsampled
                    'sensortypes', '', ...
                    'method',      'bl', ...  % DC offset correction:    x_std = x - &mu;
                    'overwrite',   1);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Amplitude threshold normal epochs %%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if (sensor_analysis == 1) || (sensor_analysis == 4) % EEG only
                disp(' ');      
                disp('-------------------------');  
                disp(['Cleaning normal epochs for ' participant{p} '_' wave{w} '_' condition{c} session{s} '_' block{b} '(EEG)']);
                disp(datetime)
                disp(' '); 
                
                % Process: Detect bad trials: Absolute threshold EEG
                sFiles_EEG = bst_process('CallProcess', 'process_CNRL_detectbad', sFiles, [], ...
                    'timewindow', [], ...
                    'meggrad',    [0, 0], ...
                    'megmag',     [0, 0], ...
                    'eeg',        reject_EEG_absolute, ...
                    'ieeg',       [0, 0], ...
                    'eog',        [0, 0], ...
                    'ecg',        [0, 0], ...
                    'rejectmode', 2);  % Reject the entire trial
                
                % save and reset bad trials
                trials_file = [root_dir_bs '/data/' participant{p} '/' folders(infolder).name '/brainstormstudy.mat'];
                eval(['load ' trials_file]); % load them
                % save current ones in events folder
                if ~exist([root_dir '/Events/BadTrials/' ProtocolName(7:9) '/EEG/' participant{p} '/' folders(infolder).name], 'dir')
                    mkdir([root_dir '/Events/BadTrials/'], [ProtocolName(7:9) '/EEG/' participant{p} '/' folders(infolder).name]);
                end
                save([root_dir '/Events/BadTrials/' ProtocolName(7:9) '/EEG/' participant{p} '/' folders(infolder).name '/BadTrials.mat'],'BadTrials'); % save bad trials
                % reset them and save back to original file
                BadTrials = cell([], 1); % exact structure it has when empty
                save(trials_file,'BadTrials', 'DateOfStudy', 'Name');
                end
                
                % BadTrials, DateOfStudy, Name
                
                if (sensor_analysis == 2) || (sensor_analysis == 4) % MEG only
                    
                disp(' ');      
                disp('-------------------------');  
                disp(['Cleaning normal epochs for ' participant{p} '_' wave{w} '_' condition{c} session{s} '_' block{b} '(MEG)']);
                disp(datetime)
                disp(' '); 
                
                % Process: Detect bad trials: Peak-to-peak  MEG ONLY
                sFiles_MEG = bst_process('CallProcess', 'process_CNRL_detectbad',sFiles, [], ...
                    'timewindow', [], ...
                    'meggrad',    reject_MEG_GRAD_absolute, ...
                    'megmag',     reject_MEG_MAG_absolute, ...
                    'eeg',        [0, 0], ...
                    'ieeg',       [0, 0], ...
                    'eog',        [0, 0], ...
                    'ecg',        [0, 0], ...
                    'rejectmode', 2);  % Reject the entire trial
                
                % save and reset bad trials
                trials_file = [root_dir_bs '/data/' participant{p} '/' folders(infolder).name '/brainstormstudy.mat'];
                eval(['load ' trials_file]); % load them
                % save current ones in events folder
                if ~exist([root_dir '/Events/BadTrials/' ProtocolName(7:9) '/MEG/' participant{p} '/' folders(infolder).name], 'dir')
                    mkdir([root_dir '/Events/BadTrials/'], [ProtocolName(7:9) '/MEG/' participant{p} '/' folders(infolder).name]);
                end
                save([root_dir '/Events/BadTrials/' ProtocolName(7:9) '/MEG/' participant{p} '/' folders(infolder).name '/BadTrials.mat'],'BadTrials'); % save bad trials
                % reset them and save back to original file
                BadTrials = cell([], 1); % exact structure it has when empty
                save(trials_file,'BadTrials', 'DateOfStudy', 'Name');
                end
                
                if (sensor_analysis == 3) || (sensor_analysis == 4) % Combined EEG and MEG
                
                disp(' ');      
                disp('-------------------------');  
                disp(['Cleaning normal epochs for ' participant{p} '_' wave{w} '_' condition{c} session{s} '_' block{b} '(both_mod)']);
                disp(datetime)
                disp(' '); 
                
                % Process: Detect bad trials: Peak-to-peak  BOTH EEG AND MEG
                sFiles_both_mod = bst_process('CallProcess', 'process_CNRL_detectbad',sFiles, [], ...
                    'timewindow', [], ...
                    'meggrad',    reject_MEG_GRAD_absolute, ...
                    'megmag',     reject_MEG_MAG_absolute, ...
                    'eeg',        reject_EEG_absolute, ...
                    'ieeg',       [0, 0], ...
                    'eog',        [0, 0], ...
                    'ecg',        [0, 0], ...
                    'rejectmode', 2);  % Reject the entire trial
                
                % save and reset bad trials
                trials_file = [root_dir_bs '/data/' participant{p} '/' folders(infolder).name '/brainstormstudy.mat'];
                eval(['load ' trials_file]); % load them
                % save current ones in events folder
                if ~exist([root_dir '/Events/BadTrials/' ProtocolName(7:9) '/both_mod/' participant{p} '/' folders(infolder).name], 'dir')
                    mkdir([root_dir '/Events/BadTrials/'], [ProtocolName(7:9) '/both_mod/' participant{p} '/' folders(infolder).name]);
                end
                save([root_dir '/Events/BadTrials/' ProtocolName(7:9) '/both_mod/' participant{p} '/' folders(infolder).name '/BadTrials.mat'],'BadTrials'); % save bad trials
                % reset them and save back to original file
                BadTrials = cell([], 1); % exact structure it has when empty
                save(trials_file,'BadTrials', 'DateOfStudy', 'Name');
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Average of sweeps normal epochs %%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if (sensor_analysis == 1) || (sensor_analysis == 4) % EEG only
                    
                % SENSOR AVERAGE EEG (WITH INTERPOLATION)   
                disp(' ');      
                disp('-------------------------');  
                disp(['Averaging normal epochs for ' participant{p} '_' wave{w} '_' condition{c} session{s} '_' block{b} '(EEG sensor)']);
                disp(datetime)
                disp(' '); 
                
                sFiles_EEG_averaged_sensor = bst_process('CallProcess', 'process_average', sFiles_EEG, [], ...
                    'avgtype',         5, ...  % By trial group (folder average)
                    'avg_func',        1, ...  % Arithmetic average:  mean(x)
                    'weighted',        0, ...
                    'keepevents', 0);
                
                % Interpolate channels
                sFiles_EEG_averaged_sensor = bst_process('CallProcess', 'process_eeg_interpbad', sFiles_EEG_averaged_sensor, [], ...
                'maxdist',     5, ...
                'sensortypes', 'EEG', ...
                'overwrite',   1);

                % Process: Add tag
                sFiles_EEG_averaged_sensor = bst_process('CallProcess', 'process_add_tag', sFiles_EEG_averaged_sensor, [], ...
                    'tag',           'EEG_average_sensor', ...
                    'output',        2);  % Add to file name (1 to add a tag)

                % Process: Set name
                sFiles_EEG_averaged_sensor = bst_process('CallProcess', 'process_set_comment', sFiles_EEG_averaged_sensor, [], ...
                    'tag',           'EEG_average_sensor', ...
                    'isindex',       1);
                

                % SOURCE AVERAGE EEG (NO INTERPOLATION)  
                disp(' ');      
                disp('-------------------------');  
                disp(['Averaging normal epochs for ' participant{p} '_' wave{w} '_' condition{c} session{s} '_' block{b} '(EEG source)']);
                disp(datetime)
                disp(' ');
                
                sFiles_EEG_averaged_source = bst_process('CallProcess', 'process_average', sFiles_EEG, [], ...
                    'avgtype',         5, ...  % By trial group (folder average)
                    'avg_func',        1, ...  % Arithmetic average:  mean(x)
                    'weighted',        0, ...
                    'keepevents', 0);
                
                % Process: Add tag
                sFiles_EEG_averaged_source = bst_process('CallProcess', 'process_add_tag', sFiles_EEG_averaged_source, [], ...
                    'tag',           'EEG_average_source', ...
                    'output',        2);  % Add to file name (1 to add a tag)

                % Process: Set name
                sFiles_EEG_averaged_source = bst_process('CallProcess', 'process_set_comment', sFiles_EEG_averaged_source, [], ...
                    'tag',           'EEG_average_source', ...
                    'isindex',       1);
                end
                
                if (sensor_analysis == 2) || (sensor_analysis == 4) % MEG only
                    
                % SOURCE AVERAGE MEG (NO INTERPOLATION)                  
                disp(' ');      
                disp('-------------------------');  
                disp(['Averaging normal epochs for ' participant{p} '_' wave{w} '_' condition{c} session{s} '_' block{b} '(MEG sensor)']);
                disp(datetime)
                disp(' '); 
                
                sFiles_MEG_averaged_sensor = bst_process('CallProcess', 'process_average', sFiles_MEG, [], ...
                    'avgtype',         5, ...  % By trial group (folder average)
                    'avg_func',        1, ...  % Arithmetic average:  mean(x)
                    'weighted',        0, ...
                    'keepevents', 0);
                
                % FIND ONE THAT WORKS FOR MEG
%                 % Interpolate channels
%                 sFiles_MEG_averaged_sensor = bst_process('CallProcess', 'process_eeg_interpbad', sFiles_MEG_averaged_sensor, [], ...
%                 'maxdist',     5, ...
%                 'sensortypes', 'MEG', ...
%                 'overwrite',   1);

                % Process: Add tag
                sFiles_MEG_averaged_sensor = bst_process('CallProcess', 'process_add_tag', sFiles_MEG_averaged_sensor, [], ...
                    'tag',           'MEG_average_sensor', ...
                    'output',        2);  % Add to file name (1 to add a tag)

                % Process: Set name
                sFiles_MEG_averaged_sensor = bst_process('CallProcess', 'process_set_comment', sFiles_MEG_averaged_sensor, [], ...
                    'tag',           'MEG_average_sensor', ...
                    'isindex',       1);
                

                % SOURCE AVERAGE MEG (NO INTERPOLATION)  
                disp(' ');      
                disp('-------------------------');  
                disp(['Averaging normal epochs for ' participant{p} '_' wave{w} '_' condition{c} session{s} '_' block{b} '(MEG source)']);
                disp(datetime)
                disp(' ');
                
                sFiles_MEG_averaged_source = bst_process('CallProcess', 'process_average', sFiles_MEG, [], ...
                    'avgtype',         5, ...  % By trial group (folder average)
                    'avg_func',        1, ...  % Arithmetic average:  mean(x)
                    'weighted',        0, ...
                    'keepevents', 0);
                
                % Process: Add tag
                sFiles_MEG_averaged_source = bst_process('CallProcess', 'process_add_tag', sFiles_MEG_averaged_source, [], ...
                    'tag',           'MEG_average_source', ...
                    'output',        2);  % Add to file name (1 to add a tag)

                % Process: Set name
                sFiles_MEG_averaged_source = bst_process('CallProcess', 'process_set_comment', sFiles_MEG_averaged_source, [], ...
                    'tag',           'MEG_average_source', ...
                    'isindex',       1);
                end
                
                if (sensor_analysis == 3) || (sensor_analysis == 4) % Combined EEG and MEG
                    
                % SOURCE AVERAGE BOTH MOD (NO INTERPOLATION)                  
                disp(' ');      
                disp('-------------------------');  
                disp(['Averaging normal epochs for ' participant{p} '_' wave{w} '_' condition{c} session{s} '_' block{b} '(both mod sensor)']);
                disp(datetime)
                disp(' '); 
                
                sFiles_both_mod_averaged_sensor = bst_process('CallProcess', 'process_average', sFiles_both_mod, [], ...
                    'avgtype',         5, ...  % By trial group (folder average)
                    'avg_func',        1, ...  % Arithmetic average:  mean(x)
                    'weighted',        0, ...
                    'keepevents', 0);
                
                % FIND ONE THAT WORKS FOR MEG and EEG combined
%                 % Interpolate channels
%                 sFiles_MEG_averaged_sensor = bst_process('CallProcess', 'process_eeg_interpbad', sFiles_MEG_averaged_sensor, [], ...
%                 'maxdist',     5, ...
%                 'sensortypes', 'MEG', ...
%                 'overwrite',   1);

                % Process: Add tag
                sFiles_both_mod_averaged_sensor = bst_process('CallProcess', 'process_add_tag', sFiles_both_mod_averaged_sensor, [], ...
                    'tag',           'both_mod_average_sensor', ...
                    'output',        2);  % Add to file name (1 to add a tag)

                % Process: Set name
                sFiles_both_mod_averaged_sensor = bst_process('CallProcess', 'process_set_comment', sFiles_both_mod_averaged_sensor, [], ...
                    'tag',           'both_mod_average_sensor', ...
                    'isindex',       1);
                

                % SOURCE AVERAGE BOTH MOD (NO INTERPOLATION)  
                disp(' ');      
                disp('-------------------------');  
                disp(['Averaging normal epochs for ' participant{p} '_' wave{w} '_' condition{c} session{s} '_' block{b} '(both mod source)']);
                disp(datetime)
                disp(' ');
                
                sFiles_both_mod_averaged_source = bst_process('CallProcess', 'process_average', sFiles_both_mod, [], ...
                    'avgtype',         5, ...  % By trial group (folder average)
                    'avg_func',        1, ...  % Arithmetic average:  mean(x)
                    'weighted',        0, ...
                    'keepevents', 0);
                
                % Process: Add tag
                sFiles_both_mod_averaged_source = bst_process('CallProcess', 'process_add_tag', sFiles_both_mod_averaged_source, [], ...
                    'tag',           'both_mod_average_source', ...
                    'output',        2);  % Add to file name (1 to add a tag)

                % Process: Set name
                sFiles_both_mod_averaged_source = bst_process('CallProcess', 'process_set_comment', sFiles_both_mod_averaged_source, [], ...
                    'tag',           'both_mod_average_source', ...
                    'isindex',       1);
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End Average of sweeps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                end
            end
        end
    end
end
            
clearvars('-except', initialVars{:});
disp 'DONE WITH NORMAL EPOCHING (PART 2: bas corr, clean, average)!!!'
disp(datetime)
toc

%% Count trials for NORMAL epochs FOR LLR
% May consider one for templates, although info is already there for that
% LLR and MLR ones will be in separate files now
% If you run first EEG and later MEG, for instance, files will overwrite.
% You run this section of the scipt once all MEG/EEG files are present so
% that you have one unique count with all three

tic
disp(' ');      
disp('-------------------------');  
disp('COUNTING TRIALS FOR NORMAL EPOCHS LLR'); 
disp(datetime)
disp('-------------------------');     
disp(' ');

if sensor_analysis == 1
    colNames = {'PRE_LLR','EEG_Post'};
    modality_data = {'EEG'};
elseif sensor_analysis == 2
    colNames = {'PRE_LLR','MEG_Post'};
    modality_data = {'MEG'};
elseif sensor_analysis == 3
    colNames = {'PRE_LLR','Combined_Post'};
    modality_data = {'both_mod'};
elseif sensor_analysis == 4
    colNames = {'PRE_LLR','EEG_Post','MEG_Post','Combined_Post'};
end

rowNames_summary = horzcat(condition, Exp_cond);
for p = 1:length(participant)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Normal_trials_detailed = [];
    Normal_trials_summary = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for w = 1:length(modality_data) % Just because of inheritance from previous script, we keep 'w' instead of 'mode'
        iter_modality = 0; % first time iteration for position purposes
        current_size_matrix = []; % just for size purposes
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        sum_orig = zeros(15,1);
        sum_final = zeros(15,1);
        rowNames_detailed = {}; % create a new for every subject (bec of diff num of sessions) 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Define sessions within this participant
        pos_par = find(strcmp(session_block_array(:,1), participant{p}(1:4)));
        session = session_block_array{pos_par,2}(1,:);
        
        for s = 1:length(session) 
            
            % Define blocks within this session
            pos_ses = find(strcmp(session_block_array{pos_par,2}(1,:), session{s}));
            block = session_block_array{pos_par,2}{2,pos_ses}(1,:);
            
            for b = 1:length(block)
                for c =1:length(condition)
                    folders = dir([root_dir_bs '/data/' participant{p} '/']);
                    results = contains({folders.name},['LLR_' condition{c} session{s} '_' block{b} '_normal']);
                    infolder = find(results);
                    if isempty(infolder)
                        % error(['No LLR_' condition{c} session{s} '_' block{b} '_normal for participant ' participant{p}])
                        continue
                    end  
                    
                    disp(' ');      
                    disp('-------------------------');
                    disp(['Counting LLR normal sweeps for ' participant{p} '_LLR_' modality_data{w} '_' session{s} '_' block{b} '_' condition{c}]);
                    disp(' '); 
                    
                    dir_sweeps = dir([root_dir_bs '/data/' participant{p} '/' folders(infolder).name]);
                    sweeps_cov = contains({dir_sweeps.name},'_trial');
                    sweeps_count = find(sweeps_cov);
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % COUNT REJECTED TRIALS
                    try 
                        load([root_dir '/Events/BadTrials/' ProtocolName(7:9) '/' modality_data{w} '/' participant{p} '/' folders(infolder).name '/BadTrials.mat']); 
                        if ~isempty(BadTrials)
                            count_rej = length(BadTrials);
                        else
                            count_rej = 0;
                        end
                    catch
                        % If no BadTrials variable is ready to load, countas 0
                        count_rej = 0;
                    end
                    % COUNT ORIGINAL TRIALS
                    count_orig = length(sweeps_count);
                    % COUNT NON REJECTED TRIALS
                    count_final = length(sweeps_count) - count_rej;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
                    
                    if iter_modality == 0 % first iteration in this modality
                        pos = 1;
                        iter_modality = 1; % this won't change until we go to next modality
                        % length(Normal_trials_detailed) == 0  % how it was
                        % before %#ok<ISMT> % DO NOT CHANGE BY ISEMPTY
                    else % determine last position and put one after
                        % size_matrix = size(Normal_trials_detailed,1);
                        % pos = size_matrix +1; 
                        pos = size(current_size_matrix,1) + 1;
                    end
                    % pos = ((s-1)*45) + ((b-1)*15) + c; % index position in table
                    % assuming each block has 15 conditions, and every session has 3 blocks
                    
                    rowNames_detailed{1,pos} = [session{s} '_' block{b} '_' condition{c}];
                    Normal_trials_detailed(pos,1) = count_orig;
                    Normal_trials_detailed(pos,w+1) = count_final;
                    current_size_matrix(pos,1) = 1; % just to keep track of size
                    sum_orig(c,1) = sum_orig(c,1) + count_orig;
                    sum_final(c,1) = sum_final(c,1) + count_final;
                end
            end
        end
        %%%%%%%%%%%%%%%% add averages across sessions and blocks %%%%%%%%%%%%%%%%%%%%%%%%%%
        for ct = 1:15
            Normal_trials_summary(ct,1) = sum_orig(ct,1);
            Normal_trials_summary(ct,w+1) = sum_final(ct,1);
        end
        % Quietest
        Normal_trials_summary(16,1) = sum(Normal_trials_summary([1,4,7,10,13],1));
        Normal_trials_summary(16,w+1) = sum(Normal_trials_summary([1,4,7,10,13],w+1));
        % Medium_dB
        Normal_trials_summary(17,1) = sum(Normal_trials_summary([2,5,8,11,14],1));
        Normal_trials_summary(17,w+1) = sum(Normal_trials_summary([2,5,8,11,14],w+1));
        % Loudest
        Normal_trials_summary(18,1) = sum(Normal_trials_summary([3,6,9,12,15],1));
        Normal_trials_summary(18,w+1) = sum(Normal_trials_summary([3,6,9,12,15],w+1)); %#ok<*SAGROW>
        
        % Shortest
        Normal_trials_summary(19,1) = sum(Normal_trials_summary([1,2,3],1));
        Normal_trials_summary(19,w+1) = sum(Normal_trials_summary([1,2,3],w+1));
        % Short
        Normal_trials_summary(20,1) = sum(Normal_trials_summary([4,5,6],1));
        Normal_trials_summary(20,w+1) = sum(Normal_trials_summary([4,5,6],w+1));
        % Medium_ISI
        Normal_trials_summary(21,1) = sum(Normal_trials_summary([7,8,9],1));
        Normal_trials_summary(21,w+1) = sum(Normal_trials_summary([7,8,9],w+1));
        % Slow
        Normal_trials_summary(22,1) = sum(Normal_trials_summary([10,11,12],1));
        Normal_trials_summary(22,w+1) = sum(Normal_trials_summary([10,11,12],w+1));
        % Slowest
        Normal_trials_summary(23,1) = sum(Normal_trials_summary([13,14,15],1));
        Normal_trials_summary(23,w+1) = sum(Normal_trials_summary([13,14,15],w+1));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    % Delete empty rows in detailed axis names and data (commented out
    % since we have fixed blocks and sessions, so if there is an empty spot
    % is worth noticing)
    % rowNames_detailed = rowNames_detailed(any(~cellfun('isempty',rowNames_detailed), 1));
    % Normal_trials_detailed = Normal_trials_detailed(any(Normal_trials_detailed,2),:);
    % Create tables from data matrix and save
    eval(['LLR_' participant{p} '_Normal_trials_detailed = array2table(Normal_trials_detailed,''RowNames'',rowNames_detailed,''VariableNames'',colNames);'])
    eval(['LLR_' participant{p} '_Normal_trials_summary = array2table(Normal_trials_summary,''RowNames'',rowNames_summary,''VariableNames'',colNames);'])
    save([root_dir '/Events/LLR_Counts_' participant{p} '_Normal_detailed.mat'],['LLR_' participant{p} '_Normal_trials_detailed']);
    save([root_dir '/Events/LLR_Counts_' participant{p} '_Normal_summary.mat'],['LLR_' participant{p} '_Normal_trials_summary']);
end

% reset modality_data variable to its original
modality_data = {'EEG','MEG','both_mod'};

clearvars('-except', initialVars{:});
disp 'DONE WITH COUNTING TRIALS FOR NORMAL EPOCHS FOR LLR!!!'
disp(datetime)
toc

%% Average waveforms from NORMAL EPOCHS across BLOCKS
% (even if sources will be computed individually for every block, 
% waveforms can already be averaged)
tic
disp(' ');      
disp('-------------------------');  
disp('AVERAGING WAVEFORMS ACROSS BLOCKS (NORMAL EPOCHS)'); 
disp(datetime)
disp('-------------------------');
disp(' ');

if sensor_analysis == 1
    modality_data = {'EEG'};
elseif sensor_analysis == 2
    modality_data = {'MEG'};
elseif sensor_analysis == 3
    modality_data = {'both_mod'};
end

for p = 1:length(participant)
    
    disp(' ');      
    disp('-------------------------');
    disp(['loading participant ' participant{p}]);
    disp(datetime)
    disp(' '); 
    
    prot_subs = bst_get('ProtocolSubjects');
    current_sub = find(strcmp({prot_subs.Subject.Name}, participant{p}));
    db_reload_conditions(current_sub);
    
    % Define sessions within this participant
    pos_par = find(strcmp(session_block_array(:,1), participant{p}(1:4)));
    session = session_block_array{pos_par,2}(1,:);
    
    for s = 1:length(session)
        
        % Define blocks within this session
        pos_ses = find(strcmp(session_block_array{pos_par,2}(1,:), session{s}));
        block = session_block_array{pos_par,2}{2,pos_ses}(1,:);
        
        for w = 1:length(wave)
            for c = 1:length(condition)
                folders = dir([root_dir_bs '/data/' participant{p} '/']);
                results =  contains({folders.name},[wave{w} '_' condition{c} session{s}]) & contains({folders.name},'_normal');
                infolder = find(results); % So get all blocks you find for that session
                if isempty(infolder)
                    Non_averaged_across_blocks{naab,1} = [wave{w} '_' participant{p}  '_' condition{c} session{s} '_normal'];
                    naab = naab +1;
                end
                for mode = 1:length(modality_data) % EEG, MEG and both mod
                    file_names = {};
                    for l= 1:length(infolder) % for each block found
                        line = infolder(l);
                        % Check if that block is in the block list, if not, go to next block
                        init_block_name = regexp(folders(line).name, '_b', 'once');
                        fol_nam = contains(block, folders(line).name(init_block_name+1:init_block_name+2));
                        in_fol_nam = find(fol_nam); % index of position in logical
                        if isempty(in_fol_nam); continue; end % the block in the folder name is not present in the list, so go to next
                        % End checking if block exist in list
                        sub_folders = dir([root_dir_bs '/data/' participant{p} '/' folders(line).name '/']);
                        sub_results =  contains({sub_folders.name},[modality_data{mode} '_average_sensor']); % Only those averages that are interpolated
                        sub_infolder = find(sub_results);
                        if isempty(sub_infolder) % if average is not present despite the folder existing
                            continue
                        end
                        file_names{l} = [participant{p} '/' folders(line).name '/' sub_folders(sub_infolder).name]; % sub_results (using logicals)
                    end
                    file_names = file_names(~cellfun('isempty', file_names')); % to avoid empty cells
                    sFiles = file_names;

                    if isempty(sFiles)
                        % It could be that all blocks from a particular
                        % modality have no average because of rejected
                        % trials in all blocks of the session, in which
                        % casse better to continue
                        continue
                    end
                    
                    disp(' ');      
                    disp('-------------------------');  
                    disp(['Averaging normal epochs for ' participant{p} '_' modality_data{mode} '_' wave{w} '_' condition{c} session{s}]);
                    disp(datetime)
                    disp(' ');

                    % If stated, find and delete previous normal averages across blocks
                    if delete_previous_file == 1
                        % check if there is already a normal average for
                        % each session in intra folder and delete it
                        folders_delete = dir([root_dir_bs '/data/' participant{p} '/@intra/']);
                        % data_LLR_11B_average_200727_1912_Normal_EEG_LLR_11B.mat
                        results_delete =  endsWith({folders_delete.name}, ['Normal_' modality_data{mode} '_' wave{w} '_' condition{c} session{s} '.mat']); 
                        infolder_delete = find(results_delete);
                        if ~isempty(infolder_delete) % file exists, therefore delete it
                           delete([root_dir_bs '/data/' participant{p} '/@intra/' folders_delete(infolder_delete).name]);
                        end
                    end
                    
                    % Process: Weighted Average: Everything
                    sFiles = bst_process('CallProcess', 'process_average', sFiles, [], ...
                        'avgtype',       1, ...  % Everything
                        'avg_func',      1, ...  % Arithmetic average:  mean(x)
                        'weighted',      1, ...
                        'keepevents',    0);

                    % Process: Add tag
                    sFiles = bst_process('CallProcess', 'process_add_tag', sFiles, [], ...
                        'tag',           ['Normal_' modality_data{mode} '_' wave{w} '_' condition{c} session{s}], ...
                        'output',        2);  % Add to file name (1 to add a tag)

                    % Process: Set name
                    sFiles = bst_process('CallProcess', 'process_set_comment', sFiles, [], ...
                        'tag',           ['Normal_' modality_data{mode} '_' wave{w} '_' condition{c} session{s}], ...
                        'isindex',       1);
                end
            end
        end
    end
end

% reset modality_data variable to its original
modality_data = {'EEG','MEG','both_mod'};
save([root_dir '/Averages/Non_averaged_across_blocks.mat'],'Non_averaged_across_blocks');

clearvars('-except', initialVars{:});
disp 'DONE WITH AVERAGING WAVEFORMS ACROSS BLOCKS (NORMAL EPOCHS)!!!'
disp(datetime)
toc

%% Average waveforms from NORMAL EPOCHS across blocks from all sessions
tic
disp(' ');      
disp('-------------------------');  
disp('AVERAGING WAVEFORMS ACROSS BLOCKS FROM ALL SESSIONS (NORMAL EPOCHS)'); 
disp(datetime)
disp('-------------------------');
disp(' ');

if sensor_analysis == 1
    modality_data = {'EEG'};
elseif sensor_analysis == 2
    modality_data = {'MEG'};
elseif sensor_analysis == 3
    modality_data = {'both_mod'};
end

for p = 1:length(participant)
    
    disp(' ');      
    disp('-------------------------');
    disp(['loading participant ' participant{p}]);
    disp(datetime)
    disp(' '); 

    % reload subject normally
    prot_subs = bst_get('ProtocolSubjects');
    current_sub = find(strcmp({prot_subs.Subject.Name}, participant{p}));
    db_reload_conditions(current_sub);
    
    % Define sessions within this participant
    pos_par = find(strcmp(session_block_array(:,1), participant{p}(1:4)));
    session = session_block_array{pos_par,2}(1,:);
    
    for w = 1:length(wave)
        for c = 1:length(condition)
            for mode =1:length(modality_data)
            file_names = {}; % create one for each p, c and w, so that it includes all sessions and blocks
            list_order = 1;
            for s = 1:length(session)
                
                % Define blocks within this session
                pos_ses = find(strcmp(session_block_array{pos_par,2}(1,:), session{s}));
                block = session_block_array{pos_par,2}{2,pos_ses}(1,:);
                
                folders = dir([root_dir_bs '/data/' participant{p} '/']);
                results =  contains({folders.name},[wave{w} '_' condition{c} session{s}]) & contains({folders.name},'_normal');
                infolder = find(results); % So get all blocks you find for that session
                if isempty(infolder)
                    continue
                    % error(['no ' wave{w} '_' condition{c} session{s} '_normal blocks for ' participant{p}]);
                end
                for l= 1:length(infolder) % for each block found for that session
                    line = infolder(l);
                    % Check if that block is in the block list, if not, go to next block
                    init_block_name = regexp(folders(line).name, '_b', 'once');
                    fol_nam = contains(block, folders(line).name(init_block_name+1:init_block_name+2));
                    in_fol_nam = find(fol_nam); % index of position in logical
                    if isempty(in_fol_nam); continue; end % the block in the folder name is not present in the list, so go to next
                    % End checking if block exist in list
                    sub_folders = dir([root_dir_bs '/data/' participant{p} '/' folders(line).name '/']);
                    sub_results =  contains({sub_folders.name},[modality_data{mode} '_average_sensor']);
                    sub_infolder = find(sub_results);
                    if isempty(sub_infolder) % if average is not present despite the folder existing
                        continue
                    end
                    file_names{list_order} = [participant{p} '/' folders(line).name '/' sub_folders(sub_infolder).name]; % sub_results (using logicals)
                    list_order = list_order + 1;
                end
            end
            file_names = file_names(~cellfun('isempty', file_names')); % to avoid empty cells
            sFiles = file_names;

            if isempty(sFiles)
            error('for whatever reason, sFiles is empty, my friend...')
            end
            
            disp(' ');      
            disp('-------------------------');  
            disp(['Averaging blocks from all sessions (normal epochs) for ' participant{p} '_' modality_data{mode} '_' wave{w} '_' condition{c}]);
            disp(datetime)
            disp(' ');

            % If stated, find and delete previous normal averages across
            % blocks from all sessions
            if delete_previous_file == 1
                % check if there is already a normal average for all sessions
                folders_delete = dir([root_dir_bs '/data/' participant{p} '/@intra/']);
                % data_LLR_11B_average_200727_1912_Normal_EEG_LLR_11B.mat
                results_delete =  endsWith({folders_delete.name}, ['Normal_' modality_data{mode} '_' wave{w} '_' condition{c} '.mat']); 
                infolder_delete = find(results_delete);
                if ~isempty(infolder_delete) % file exists, therefore delete it
                   delete([root_dir_bs '/data/' participant{p} '/@intra/' folders_delete(infolder_delete).name]);
                end
            end

            % Process: Weighted Average: Everything
            sFiles = bst_process('CallProcess', 'process_average', sFiles, [], ...
                'avgtype',       1, ...  % Everything
                'avg_func',      1, ...  % Arithmetic average:  mean(x)
                'weighted',      1, ...
                'keepevents',    0);


            % Process: Add tag
            sFiles = bst_process('CallProcess', 'process_add_tag', sFiles, [], ...
                'tag',           ['Normal_' modality_data{mode} '_' wave{w} '_' condition{c}], ...
                'output',        2);  % Add to file name (1 to add a tag)

            % Process: Set name
            sFiles = bst_process('CallProcess', 'process_set_comment', sFiles, [], ...
                'tag',           ['Normal_' modality_data{mode} '_' wave{w} '_' condition{c}], ...
                'isindex',       1);   
            end
        end     
    end
end

% reset modality_data variable to its original
modality_data = {'EEG','MEG','both_mod'};

clearvars('-except', initialVars{:});
disp 'DONE WITH AVERAGING WAVEFORMS ACROSS BLOCKS FROM ALL SESSIONS (NORMAL EPOCHS)!!!'
disp(datetime)
toc

%% Average ISI and dB per each session (NORMAL epochs)
tic
disp(' ');      
disp('-------------------------');  
disp('OBTAINING ISI and DB AVERAGES PER SESSION (NORMAL EPOCHS)'); 
disp(datetime)
disp('-------------------------');
disp(' ');

if sensor_analysis == 1
    modality_data = {'EEG'};
elseif sensor_analysis == 2
    modality_data = {'MEG'};
elseif sensor_analysis == 3
    modality_data = {'both_mod'};
end

for p = 1:length(participant)
    
    disp(' ');      
    disp('-------------------------');
    disp(['loading participant ' participant{p}]);
    disp(datetime)
    disp(' '); 

    % reload subject normally
    prot_subs = bst_get('ProtocolSubjects');
    current_sub = find(strcmp({prot_subs.Subject.Name}, participant{p}));
    db_reload_conditions(current_sub);
    
    % Define sessions within this participant
    pos_par = find(strcmp(session_block_array(:,1), participant{p}(1:4)));
    session = session_block_array{pos_par,2}(1,:);
    for mode = 1:length(modality_data)
        for w = 1:length(wave)        
            for s = 1:length(session)
                
                % Define blocks within this session
                pos_ses = find(strcmp(session_block_array{pos_par,2}(1,:), session{s}));
                block = session_block_array{pos_par,2}{2,pos_ses}(1,:);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Quietest %%%%%%%%%%%%%%%%%%%%%%%%
                folders = dir([root_dir_bs '/data/' participant{p} '/']);
                results =  (contains({folders.name},[wave{w} '_11' session{s}])... 
                | contains({folders.name},[wave{w} '_31' session{s}]) ...
                | contains({folders.name},[wave{w} '_51' session{s}]) ...
                | contains({folders.name},[wave{w} '_71' session{s}]) ...
                | contains({folders.name},[wave{w} '_91' session{s}])) ...
                & contains({folders.name},'_normal');   
                infolder = find(results); % So get all blocks you find for that session and condition
                if isempty(infolder)
                    continue;
                end            

                for l= 1:length(infolder) % for each block found for that session
                    line = infolder(l);
                    % Check if that block is in the block list, if not, go to next block
                    init_block_name = regexp(folders(line).name, '_b', 'once');
                    fol_nam = contains(block, folders(line).name(init_block_name+1:init_block_name+2));
                    in_fol_nam = find(fol_nam); % index of position in logical
                    if isempty(in_fol_nam); continue; end % the block in the folder name is not present in the list, so go to next
                    % End checking if block exist in list
                    sub_folders = dir([root_dir_bs '/data/' participant{p} '/' folders(line).name '/']);
                    sub_results =  contains({sub_folders.name},[modality_data{mode} '_average_sensor']);
                    sub_infolder = find(sub_results);
                    if isempty(sub_infolder) % if average is not present despite the folder existing
                        continue
                    end
                    file_names_Quietest{l} = [participant{p} '/' folders(line).name '/' sub_folders(sub_infolder).name]; % sub_results (using logicals)
                end        
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Medium_dB %%%%%%%%%%%%%%%%%%%%%%%%
                folders = dir([root_dir_bs '/data/' participant{p} '/']);
                results =  (contains({folders.name},[wave{w} '_12' session{s}])... 
                | contains({folders.name},[wave{w} '_32' session{s}]) ...
                | contains({folders.name},[wave{w} '_52' session{s}]) ...
                | contains({folders.name},[wave{w} '_72' session{s}]) ...
                | contains({folders.name},[wave{w} '_92' session{s}])) ...
                & contains({folders.name},'_normal');   
                infolder = find(results); % So get all blocks you find for that session and condition
                if isempty(infolder)
                    error('no coincidences');
                end            

                for l= 1:length(infolder) % for each block found for that session
                    line = infolder(l);
                    % Check if that block is in the block list, if not, go to next block
                    init_block_name = regexp(folders(line).name, '_b', 'once');
                    fol_nam = contains(block, folders(line).name(init_block_name+1:init_block_name+2));
                    in_fol_nam = find(fol_nam); % index of position in logical
                    if isempty(in_fol_nam); continue; end % the block in the folder name is not present in the list, so go to next
                    % End checking if block exist in list
                    sub_folders = dir([root_dir_bs '/data/' participant{p} '/' folders(line).name '/']);
                    sub_results =  contains({sub_folders.name},[modality_data{mode} '_average_sensor']);
                    sub_infolder = find(sub_results);
                    if isempty(sub_infolder) % if average is not present despite the folder existing
                        continue
                    end
                    file_names_Medium_dB{l} = [participant{p} '/' folders(line).name '/' sub_folders(sub_infolder).name]; % sub_results (using logicals)
                end        
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Loudest %%%%%%%%%%%%%%%%%%%%%%%%
                folders = dir([root_dir_bs '/data/' participant{p} '/']);
                results =  (contains({folders.name},[wave{w} '_13' session{s}])... 
                | contains({folders.name},[wave{w} '_33' session{s}]) ...
                | contains({folders.name},[wave{w} '_53' session{s}]) ...
                | contains({folders.name},[wave{w} '_73' session{s}]) ...
                | contains({folders.name},[wave{w} '_93' session{s}])) ...
                & contains({folders.name},'_normal');   
                infolder = find(results); % So get all blocks you find for that session and condition
                if isempty(infolder)
                    error('no coincidences');
                end            

                for l= 1:length(infolder) % for each block found for that session
                    line = infolder(l);
                    % Check if that block is in the block list, if not, go to next block
                    init_block_name = regexp(folders(line).name, '_b', 'once');
                    fol_nam = contains(block, folders(line).name(init_block_name+1:init_block_name+2));
                    in_fol_nam = find(fol_nam); % index of position in logical
                    if isempty(in_fol_nam); continue; end % the block in the folder name is not present in the list, so go to next
                    % End checking if block exist in list
                    sub_folders = dir([root_dir_bs '/data/' participant{p} '/' folders(line).name '/']);
                    sub_results =  contains({sub_folders.name},[modality_data{mode} '_average_sensor']);
                    sub_infolder = find(sub_results);
                    if isempty(sub_infolder) % if average is not present despite the folder existing
                        continue
                    end
                    file_names_Loudest{l} = [participant{p} '/' folders(line).name '/' sub_folders(sub_infolder).name]; % sub_results (using logicals)
                end        
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fastest %%%%%%%%%%%%%%%%%%%%%%%%
                folders = dir([root_dir_bs '/data/' participant{p} '/']);
                results =  (contains({folders.name},[wave{w} '_11' session{s}])... 
                | contains({folders.name},[wave{w} '_12' session{s}]) ...
                | contains({folders.name},[wave{w} '_13' session{s}])) ...
                & contains({folders.name},'_normal');   
                infolder = find(results); % So get all blocks you find for that session and condition
                if isempty(infolder)
                    error('no coincidences');
                end            

                for l= 1:length(infolder) % for each block found for that session
                    line = infolder(l);
                    % Check if that block is in the block list, if not, go to next block
                    init_block_name = regexp(folders(line).name, '_b', 'once');
                    fol_nam = contains(block, folders(line).name(init_block_name+1:init_block_name+2));
                    in_fol_nam = find(fol_nam); % index of position in logical
                    if isempty(in_fol_nam); continue; end % the block in the folder name is not present in the list, so go to next
                    % End checking if block exist in list
                    sub_folders = dir([root_dir_bs '/data/' participant{p} '/' folders(line).name '/']);
                    sub_results =  contains({sub_folders.name},[modality_data{mode} '_average_sensor']);
                    sub_infolder = find(sub_results);
                    if isempty(sub_infolder) % if average is not present despite the folder existing
                        continue
                    end
                    file_names_Fastest{l} = [participant{p} '/' folders(line).name '/' sub_folders(sub_infolder).name]; % sub_results (using logicals)
                end        
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fast %%%%%%%%%%%%%%%%%%%%%%%%
                folders = dir([root_dir_bs '/data/' participant{p} '/']);
                results =  (contains({folders.name},[wave{w} '_31' session{s}])... 
                | contains({folders.name},[wave{w} '_32' session{s}]) ...
                | contains({folders.name},[wave{w} '_33' session{s}])) ...
                & contains({folders.name},'_normal');   
                infolder = find(results); % So get all blocks you find for that session and condition
                if isempty(infolder)
                    error('no coincidences');
                end            

                for l= 1:length(infolder) % for each block found for that session
                    line = infolder(l);
                    % Check if that block is in the block list, if not, go to next block
                    init_block_name = regexp(folders(line).name, '_b', 'once');
                    fol_nam = contains(block, folders(line).name(init_block_name+1:init_block_name+2));
                    in_fol_nam = find(fol_nam); % index of position in logical
                    if isempty(in_fol_nam); continue; end % the block in the folder name is not present in the list, so go to next
                    % End checking if block exist in list
                    sub_folders = dir([root_dir_bs '/data/' participant{p} '/' folders(line).name '/']);
                    sub_results =  contains({sub_folders.name},[modality_data{mode} '_average_sensor']);
                    sub_infolder = find(sub_results);
                    if isempty(sub_infolder) % if average is not present despite the folder existing
                        continue
                    end
                    file_names_Fast{l} = [participant{p} '/' folders(line).name '/' sub_folders(sub_infolder).name]; % sub_results (using logicals)
                end        
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Medium_ISI %%%%%%%%%%%%%%%%%%%%%%%%
                folders = dir([root_dir_bs '/data/' participant{p} '/']);
                results =  (contains({folders.name},[wave{w} '_51' session{s}])... 
                | contains({folders.name},[wave{w} '_52' session{s}]) ...
                | contains({folders.name},[wave{w} '_53' session{s}])) ...
                & contains({folders.name},'_normal');   
                infolder = find(results); % So get all blocks you find for that session and condition
                if isempty(infolder)
                    error('no coincidences');
                end            

                for l= 1:length(infolder) % for each block found for that session
                    line = infolder(l);
                    % Check if that block is in the block list, if not, go to next block
                    init_block_name = regexp(folders(line).name, '_b', 'once');
                    fol_nam = contains(block, folders(line).name(init_block_name+1:init_block_name+2));
                    in_fol_nam = find(fol_nam); % index of position in logical
                    if isempty(in_fol_nam); continue; end % the block in the folder name is not present in the list, so go to next
                    % End checking if block exist in list
                    sub_folders = dir([root_dir_bs '/data/' participant{p} '/' folders(line).name '/']);
                    sub_results =  contains({sub_folders.name},[modality_data{mode} '_average_sensor']);
                    sub_infolder = find(sub_results);
                    if isempty(sub_infolder) % if average is not present despite the folder existing
                        continue
                    end
                    file_names_Medium_ISI{l} = [participant{p} '/' folders(line).name '/' sub_folders(sub_infolder).name]; % sub_results (using logicals)
                end        
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Slow %%%%%%%%%%%%%%%%%%%%%%%%
                folders = dir([root_dir_bs '/data/' participant{p} '/']);
                results =  (contains({folders.name},[wave{w} '_71' session{s}])... 
                | contains({folders.name},[wave{w} '_72' session{s}]) ...
                | contains({folders.name},[wave{w} '_73' session{s}])) ...
                & contains({folders.name},'_normal');   
                infolder = find(results); % So get all blocks you find for that session and condition
                if isempty(infolder)
                    error('no coincidences');
                end            

                for l= 1:length(infolder) % for each block found for that session
                    line = infolder(l);
                    % Check if that block is in the block list, if not, go to next block
                    init_block_name = regexp(folders(line).name, '_b', 'once');
                    fol_nam = contains(block, folders(line).name(init_block_name+1:init_block_name+2));
                    in_fol_nam = find(fol_nam); % index of position in logical
                    if isempty(in_fol_nam); continue; end % the block in the folder name is not present in the list, so go to next
                    % End checking if block exist in list
                    sub_folders = dir([root_dir_bs '/data/' participant{p} '/' folders(line).name '/']);
                    sub_results =  contains({sub_folders.name},[modality_data{mode} '_average_sensor']);
                    sub_infolder = find(sub_results);
                    if isempty(sub_infolder) % if average is not present despite the folder existing
                        continue
                    end
                    file_names_Slow{l} = [participant{p} '/' folders(line).name '/' sub_folders(sub_infolder).name]; % sub_results (using logicals)
                end        
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Slowest %%%%%%%%%%%%%%%%%%%%%%%%
                folders = dir([root_dir_bs '/data/' participant{p} '/']);
                results =  (contains({folders.name},[wave{w} '_91' session{s}])... 
                | contains({folders.name},[wave{w} '_92' session{s}]) ...
                | contains({folders.name},[wave{w} '_93' session{s}])) ...
                & contains({folders.name},'_normal');   
                infolder = find(results); % So get all blocks you find for that session and condition
                if isempty(infolder)
                    error('no coincidences');
                end            

                for l= 1:length(infolder) % for each block found for that session
                    line = infolder(l);
                    % Check if that block is in the block list, if not, go to next block
                    init_block_name = regexp(folders(line).name, '_b', 'once');
                    fol_nam = contains(block, folders(line).name(init_block_name+1:init_block_name+2));
                    in_fol_nam = find(fol_nam); % index of position in logical
                    if isempty(in_fol_nam); continue; end % the block in the folder name is not present in the list, so go to next
                    % End checking if block exist in list
                    sub_folders = dir([root_dir_bs '/data/' participant{p} '/' folders(line).name '/']);
                    sub_results =  contains({sub_folders.name},[modality_data{mode} '_average_sensor']);
                    sub_infolder = find(sub_results);
                    if isempty(sub_infolder) % if average is not present despite the folder existing
                        continue
                    end
                    file_names_Slowest{l} = [participant{p} '/' folders(line).name '/' sub_folders(sub_infolder).name]; % sub_results (using logicals)
                end        
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                file_names_Quietest = file_names_Quietest(~cellfun('isempty', file_names_Quietest'));
                file_names_Medium_dB = file_names_Medium_dB(~cellfun('isempty', file_names_Medium_dB'));
                file_names_Loudest = file_names_Loudest(~cellfun('isempty', file_names_Loudest'));
                file_names_Fastest = file_names_Fastest(~cellfun('isempty', file_names_Fastest'));
                file_names_Fast = file_names_Fast(~cellfun('isempty', file_names_Fast'));
                file_names_Medium_ISI = file_names_Medium_ISI(~cellfun('isempty', file_names_Medium_ISI'));
                file_names_Slow = file_names_Slow(~cellfun('isempty', file_names_Slow'));
                file_names_Slowest = file_names_Slowest(~cellfun('isempty', file_names_Slowest'));

                for ex = 1:length(Exp_cond)
                    eval(['sFiles = file_names_' Exp_cond{ex} ';'])

                    if isempty(sFiles)
                        % It could be that all blocks from a particular
                        % modality have no average because of rejected
                        % trials in all blocks of the session, in which
                        % case better to continue
                        continue
                    end
                    
                    disp(' ');      
                    disp('-------------------------');  
                    disp(['Averaging dB and ISI across blocks (normal epochs) for ' participant{p} '_' modality_data{mode} '_' wave{w} '_' Exp_cond{ex} '_' session{s}]);
                    disp(datetime)
                    disp(' ');

                    % If stated, find and delete previous normal averages
                    if delete_previous_file == 1
                        % check if there is already a normal average for this session
                        folders_delete = dir([root_dir_bs '/data/' participant{p} '/@intra/']);
                        % data_LLR_average_200727_2009_Normal_both_mod_LLR_Fastest_A.mat
                        results_delete =  endsWith({folders_delete.name}, ['Normal_' modality_data{mode} '_' wave{w} '_' Exp_cond{ex} '_' session{s} '.mat']); 
                        infolder_delete = find(results_delete);
                        if ~isempty(infolder_delete) % file exists, therefore delete it
                           delete([root_dir_bs '/data/' participant{p} '/@intra/' folders_delete(infolder_delete).name]);
                        end
                    end
                    
                    % Process: Weighted Average: Everything
                    sFiles = bst_process('CallProcess', 'process_average', sFiles, [], ...
                        'avgtype',       1, ...  % Everything
                        'avg_func',      1, ...  % Arithmetic average:  mean(x)
                        'weighted',      1, ...
                        'keepevents',    0);

                    % Process: Add tag
                    sFiles = bst_process('CallProcess', 'process_add_tag', sFiles, [], ...
                        'tag',           ['Normal_' modality_data{mode} '_' wave{w} '_' Exp_cond{ex} '_' session{s}], ...
                        'output',        2);  % Add to file name (1 to add a tag)

                    % Process: Set name
                    sFiles = bst_process('CallProcess', 'process_set_comment', sFiles, [], ...
                        'tag',           ['Normal_' modality_data{mode} '_' wave{w} '_' Exp_cond{ex} '_' session{s}], ...
                        'isindex',       1);   
                end
            end
        end
    end
end

% reset modality_data variable to its original
modality_data = {'EEG','MEG','both_mod'};

clearvars('-except', initialVars{:});
disp 'DONE WITH AVERAGING dB AND ISI CONDITIONS PER SESSION (NORMAL EPOCHS)!!!'
disp(datetime)
toc

%% Average ALL conditions per each session (NORMAL epochs)
tic
disp(' ');      
disp('-------------------------');  
disp('OBTAINING ALL_COND AVERAGES PER SESSION (NORMAL EPOCHS)'); 
disp(datetime)
disp('-------------------------');
disp(' ');

if sensor_analysis == 1
    modality_data = {'EEG'};
elseif sensor_analysis == 2
    modality_data = {'MEG'};
elseif sensor_analysis == 3
    modality_data = {'both_mod'};
end

for p = 1:length(participant)
    
    disp(' ');      
    disp('-------------------------');
    disp(['loading participant ' participant{p}]);
    disp(datetime)
    disp(' '); 

    % Reload subject normally
    prot_subs = bst_get('ProtocolSubjects');
    current_sub = find(strcmp({prot_subs.Subject.Name}, participant{p}));
    db_reload_conditions(current_sub);
    
    % Define sessions within this participant
    pos_par = find(strcmp(session_block_array(:,1), participant{p}(1:4)));
    session = session_block_array{pos_par,2}(1,:);
    for mode = 1:length(modality_data)
        for w = 1:length(wave)        
            for s = 1:length(session)
                
                % Define blocks within this session
                pos_ses = find(strcmp(session_block_array{pos_par,2}(1,:), session{s}));
                block = session_block_array{pos_par,2}{2,pos_ses}(1,:);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ALL %%%%%%%%%%%%%%%%%%%%%%%%
                folders = dir([root_dir_bs '/data/' participant{p} '/']);
                results =  (contains({folders.name},[wave{w} '_11' session{s}])... 
                | contains({folders.name},[wave{w} '_12' session{s}]) ...
                | contains({folders.name},[wave{w} '_13' session{s}]) ...
                | contains({folders.name},[wave{w} '_31' session{s}]) ...
                | contains({folders.name},[wave{w} '_32' session{s}]) ...
                | contains({folders.name},[wave{w} '_33' session{s}]) ...
                | contains({folders.name},[wave{w} '_51' session{s}]) ...
                | contains({folders.name},[wave{w} '_52' session{s}]) ...
                | contains({folders.name},[wave{w} '_53' session{s}]) ...
                | contains({folders.name},[wave{w} '_71' session{s}]) ...
                | contains({folders.name},[wave{w} '_72' session{s}]) ...
                | contains({folders.name},[wave{w} '_73' session{s}]) ...
                | contains({folders.name},[wave{w} '_91' session{s}]) ...
                | contains({folders.name},[wave{w} '_92' session{s}]) ...
                | contains({folders.name},[wave{w} '_93' session{s}])) ...
                & contains({folders.name},'_normal');   
                infolder = find(results); % So get all blocks you find for that session and gather all conditions
                if isempty(infolder)
                    continue;
                end            

                for l= 1:length(infolder) % for each block found for that session
                    line = infolder(l);
                    % Check if that block is in the block list, if not, go to next block
                    init_block_name = regexp(folders(line).name, '_b', 'once');
                    fol_nam = contains(block, folders(line).name(init_block_name+1:init_block_name+2));
                    in_fol_nam = find(fol_nam); % index of position in logical
                    if isempty(in_fol_nam); continue; end % the block in the folder name is not present in the list, so go to next
                    % End checking if block exist in list
                    sub_folders = dir([root_dir_bs '/data/' participant{p} '/' folders(line).name '/']);
                    sub_results =  contains({sub_folders.name},[modality_data{mode} '_average_sensor']);
                    sub_infolder = find(sub_results);
                    if isempty(sub_infolder) % if average is not present despite the folder existing
                        continue
                    end
                    file_names_ALL{l} = [participant{p} '/' folders(line).name '/' sub_folders(sub_infolder).name]; % sub_results (using logicals)
                end        
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                file_names_ALL = file_names_ALL(~cellfun('isempty', file_names_ALL'));
                sFiles = file_names_ALL;
                
                if isempty(sFiles)
                    % It could be that all blocks from a particular
                    % modality have no average because of rejected
                    % trials in all blocks of the session, in which
                    % case better to continue
                    continue
                end

                disp(' ');      
                disp('-------------------------');  
                disp(['Averaging ALL conditions across blocks (normal epochs) for ' participant{p} '_' modality_data{mode} '_' wave{w} '_' session{s}]);
                disp(datetime)
                disp(' ');

                % If stated, find and delete previous normal averages
                if delete_previous_file == 1
                    % check if there is already a normal average for this session
                    folders_delete = dir([root_dir_bs '/data/' participant{p} '/@intra/']);
                    % data_LLR_average_200727_2009_Normal_both_mod_LLR_Fastest_A.mat
                    results_delete =  endsWith({folders_delete.name}, ['Normal_' modality_data{mode} '_' wave{w} '_ALL_' session{s} '.mat']); 
                    infolder_delete = find(results_delete);
                    if ~isempty(infolder_delete) % file exists, therefore delete it
                       delete([root_dir_bs '/data/' participant{p} '/@intra/' folders_delete(infolder_delete).name]);
                    end
                end

                % Process: Weighted Average: Everything
                sFiles = bst_process('CallProcess', 'process_average', sFiles, [], ...
                    'avgtype',       1, ...  % Everything
                    'avg_func',      1, ...  % Arithmetic average:  mean(x)
                    'weighted',      1, ...
                    'keepevents',    0);

                % Process: Add tag
                sFiles = bst_process('CallProcess', 'process_add_tag', sFiles, [], ...
                    'tag',           ['Normal_' modality_data{mode} '_' wave{w} '_ALL_' session{s}], ...
                    'output',        2);  % Add to file name (1 to add a tag)

                % Process: Set name
                sFiles = bst_process('CallProcess', 'process_set_comment', sFiles, [], ...
                    'tag',           ['Normal_' modality_data{mode} '_' wave{w} '_ALL_' session{s}], ...
                    'isindex',       1);   
            end
        end
    end
end

% reset modality_data variable to its original
modality_data = {'EEG','MEG','both_mod'};

clearvars('-except', initialVars{:});
disp 'DONE OBTAINING ALL_COND AVERAGES PER SESSION (NORMAL EPOCHS)!!!'
disp(datetime)
toc

%% Average across ISI and dB within subjects (NORMAL epochs always)
tic
disp(' ');      
disp('-------------------------');  
disp('AVERAGING WAVEFORMS ACROSS dB AND ISI CONDITIONS (all blocks and sessions)'); 
disp(datetime)
disp('-------------------------');
disp(' ');

if sensor_analysis == 1
    modality_data = {'EEG'};
elseif sensor_analysis == 2
    modality_data = {'MEG'};
elseif sensor_analysis == 3
    modality_data = {'both_mod'};
end

for p = 1:length(participant)
    
    disp(' ');      
    disp('-------------------------');
    disp(['loading participant ' participant{p}]);
    disp(datetime)
    disp(' '); 

    % reload subject normally
    prot_subs = bst_get('ProtocolSubjects');
    current_sub = find(strcmp({prot_subs.Subject.Name}, participant{p}));
    db_reload_conditions(current_sub);
    
    % Define sessions within this participant
    pos_par = find(strcmp(session_block_array(:,1), participant{p}(1:4)));
    session = session_block_array{pos_par,2}(1,:);
    for mode = 1:length(modality_data)
        for w = 1:length(wave)        
            file_names_Quietest = {}; % create one for each p and w, so that it includes all sessions and blocks
            list_order_Quietest = 1;
            file_names_Medium_dB = {};
            list_order_Medium_dB = 1;
            file_names_Loudest = {};
            list_order_Loudest = 1;
            file_names_Fastest= {};
            list_order_Fastest = 1;
            file_names_Fast = {};
            list_order_Fast = 1;
            file_names_Medium_ISI = {};
            list_order_Medium_ISI = 1;
            file_names_Slow = {};
            list_order_Slow= 1;
            file_names_Slowest = {};
            list_order_Slowest = 1;
            for s = 1:length(session)
                
                % Define blocks within this session
                pos_ses = find(strcmp(session_block_array{pos_par,2}(1,:), session{s}));
                block = session_block_array{pos_par,2}{2,pos_ses}(1,:);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Quietest %%%%%%%%%%%%%%%%%%%%%%%%
                folders = dir([root_dir_bs '/data/' participant{p} '/']);
                results =  (contains({folders.name},[wave{w} '_11' session{s}])... 
                | contains({folders.name},[wave{w} '_31' session{s}]) ...
                | contains({folders.name},[wave{w} '_51' session{s}]) ...
                | contains({folders.name},[wave{w} '_71' session{s}]) ...
                | contains({folders.name},[wave{w} '_91' session{s}])) ...
                & contains({folders.name},'_normal');   
                infolder = find(results); % So get all blocks you find for that session and condition
                if isempty(infolder)
                    continue
                    % error('no coincidences');
                end            

                for l= 1:length(infolder) % for each block found for that session
                    line = infolder(l);
                    % Check if that block is in the block list, if not, go to next block
                    init_block_name = regexp(folders(line).name, '_b', 'once');
                    fol_nam = contains(block, folders(line).name(init_block_name+1:init_block_name+2));
                    in_fol_nam = find(fol_nam); % index of position in logical
                    if isempty(in_fol_nam); continue; end % the block in the folder name is not present in the list, so go to next
                    % End checking if block exist in list
                    sub_folders = dir([root_dir_bs '/data/' participant{p} '/' folders(line).name '/']);
                    sub_results =  contains({sub_folders.name},[modality_data{mode} '_average_sensor']);
                    sub_infolder = find(sub_results);
                    if isempty(sub_infolder) % if average is not present despite the folder existing
                        continue
                    end
                    file_names_Quietest{list_order_Quietest} = [participant{p} '/' folders(line).name '/' sub_folders(sub_infolder).name]; % sub_results (using logicals)
                    list_order_Quietest = list_order_Quietest + 1;
                end        
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Medium_dB %%%%%%%%%%%%%%%%%%%%%%%%
                folders = dir([root_dir_bs '/data/' participant{p} '/']);
                results =  (contains({folders.name},[wave{w} '_12' session{s}])... 
                | contains({folders.name},[wave{w} '_32' session{s}]) ...
                | contains({folders.name},[wave{w} '_52' session{s}]) ...
                | contains({folders.name},[wave{w} '_72' session{s}]) ...
                | contains({folders.name},[wave{w} '_92' session{s}])) ...
                & contains({folders.name},'_normal');   
                infolder = find(results); % So get all blocks you find for that session and condition
                if isempty(infolder)
                    error('no coincidences');
                end            

                for l= 1:length(infolder) % for each block found for that session
                    line = infolder(l);
                    % Check if that block is in the block list, if not, go to next block
                    init_block_name = regexp(folders(line).name, '_b', 'once');
                    fol_nam = contains(block, folders(line).name(init_block_name+1:init_block_name+2));
                    in_fol_nam = find(fol_nam); % index of position in logical
                    if isempty(in_fol_nam); continue; end % the block in the folder name is not present in the list, so go to next
                    % End checking if block exist in list
                    sub_folders = dir([root_dir_bs '/data/' participant{p} '/' folders(line).name '/']);
                    sub_results =  contains({sub_folders.name},[modality_data{mode} '_average_sensor']);
                    sub_infolder = find(sub_results);
                    if isempty(sub_infolder) % if average is not present despite the folder existing
                        continue
                    end
                    file_names_Medium_dB{list_order_Medium_dB} = [participant{p} '/' folders(line).name '/' sub_folders(sub_infolder).name]; % sub_results (using logicals)
                    list_order_Medium_dB = list_order_Medium_dB + 1;
                end        
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Loudest %%%%%%%%%%%%%%%%%%%%%%%%
                folders = dir([root_dir_bs '/data/' participant{p} '/']);
                results =  (contains({folders.name},[wave{w} '_13' session{s}])... 
                | contains({folders.name},[wave{w} '_33' session{s}]) ...
                | contains({folders.name},[wave{w} '_53' session{s}]) ...
                | contains({folders.name},[wave{w} '_73' session{s}]) ...
                | contains({folders.name},[wave{w} '_93' session{s}])) ...
                & contains({folders.name},'_normal');   
                infolder = find(results); % So get all blocks you find for that session and condition
                if isempty(infolder)
                    error('no coincidences');
                end            

                for l= 1:length(infolder) % for each block found for that session
                    line = infolder(l);
                    % Check if that block is in the block list, if not, go to next block
                    init_block_name = regexp(folders(line).name, '_b', 'once');
                    fol_nam = contains(block, folders(line).name(init_block_name+1:init_block_name+2));
                    in_fol_nam = find(fol_nam); % index of position in logical
                    if isempty(in_fol_nam); continue; end % the block in the folder name is not present in the list, so go to next
                    % End checking if block exist in list
                    sub_folders = dir([root_dir_bs '/data/' participant{p} '/' folders(line).name '/']);
                    sub_results =  contains({sub_folders.name},[modality_data{mode} '_average_sensor']);
                    sub_infolder = find(sub_results);
                    if isempty(sub_infolder) % if average is not present despite the folder existing
                        continue
                    end
                    file_names_Loudest{list_order_Loudest} = [participant{p} '/' folders(line).name '/' sub_folders(sub_infolder).name]; % sub_results (using logicals)
                    list_order_Loudest = list_order_Loudest + 1;
                end        
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fastest %%%%%%%%%%%%%%%%%%%%%%%%
                folders = dir([root_dir_bs '/data/' participant{p} '/']);
                results =  (contains({folders.name},[wave{w} '_11' session{s}])... 
                | contains({folders.name},[wave{w} '_12' session{s}]) ...
                | contains({folders.name},[wave{w} '_13' session{s}])) ...
                & contains({folders.name},'_normal');   
                infolder = find(results); % So get all blocks you find for that session and condition
                if isempty(infolder)
                    error('no coincidences');
                end            

                for l= 1:length(infolder) % for each block found for that session
                    line = infolder(l);
                    % Check if that block is in the block list, if not, go to next block
                    init_block_name = regexp(folders(line).name, '_b', 'once');
                    fol_nam = contains(block, folders(line).name(init_block_name+1:init_block_name+2));
                    in_fol_nam = find(fol_nam); % index of position in logical
                    if isempty(in_fol_nam); continue; end % the block in the folder name is not present in the list, so go to next
                    % End checking if block exist in list
                    sub_folders = dir([root_dir_bs '/data/' participant{p} '/' folders(line).name '/']);
                    sub_results =  contains({sub_folders.name},[modality_data{mode} '_average_sensor']);
                    sub_infolder = find(sub_results);
                    if isempty(sub_infolder) % if average is not present despite the folder existing
                        continue
                    end
                    file_names_Fastest{list_order_Fastest} = [participant{p} '/' folders(line).name '/' sub_folders(sub_infolder).name]; % sub_results (using logicals)
                    list_order_Fastest = list_order_Fastest + 1;
                end        
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fast %%%%%%%%%%%%%%%%%%%%%%%%
                folders = dir([root_dir_bs '/data/' participant{p} '/']);
                results =  (contains({folders.name},[wave{w} '_31' session{s}])... 
                | contains({folders.name},[wave{w} '_32' session{s}]) ...
                | contains({folders.name},[wave{w} '_33' session{s}])) ...
                & contains({folders.name},'_normal');   
                infolder = find(results); % So get all blocks you find for that session and condition
                if isempty(infolder)
                    error('no coincidences');
                end            

                for l= 1:length(infolder) % for each block found for that session
                    line = infolder(l);
                    % Check if that block is in the block list, if not, go to next block
                    init_block_name = regexp(folders(line).name, '_b', 'once');
                    fol_nam = contains(block, folders(line).name(init_block_name+1:init_block_name+2));
                    in_fol_nam = find(fol_nam); % index of position in logical
                    if isempty(in_fol_nam); continue; end % the block in the folder name is not present in the list, so go to next
                    % End checking if block exist in list
                    sub_folders = dir([root_dir_bs '/data/' participant{p} '/' folders(line).name '/']);
                    sub_results =  contains({sub_folders.name},[modality_data{mode} '_average_sensor']);
                    sub_infolder = find(sub_results);
                    if isempty(sub_infolder) % if average is not present despite the folder existing
                        continue
                    end
                    file_names_Fast{list_order_Fast} = [participant{p} '/' folders(line).name '/' sub_folders(sub_infolder).name]; % sub_results (using logicals)
                    list_order_Fast = list_order_Fast + 1;
                end        
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Medium_ISI %%%%%%%%%%%%%%%%%%%%%%%%
                folders = dir([root_dir_bs '/data/' participant{p} '/']);
                results =  (contains({folders.name},[wave{w} '_51' session{s}])... 
                | contains({folders.name},[wave{w} '_52' session{s}]) ...
                | contains({folders.name},[wave{w} '_53' session{s}])) ...
                & contains({folders.name},'_normal');   
                infolder = find(results); % So get all blocks you find for that session and condition
                if isempty(infolder)
                    error('no coincidences');
                end            

                for l= 1:length(infolder) % for each block found for that session
                    line = infolder(l);
                    % Check if that block is in the block list, if not, go to next block
                    init_block_name = regexp(folders(line).name, '_b', 'once');
                    fol_nam = contains(block, folders(line).name(init_block_name+1:init_block_name+2));
                    in_fol_nam = find(fol_nam); % index of position in logical
                    if isempty(in_fol_nam); continue; end % the block in the folder name is not present in the list, so go to next
                    % End checking if block exist in list
                    sub_folders = dir([root_dir_bs '/data/' participant{p} '/' folders(line).name '/']);
                    sub_results =  contains({sub_folders.name},[modality_data{mode} '_average_sensor']);
                    sub_infolder = find(sub_results);
                    if isempty(sub_infolder) % if average is not present despite the folder existing
                        continue
                    end
                    file_names_Medium_ISI{list_order_Medium_ISI} = [participant{p} '/' folders(line).name '/' sub_folders(sub_infolder).name]; % sub_results (using logicals)
                    list_order_Medium_ISI = list_order_Medium_ISI + 1;
                end        
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Slow %%%%%%%%%%%%%%%%%%%%%%%%
                folders = dir([root_dir_bs '/data/' participant{p} '/']);
                results =  (contains({folders.name},[wave{w} '_71' session{s}])... 
                | contains({folders.name},[wave{w} '_72' session{s}]) ...
                | contains({folders.name},[wave{w} '_73' session{s}])) ...
                & contains({folders.name},'_normal');   
                infolder = find(results); % So get all blocks you find for that session and condition
                if isempty(infolder)
                    error('no coincidences');
                end            

                for l= 1:length(infolder) % for each block found for that session
                    line = infolder(l);
                    % Check if that block is in the block list, if not, go to next block
                    init_block_name = regexp(folders(line).name, '_b', 'once');
                    fol_nam = contains(block, folders(line).name(init_block_name+1:init_block_name+2));
                    in_fol_nam = find(fol_nam); % index of position in logical
                    if isempty(in_fol_nam); continue; end % the block in the folder name is not present in the list, so go to next
                    % End checking if block exist in list
                    sub_folders = dir([root_dir_bs '/data/' participant{p} '/' folders(line).name '/']);
                    sub_results =  contains({sub_folders.name},[modality_data{mode} '_average_sensor']);
                    sub_infolder = find(sub_results);
                    if isempty(sub_infolder) % if average is not present despite the folder existing
                        continue
                    end
                    file_names_Slow{list_order_Slow} = [participant{p} '/' folders(line).name '/' sub_folders(sub_infolder).name]; % sub_results (using logicals)
                    list_order_Slow = list_order_Slow + 1;
                end        
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Slowest %%%%%%%%%%%%%%%%%%%%%%%%
                folders = dir([root_dir_bs '/data/' participant{p} '/']);
                results =  (contains({folders.name},[wave{w} '_91' session{s}])... 
                | contains({folders.name},[wave{w} '_92' session{s}]) ...
                | contains({folders.name},[wave{w} '_93' session{s}])) ...
                & contains({folders.name},'_normal');   
                infolder = find(results); % So get all blocks you find for that session and condition
                if isempty(infolder)
                    error('no coincidences');
                end            

                for l= 1:length(infolder) % for each block found for that session
                    line = infolder(l);
                    % Check if that block is in the block list, if not, go to next block
                    init_block_name = regexp(folders(line).name, '_b', 'once');
                    fol_nam = contains(block, folders(line).name(init_block_name+1:init_block_name+2));
                    in_fol_nam = find(fol_nam); % index of position in logical
                    if isempty(in_fol_nam); continue; end % the block in the folder name is not present in the list, so go to next
                    % End checking if block exist in list
                    sub_folders = dir([root_dir_bs '/data/' participant{p} '/' folders(line).name '/']);
                    sub_results =  contains({sub_folders.name},[modality_data{mode} '_average_sensor']);
                    sub_infolder = find(sub_results);
                    if isempty(sub_infolder) % if average is not present despite the folder existing
                        continue
                    end
                    file_names_Slowest{list_order_Slowest} = [participant{p} '/' folders(line).name '/' sub_folders(sub_infolder).name]; % sub_results (using logicals)
                    list_order_Slowest = list_order_Slowest + 1;
                end        
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            end

            file_names_Quietest = file_names_Quietest(~cellfun('isempty', file_names_Quietest'));
            file_names_Medium_dB = file_names_Medium_dB(~cellfun('isempty', file_names_Medium_dB'));
            file_names_Loudest = file_names_Loudest(~cellfun('isempty', file_names_Loudest'));
            file_names_Fastest = file_names_Fastest(~cellfun('isempty', file_names_Fastest'));
            file_names_Fast = file_names_Fast(~cellfun('isempty', file_names_Fast'));
            file_names_Medium_ISI = file_names_Medium_ISI(~cellfun('isempty', file_names_Medium_ISI'));
            file_names_Slow = file_names_Slow(~cellfun('isempty', file_names_Slow'));
            file_names_Slowest = file_names_Slowest(~cellfun('isempty', file_names_Slowest'));

            for ex = 1:length(Exp_cond)
                eval(['sFiles = file_names_' Exp_cond{ex} ';'])
                
                if isempty(sFiles)
                error('for whatever reason, sFiles is empty, my friend...')
                end
                
                disp(' ');      
                disp('-------------------------');  
                disp(['Averaging dB and ISI across blocks from all sessions for ' participant{p} '_' modality_data{mode} '_' wave{w} '_' Exp_cond{ex}]);
                disp(datetime)
                disp(' ');
                
                % If stated, find and delete previous normal averages
                if delete_previous_file == 1
                    % check if there is already a normal average
                    folders_delete = dir([root_dir_bs '/data/' participant{p} '/@intra/']);
                    % data_LLR_average_200727_2009_Normal_both_mod_LLR_Fastest.mat
                    results_delete =  endsWith({folders_delete.name}, ['Normal_' modality_data{mode} '_' wave{w} '_' Exp_cond{ex} '.mat']); 
                    infolder_delete = find(results_delete);
                    if ~isempty(infolder_delete) % file exists, therefore delete it
                       delete([root_dir_bs '/data/' participant{p} '/@intra/' folders_delete(infolder_delete).name]);
                    end
                end

                % Process: Weighted Average: Everything
                sFiles = bst_process('CallProcess', 'process_average', sFiles, [], ...
                    'avgtype',       1, ...  % Everything
                    'avg_func',      1, ...  % Arithmetic average:  mean(x)
                    'weighted',      1, ...
                    'keepevents',    0);

                % Process: Add tag
                sFiles = bst_process('CallProcess', 'process_add_tag', sFiles, [], ...
                    'tag',           ['Normal_' modality_data{mode} '_' wave{w} '_' Exp_cond{ex}], ...
                    'output',        2);  % Add to file name (1 to add a tag)

                % Process: Set name
                sFiles = bst_process('CallProcess', 'process_set_comment', sFiles, [], ...
                    'tag',           ['Normal_' modality_data{mode} '_' wave{w} '_' Exp_cond{ex}], ...
                    'isindex',       1);   
            end
        end
    end
end

% reset modality_data variable to its original
modality_data = {'EEG','MEG','both_mod'};

clearvars('-except', initialVars{:});
disp 'DONE WITH AVERAGING WAVEFORMS ACROSS dB AND ISI CONDITIONS (all blocks and sessions)!!!'
disp(datetime)
toc

%% Collapse SENSOR DATA across ALL CONDITIONS (to have representation of LLR)
tic
disp(' ');      
disp('-------------------------');  
disp('AVERAGING WAVEFORMS OF ALL CONDITIONS'); 
disp(datetime)
disp('-------------------------');
disp(' ');

if sensor_analysis == 1
    modality_data = {'EEG'};
elseif sensor_analysis == 2
    modality_data = {'MEG'};
elseif sensor_analysis == 3
    modality_data = {'both_mod'};
end

for p = 1:length(participant)

    % reload subject normally
    prot_subs = bst_get('ProtocolSubjects');
    current_sub = find(strcmp({prot_subs.Subject.Name}, participant{p}));
    db_reload_conditions(current_sub);
    
    % Define sessions within this participant
    pos_par = find(strcmp(session_block_array(:,1), participant{p}(1:4)));
    session = session_block_array{pos_par,2}(1,:);
    
    for mode = 1:length(modality_data)
        for w = 1:length(wave)   

        disp(' ');      
        disp('-------------------------');  
        disp(['Averaging waveforms of ALL conditions for ' participant{p} '_' modality_data{mode} '_' wave{w}]);
        disp(datetime)
        disp(' ');      

        file_names = {}; % create one for each p and w, so that it includes all sessions and blocks
        list_order = 1;

            for s = 1:length(session)
                
                % Define blocks within this session
                pos_ses = find(strcmp(session_block_array{pos_par,2}(1,:), session{s}));
                block = session_block_array{pos_par,2}{2,pos_ses}(1,:);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ALL %%%%%%%%%%%%%%%%%%%%%%%%
                folders = dir([root_dir_bs '/data/' participant{p} '/']);
                results =  (contains({folders.name},[wave{w} '_11' session{s}])... 
                | contains({folders.name},[wave{w} '_12' session{s}]) ...
                | contains({folders.name},[wave{w} '_13' session{s}]) ...
                | contains({folders.name},[wave{w} '_31' session{s}]) ...
                | contains({folders.name},[wave{w} '_32' session{s}]) ...
                | contains({folders.name},[wave{w} '_33' session{s}]) ...
                | contains({folders.name},[wave{w} '_51' session{s}]) ...
                | contains({folders.name},[wave{w} '_52' session{s}]) ...
                | contains({folders.name},[wave{w} '_53' session{s}]) ...
                | contains({folders.name},[wave{w} '_71' session{s}]) ...
                | contains({folders.name},[wave{w} '_72' session{s}]) ...
                | contains({folders.name},[wave{w} '_73' session{s}]) ...
                | contains({folders.name},[wave{w} '_91' session{s}]) ...
                | contains({folders.name},[wave{w} '_92' session{s}]) ...
                | contains({folders.name},[wave{w} '_93' session{s}])) ...
                & contains({folders.name},'_normal');   
                infolder = find(results); % So get all blocks you find for that session and condition
                if isempty(infolder)
                    continue
                    % error('no coincidences');
                end            

                for l= 1:length(infolder) % for each block found for that session
                    line = infolder(l);
                    % Check if that block is in the block list, if not, go to next block
                    init_block_name = regexp(folders(line).name, '_b', 'once');
                    fol_nam = contains(block, folders(line).name(init_block_name+1:init_block_name+2));
                    in_fol_nam = find(fol_nam); % index of position in logical
                    if isempty(in_fol_nam); continue; end % the block in the folder name is not present in the list, so go to next
                    % End checking if block exist in list
                    sub_folders = dir([root_dir_bs '/data/' participant{p} '/' folders(line).name '/']);
                    sub_results =  contains({sub_folders.name},[modality_data{mode} '_average_sensor']);
                    sub_infolder = find(sub_results);
                    if isempty(sub_infolder) % if average is not present despite the folder existing
                        continue
                    end
                    file_names{list_order} = [participant{p} '/' folders(line).name '/' sub_folders(sub_infolder).name]; % sub_results (using logicals)
                    list_order = list_order + 1;
                end        
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end

            sFiles = file_names;
            
            if isempty(sFiles)
            error('for whatever reason, sFiles is empty, my friend...')
            end
            
            % If stated, find and delete previous ALL sensor average
            if delete_previous_file == 1
                % check if there is already an ALL average sensor
                folders_delete = dir([root_dir_bs '/data/' participant{p} '/@intra/']);
                % data_LLR_11B_average_200727_1912_Normal_EEG_LLR_ALL.mat
                results_delete =  endsWith({folders_delete.name}, ['Normal_' modality_data{mode} '_' wave{w} '_ALL.mat']); 
                infolder_delete = find(results_delete);
                if ~isempty(infolder_delete) % file exists, therefore delete it
                   delete([root_dir_bs '/data/' participant{p} '/@intra/' folders_delete(infolder_delete).name]);
                end
            end

            % Process: Weighted Average: Everything
            sFiles = bst_process('CallProcess', 'process_average', sFiles, [], ...
                'avgtype',       1, ...  % Everything
                'avg_func',      1, ...  % Arithmetic average:  mean(x)
                'weighted',      1, ...
                'keepevents',    0);    

            % Process: Add tag
            sFiles = bst_process('CallProcess', 'process_add_tag', sFiles, [], ...
                'tag',           ['Normal_' modality_data{mode} '_' wave{w} '_ALL'], ...
                'output',        2);  % Add to file name (1 to add a tag)

            % Process: Set name
            sFiles = bst_process('CallProcess', 'process_set_comment', sFiles, [], ...
                'tag',           ['Normal_' modality_data{mode} '_' wave{w} '_ALL'], ...
                'isindex',       1);   
        end
    end
end

% reset modality_data variable to its original
modality_data = {'EEG','MEG','both_mod'};

clearvars('-except', initialVars{:});
disp 'DONE WITH AVERAGING WAVEFORMS OF ALL CONDITIONS!!!'
disp(datetime)
toc

%% Get amplitude values of SENSOR DATA out of Brainstorm (Normal epochs)

% Will need to add another one to get the GAVRs out, and maybe putting them
% in a structure that I can use for plotting (watch out with converting
% scales, which I already did in my custom GAVR function).

% Chose whether you want the sensor data extracted to be the best of each
% (MEG and EEG amplitude thresholds separated), or the combined survivor
% sweeps (by default is the best of each)
extracted_sensor_data = 1; % 1 = separated amplitude threshold 2 = combined
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sensor_analysis == 1
    modality_data = {'EEG'};
    size_mode_iteration = 1;
elseif sensor_analysis == 2
    modality_data = {'MEG'};
    size_mode_iteration = 1;
elseif sensor_analysis == 3
    modality_data = {'both_mod'};
    size_mode_iteration = 1;
elseif sensor_analysis == 4 % only then you can chosee between getting amplitudes from combined or separated
    if extracted_sensor_data == 1
        size_mode_iteration = [1:2]; % thus covering only EEG and MEG
    elseif extracted_sensor_data == 2
        size_mode_iteration = 3;
    end
end

tic
disp(' ');      
disp('-------------------------');  
disp('GETTING SENSOR DATA OUT OF BRAINSTORM (NORMAL EPOCHS)');
disp(datetime)
disp('-------------------------');
disp(' ');

% Data from every session (for quality control, invididual responses only)
for p = 1:length(participant)
    
    % reload subject normally
    prot_subs = bst_get('ProtocolSubjects');
    current_sub = find(strcmp({prot_subs.Subject.Name}, participant{p}));
    db_reload_conditions(current_sub);
    
    % Define sessions within this participant
    pos_par = find(strcmp(session_block_array(:,1), participant{p}(1:4)));
    session = session_block_array{pos_par,2}(1,:);
    
    for s = 1:length(session)
        for mode = size_mode_iteration
            for w = 1:length(wave)
                for c = 1:length(condition)  
                folders = dir([root_dir_bs '/data/' participant{p} '/@intra/']);
                results =  endsWith({folders.name},['Normal_' modality_data{mode} '_' wave{w} '_' condition{c} session{s} '.mat']);
                infolder = find(results);
                if isempty(infolder) % If there are no coincidences, continue to next block/session
                    continue
                end
                file_name = [root_dir_bs '/data/' participant{p} '/@intra/' folders(infolder).name];            
                % This way so that we can preserve the structure for the plotting script
                eval(['load ' file_name])
                    if strcmp(modality_data{mode}, 'EEG')
                        F_EEG = F(315:376,:);
                        clear F
                        F = F_EEG;
                        if ~exist([root_dir '/Averages/' participant{p} '/' session{s}], 'dir')
                            mkdir([root_dir '/Averages/'], [participant{p} '/' session{s}]);
                        end
                        save([root_dir '/Averages/' participant{p} '/' session{s} '/EEG_' wave{w} '_' condition{c} '.mat'], 'F');
                        clear F
                    elseif strcmp(modality_data{mode}, 'MEG')
                        F_MEG = F(1:306,:);
                        clear F
                        if ~exist([root_dir '/Averages/' participant{p} '/' session{s}], 'dir')
                            mkdir([root_dir '/Averages/'], [participant{p} '/' session{s}]);
                        end
                        F = F_MEG;
                        save([root_dir '/Averages/' participant{p} '/' session{s} '/MEG_' wave{w} '_' condition{c} '.mat'], 'F');
                        clear F
                    elseif strcmp(modality_data{mode}, 'both_mod')
                        F_EEG = F(315:376,:);
                        F_MEG = F(1:306,:);
                        clear F
                        F = F_EEG;
                        if ~exist([root_dir '/Averages/' participant{p} '/' session{s}], 'dir')
                            mkdir([root_dir '/Averages/'], [participant{p} '/' session{s}]);
                        end
                        save([root_dir '/Averages/' participant{p} '/' session{s} '/EEG_' wave{w} '_' condition{c} '.mat'], 'F');
                        clear F
                        F = F_MEG;
                        save([root_dir '/Averages/' participant{p} '/' session{s} '/MEG_' wave{w} '_' condition{c} '.mat'], 'F');
                    end
                end
            end
        end
    end
end

% Data from every session (dB and ISI for plot purposes)
for p = 1:length(participant)
    
    % reload subject normally
    prot_subs = bst_get('ProtocolSubjects');
    current_sub = find(strcmp({prot_subs.Subject.Name}, participant{p}));
    db_reload_conditions(current_sub);
    
    % Define sessions within this participant
    pos_par = find(strcmp(session_block_array(:,1), participant{p}(1:4)));
    session = session_block_array{pos_par,2}(1,:);
    
    for s = 1:length(session)
        for mode = size_mode_iteration
            for w = 1:length(wave)
                for ex = 1:length(Exp_cond)  
                    folders = dir([root_dir_bs '/data/' participant{p} '/@intra/']);
                    results =  endsWith({folders.name},['Normal_' modality_data{mode} '_' wave{w} '_' Exp_cond{ex} '_' session{s} '.mat']);
                    infolder = find(results);
                    if isempty(infolder) % If there are no coincidences, continue to next block/session
                        continue
                    end
                    file_name = [root_dir_bs '/data/' participant{p} '/@intra/' folders(infolder).name];            
                    % This way so that we can preserve the structure for the plotting script
                    eval(['load ' file_name])
                    if strcmp(modality_data{mode}, 'EEG')
                        F_EEG = F(315:376,:);
                        clear F
                        F = F_EEG;
                        if ~exist([root_dir '/Averages/' participant{p} '/' session{s}], 'dir')
                            mkdir([root_dir '/Averages/'], [participant{p} '/' session{s}]);
                        end
                        save([root_dir '/Averages/' participant{p} '/' session{s} '/EEG_' wave{w} '_' Exp_cond{ex} '.mat'], 'F');
                        clear F
                    elseif strcmp(modality_data{mode}, 'MEG')
                        F_MEG = F(1:306,:);
                        clear F
                        if ~exist([root_dir '/Averages/' participant{p} '/' session{s}], 'dir')
                            mkdir([root_dir '/Averages/'], [participant{p} '/' session{s}]);
                        end
                        F = F_MEG;
                        save([root_dir '/Averages/' participant{p} '/' session{s} '/MEG_' wave{w} '_' Exp_cond{ex} '.mat'], 'F');
                        clear F
                    elseif strcmp(modality_data{mode}, 'both_mod')
                        F_EEG = F(315:376,:);
                        F_MEG = F(1:306,:);
                        clear F
                        F = F_EEG;
                        if ~exist([root_dir '/Averages/' participant{p} '/' session{s}], 'dir')
                            mkdir([root_dir '/Averages/'], [participant{p} '/' session{s}]);
                        end
                        save([root_dir '/Averages/' participant{p} '/' session{s} '/EEG_' wave{w} '_' Exp_cond{ex} '.mat'], 'F');
                        clear F
                        F = F_MEG;
                        save([root_dir '/Averages/' participant{p} '/' session{s} '/MEG_' wave{w} '_' Exp_cond{ex} '.mat'], 'F');
                    end
                end
            end
        end
    end
end

% Individual responses of each participant and condition
for p = 1:length(participant)
    if ~exist([root_dir '/Averages/' participant{p}], 'dir')
       mkdir([root_dir '/Averages/'], [participant{p}]);
    end
    for mode = size_mode_iteration
        for w = 1:length(wave)
            for c = 1:length(condition)  
            folders = dir([root_dir_bs '/data/' participant{p} '/@intra/']);
            results =  contains({folders.name},['Normal_' modality_data{mode} '_' wave{w} '_' condition{c} '.mat']);
            infolder = find(results);
            file_name = [root_dir_bs '/data/' participant{p} '/@intra/' folders(infolder).name];            
            % This way so that we can preserve the structure for the plotting script
            eval(['load ' file_name])
                if strcmp(modality_data{mode}, 'EEG')
                    F_EEG = F(315:376,:);
                    clear F
                    F = F_EEG;
                    save([root_dir '/Averages/' participant{p} '/EEG_' wave{w} '_' condition{c} '.mat'], 'F');
                    clear F
                elseif strcmp(modality_data{mode}, 'MEG')
                    F_MEG = F(1:306,:);
                    clear F
                    F = F_MEG;
                    save([root_dir '/Averages/' participant{p} '/MEG_' wave{w} '_' condition{c} '.mat'], 'F');
                    clear F
                elseif strcmp(modality_data{mode}, 'both_mod')
                    F_EEG = F(315:376,:);
                    F_MEG = F(1:306,:);
                    clear F
                    F = F_EEG;
                    save([root_dir '/Averages/' participant{p} '/EEG_' wave{w} '_' condition{c} '.mat'], 'F');
                    clear F
                    F = F_MEG;
                    save([root_dir '/Averages/' participant{p} '/MEG_' wave{w} '_' condition{c} '.mat'], 'F');
                end          
            end
        end
    end
end

% Responses of each participant averaged by dB and ISI
for p = 1:length(participant)
    if ~exist([root_dir '/Averages/' participant{p}], 'dir')
       mkdir([root_dir '/Averages/'], [participant{p}]);
    end
    for mode = size_mode_iteration 
        for w = 1:length(wave)
            for ex = 1:length(Exp_cond)
            folders = dir([root_dir_bs '/data/' participant{p} '/@intra/']);
            results =  contains({folders.name},['Normal_' modality_data{mode} '_' wave{w} '_' Exp_cond{ex} '.mat']);
            infolder = find(results);
            file_name = [root_dir_bs '/data/' participant{p} '/@intra/' folders(infolder).name];            
            % This way so that we can preserve the structure for the plotting script
            eval(['load ' file_name]) 
                if strcmp(modality_data{mode}, 'EEG')
                    F_EEG = F(315:376,:);
                    clear F
                    F = F_EEG;
                    save([root_dir '/Averages/' participant{p} '/EEG_' wave{w} '_' Exp_cond{ex} '.mat'], 'F');
                    clear F
                elseif strcmp(modality_data{mode}, 'MEG')
                    F_MEG = F(1:306,:);
                    clear F
                    F = F_MEG;
                    save([root_dir '/Averages/' participant{p} '/MEG_' wave{w} '_' Exp_cond{ex} '.mat'], 'F');
                    clear F
                elseif strcmp(modality_data{mode}, 'both_mod')
                    F_EEG = F(315:376,:);
                    F_MEG = F(1:306,:);
                    clear F
                    F = F_EEG;
                    save([root_dir '/Averages/' participant{p} '/EEG_' wave{w} '_' Exp_cond{ex} '.mat'], 'F');
                    clear F
                    F = F_MEG;
                    save([root_dir '/Averages/' participant{p} '/MEG_' wave{w} '_' Exp_cond{ex} '.mat'], 'F');
                end
            end
        end
    end
end

% Responses of each participant averaged across everything (ALL)
for p = 1:length(participant)
    if ~exist([root_dir '/Averages/' participant{p}], 'dir')
       mkdir([root_dir '/Averages/'], [participant{p}]);
    end
    for mode = size_mode_iteration
        for w = 1:length(wave)
        folders = dir([root_dir_bs '/data/' participant{p} '/@intra/']);
        results =  contains({folders.name},['Normal_' modality_data{mode} '_' wave{w} '_ALL.mat']);
        infolder = find(results);
        file_name = [root_dir_bs '/data/' participant{p} '/@intra/' folders(infolder).name];            
        % This way so that we can preserve the structure for the plotting script
        eval(['load ' file_name])
            if strcmp(modality_data{mode}, 'EEG')
                F_EEG = F(315:376,:);
                clear F
                F = F_EEG;
                save([root_dir '/Averages/' participant{p} '/EEG_' wave{w} '_ALL.mat'], 'F');
                clear F
            elseif strcmp(modality_data{mode}, 'MEG')
                F_MEG = F(1:306,:);
                clear F
                F = F_MEG;
                save([root_dir '/Averages/' participant{p} '/MEG_' wave{w} '_ALL.mat'], 'F');
                clear F
            elseif strcmp(modality_data{mode}, 'both_mod')
                F_EEG = F(315:376,:);
                F_MEG = F(1:306,:);
                clear F
                F = F_EEG;
                save([root_dir '/Averages/' participant{p} '/EEG_' wave{w} '_ALL.mat'], 'F');
                clear F
                F = F_MEG;
                save([root_dir '/Averages/' participant{p} '/MEG_' wave{w} '_ALL.mat'], 'F');
            end        
        end
    end
end

% reset modality_data variable to its original
modality_data = {'EEG','MEG','both_mod'};

clearvars('-except', initialVars{:});
disp 'DONE GETTING SENSOR DATA OUT OF BRAINSTORM (NORMAL EPOCHS)!!!'
disp(datetime)
toc
     
%% Make epochs for noise covariance analysis (Part 2)

if covar_epochs == 1

tic
disp(' ');      
disp('-------------------------');  
disp('OBTAINING EPOCHS FOR COVARIANCE ANALYSES (Part 2)');  
disp(datetime)
disp('-------------------------');
disp(' ');

for p = 1:length(participant)
    
    disp(' ');      
    disp('-------------------------');
    disp(['loading participant ' participant{p}]);
    disp(datetime)
    disp(' '); 
    
    prot_subs = bst_get('ProtocolSubjects');
    current_sub = find(strcmp({prot_subs.Subject.Name}, participant{p}));
    db_reload_conditions(current_sub);
    
    % Define sessions within this participant
    pos_par = find(strcmp(session_block_array(:,1), participant{p}(1:4)));
    session = session_block_array{pos_par,2}(1,:);
    
    for s = 1:length(session)
        
        % Define blocks within this session
        pos_ses = find(strcmp(session_block_array{pos_par,2}(1,:), session{s}));
        block = session_block_array{pos_par,2}{2,pos_ses}(1,:);
        
        for b = 1:length(block)
            for c =13:length(condition) % ONLY 91, 92 AND 93
                for w = 1:length(wave)
                    
                    folders = dir([root_dir_bs '/data/' participant{p} '/']);
                    results = contains({folders.name},[wave{w} '_' condition{c} session{s} '_' block{b} '_covariance']);
                    infolder = find(results);
                    
                    if isempty(infolder)
                        continue
                    end
                    if size(infolder,2)> 1 % in case of more than one coincidence, error
                        error('More than one coincidence');
                    end
                        
                    file_names = {};
                    dir_sweeps = dir([root_dir_bs '/data/' participant{p} '/' folders(infolder).name]);
                    sweeps_norm = contains({dir_sweeps.name},'_trial');
                    sweeps_list = find(sweeps_norm);
                    for j= 1:length(sweeps_list)
                        position = sweeps_list(j);
                        file_names{j} = [participant{p} '/' folders(infolder).name '/' dir_sweeps(position).name];
                    end
                    sFiles = file_names; % all sweeps from that block and condition contained here now
                    
                    if isempty(sFiles)
                    error('for whatever reason, sFiles is empty, my friend...')
                    end
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Amplitude threshold normal epochs %%%%%%%%%%%%%%%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    if (sensor_analysis == 1) || (sensor_analysis == 4) % EEG only
                    disp(' ');      
                    disp('-------------------------');  
                    disp(['Cleaning covariance epochs for ' participant{p} '_' wave{w} '_' condition{c} session{s} '_' block{b} '(EEG)']);
                    disp(datetime)
                    disp(' '); 

                    % Process: Detect bad trials: Absolute threshold EEG
                    sFiles_EEG = bst_process('CallProcess', 'process_CNRL_detectbad', sFiles, [], ...
                        'timewindow', [], ...
                        'meggrad',    [0, 0], ...
                        'megmag',     [0, 0], ...
                        'eeg',        reject_EEG_absolute, ...
                        'ieeg',       [0, 0], ...
                        'eog',        [0, 0], ...
                        'ecg',        [0, 0], ...
                        'rejectmode', 2);  % Reject the entire trial

                    % save and reset bad trials
                    trials_file = [root_dir_bs '/data/' participant{p} '/' folders(infolder).name '/brainstormstudy.mat'];
                    eval(['load ' trials_file]); % load them
                    % save current ones in events folder
                    if ~exist([root_dir '/Events/BadTrials_Cov/' ProtocolName(7:9) '/EEG/' participant{p} '/' folders(infolder).name], 'dir')
                        mkdir([root_dir '/Events/BadTrials_Cov/'], [ProtocolName(7:9) '/EEG/' participant{p} '/' folders(infolder).name]);
                    end
                    save([root_dir '/Events/BadTrials_Cov/' ProtocolName(7:9) '/EEG/' participant{p} '/' folders(infolder).name '/BadTrials.mat'],'BadTrials'); % save bad trials
                    % reset them and save back to original file
                    BadTrials = cell([], 1); % exact structure it has when empty
                    save(trials_file,'BadTrials', 'DateOfStudy', 'Name');
                    end

                    % BadTrials, DateOfStudy, Name

                    if (sensor_analysis == 2) || (sensor_analysis == 4) % MEG only

                    disp(' ');      
                    disp('-------------------------');  
                    disp(['Cleaning covariance epochs for ' participant{p} '_' wave{w} '_' condition{c} session{s} '_' block{b} '(MEG)']);
                    disp(datetime)
                    disp(' '); 

                    % Process: Detect bad trials: Peak-to-peak  MEG ONLY
                    sFiles_MEG = bst_process('CallProcess', 'process_CNRL_detectbad',sFiles, [], ...
                        'timewindow', [], ...
                        'meggrad',    reject_MEG_GRAD_absolute, ...
                        'megmag',     reject_MEG_MAG_absolute, ...
                        'eeg',        [0, 0], ...
                        'ieeg',       [0, 0], ...
                        'eog',        [0, 0], ...
                        'ecg',        [0, 0], ...
                        'rejectmode', 2);  % Reject the entire trial

                    % save and reset bad trials
                    trials_file = [root_dir_bs '/data/' participant{p} '/' folders(infolder).name '/brainstormstudy.mat'];
                    eval(['load ' trials_file]); % load them
                    % save current ones in events folder
                    if ~exist([root_dir '/Events/BadTrials_Cov/' ProtocolName(7:9) '/MEG/' participant{p} '/' folders(infolder).name], 'dir')
                        mkdir([root_dir '/Events/BadTrials_Cov/'], [ProtocolName(7:9) '/MEG/' participant{p} '/' folders(infolder).name]);
                    end
                    save([root_dir '/Events/BadTrials_Cov/' ProtocolName(7:9) '/MEG/' participant{p} '/' folders(infolder).name '/BadTrials.mat'],'BadTrials'); % save bad trials
                    % reset them and save back to original file
                    BadTrials = cell([], 1); % exact structure it has when empty
                    save(trials_file,'BadTrials', 'DateOfStudy', 'Name');
                    end

                    if (sensor_analysis == 3) || (sensor_analysis == 4) % Combined EEG and MEG

                    disp(' ');      
                    disp('-------------------------');  
                    disp(['Cleaning covariance epochs for ' participant{p} '_' wave{w} '_' condition{c} session{s} '_' block{b} '(both_mod)']);
                    disp(datetime)
                    disp(' '); 

                    % Process: Detect bad trials: Peak-to-peak  BOTH EEG AND MEG
                    sFiles_both_mod = bst_process('CallProcess', 'process_CNRL_detectbad',sFiles, [], ...
                        'timewindow', [], ...
                        'meggrad',    reject_MEG_GRAD_absolute, ...
                        'megmag',     reject_MEG_MAG_absolute, ...
                        'eeg',        reject_EEG_absolute, ...
                        'ieeg',       [0, 0], ...
                        'eog',        [0, 0], ...
                        'ecg',        [0, 0], ...
                        'rejectmode', 2);  % Reject the entire trial

                    % save and reset bad trials
                    trials_file = [root_dir_bs '/data/' participant{p} '/' folders(infolder).name '/brainstormstudy.mat'];
                    eval(['load ' trials_file]); % load them
                    % save current ones in events folder
                    if ~exist([root_dir '/Events/BadTrials_Cov/' ProtocolName(7:9) '/both_mod/' participant{p} '/' folders(infolder).name], 'dir')
                        mkdir([root_dir '/Events/BadTrials_Cov/'], [ProtocolName(7:9) '/both_mod/' participant{p} '/' folders(infolder).name]);
                    end
                    save([root_dir '/Events/BadTrials_Cov/' ProtocolName(7:9) '/both_mod/' participant{p} '/' folders(infolder).name '/BadTrials.mat'],'BadTrials'); % save bad trials
                    % reset them and save back to original file
                    BadTrials = cell([], 1); % exact structure it has when empty
                    save(trials_file,'BadTrials', 'DateOfStudy', 'Name');
                    end
                end
            end
        end
    end
end

clearvars('-except', initialVars{:});
disp 'DONE WITH OBTAINING EPOCHS FOR COVARIANCE ANALYSES (Part 2)!!!'
disp(datetime)
toc

end

%% Count trials for Covariance epochs

if covar_epochs == 1

tic
disp(' ');      
disp('-------------------------');  
disp('COUNTING TRIALS FOR COVARIANCE EPOCHS (LLR)'); 
disp(datetime)
disp('-------------------------');     
disp(' ');

if sensor_analysis == 1
    colNames = {'PRE_LLR','EEG_Post'};
    modality_data = {'EEG'};
elseif sensor_analysis == 2
    colNames = {'PRE_LLR','MEG_Post'};
    modality_data = {'MEG'};
elseif sensor_analysis == 3
    colNames = {'PRE_LLR','Combined_Post'};
    modality_data = {'both_mod'};
elseif sensor_analysis == 4
    colNames = {'PRE_LLR','EEG_Post','MEG_Post','Combined_Post'};
end

cov_condition = {'91' '92' '93'};

% For covariance 
% A detailed one for covariance trials with 1(91,92 + 93)
for p = 1:length(participant)
    
    % Define sessions within this participant
    pos_par = find(strcmp(session_block_array(:,1), participant{p}(1:4)));
    session = session_block_array{pos_par,2}(1,:);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Cov_trials = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for w = 1:length(modality_data) % Just because of inheritance from previous script, we keep 'w' instead of 'mode'
        iter_modality = 0; % first time iteration for position purposes
        current_size_matrix = []; % just for size purposes
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        sum_orig = zeros(15,1);
        sum_final = zeros(15,1);
        rowNames_detailed = {}; % create a new for every subject (bec of diff num of sessions) 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for s = 1:length(session)
            
            % Define blocks within this session
            pos_ses = find(strcmp(session_block_array{pos_par,2}(1,:), session{s}));
            block = session_block_array{pos_par,2}{2,pos_ses}(1,:);
            
            for b = 1:length(block) 
                folders = dir([root_dir_bs '/data/' participant{p} '/']);
                results = contains({folders.name},['LLR_91' session{s} '_' block{b} '_covariance'])...
                    | contains({folders.name},['LLR_92' session{s} '_' block{b} '_covariance'])...
                    | contains({folders.name},['LLR_93' session{s} '_' block{b} '_covariance']);
                infolder = find(results);
                if isempty(infolder)
                    continue
                end
                
                disp(' ');      
                disp('-------------------------');
                disp(['Counting covariance sweeps for ' participant{p} 'LLR_' modality_data{w} '_' session{s} '_' block{b}]);
                disp(' '); 
                
                for l= 1:length(infolder) % For each of the folders (91, 92 and 93 sweeps togheter)
                    line = infolder(l);
                    dir_sweeps = dir([root_dir_bs '/data/' participant{p} '/' folders(line).name]);
                    sweeps_cov = contains({dir_sweeps.name},'_trial');
                    sweeps_count = find(sweeps_cov);

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % COUNT REJECTED TRIALS
                    try 
                        load([root_dir '/Events/BadTrials_Cov/' ProtocolName(7:9) '/' modality_data{w} '/' participant{p} '/' folders(line).name '/BadTrials.mat']); 
                        if ~isempty(BadTrials)
                            count_rej = length(BadTrials);
                        else
                            count_rej = 0;
                        end
                    catch
                        % If no BadTrials variable is ready to load, countas 0
                        count_rej = 0;
                    end

                    % COUNT ORIGINAL TRIALS
                    count_orig = length(sweeps_count);
                    % COUNT NON REJECTED TRIALS
                    count_final = length(sweeps_count) - count_rej;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

                    if length(Cov_trials) == 0  %#ok<ISMT>
                        pos = 1;
                    else % determine last position and put one after
                        size_matrix = size(Cov_trials,1);
                        pos = size_matrix +1; 
                    end
                    
                    if iter_modality == 0 % first iteration in this modality
                        pos = 1;
                        iter_modality = 1; % this won't change until we go to next modality
                        % length(Normal_trials_detailed) == 0  % how it was
                        % before %#ok<ISMT> % DO NOT CHANGE BY ISEMPTY
                    else % determine last position and put one after
                        % size_matrix = size(Normal_trials_detailed,1);
                        % pos = size_matrix +1; 
                        pos = size(current_size_matrix,1) + 1;
                    end
                    
                    % pos = ((s-1)*45) + ((b-1)*15) + l; % index position in table
                    rowNames_detailed{1,pos} = [session{s} '_' block{b} '_' cov_condition{l}];
                    Cov_trials(pos,1) = count_orig;
                    Cov_trials(pos,w+1) = count_final;
                    current_size_matrix(pos,1) = 1; % just to keep track of size
                end
            end
        end
    end
    % Delete empty rows in detailed axis names and data
    % since we have fixed blocks and sessions, so if there is an empty spot
    % is worth noticing)
    % rowNames_detailed = rowNames_detailed(any(~cellfun('isempty',rowNames_detailed), 1));
    % Cov_trials = Cov_trials(any(Cov_trials,2),:);
    % Create tables from data matrix and save
    eval(['LLR_' participant{p} '_Cov_trials = array2table(Cov_trials,''RowNames'',rowNames_detailed,''VariableNames'',colNames);'])
    save([root_dir '/Events/LLR_Counts_' participant{p} '_Covariance.mat'],['LLR_' participant{p} '_Cov_trials']);
end

% reset modality_data variable to its original
modality_data = {'EEG','MEG','both_mod'};

clearvars('-except', initialVars{:});
disp 'DONE WITH COUNTING TRIALS FOR COVARIANCE EPOCHS (LLR)!!!'
disp(datetime)
toc

end

%% GAVR sensor level (ALWAYS BE SURE TO HAVEN'T SHORTENED THE PARTICIPANT VARIABLE)

tic
disp(' ');      
disp('-------------------------');  
disp('GAVR SENSOR DATA (LLR)'); 
disp(datetime)
disp('-------------------------');
disp(' ');

% To correct if we used a loop through sections
participant = participant_general_list;

% Reload all subjects before anything
for p = 1:length(participant)
    
    disp(' ');      
    disp('-------------------------');
    disp(['loading participant ' participant{p}]);
    disp(datetime)
    disp(' '); 
    
    prot_subs = bst_get('ProtocolSubjects');
    current_sub = find(strcmp({prot_subs.Subject.Name}, participant{p}));
    db_reload_conditions(current_sub);
end


if sensor_analysis == 1
    modality_data = {'EEG'};
elseif sensor_analysis == 2
    modality_data = {'MEG'};
elseif sensor_analysis == 3
    modality_data = {'both_mod'};
end

% Reload Group analysis folder too (destiny files)
prot_subs = bst_get('ProtocolSubjects');
current_sub = find(strcmp({prot_subs.Subject.Name}, 'Group_analysis'));
db_reload_conditions(current_sub);

% Average source data individual condition (11, 12, etc) across subjects
for mode = 1:length(modality_data)
for w = 1:length(wave)
    for c = 1:length(condition)
         file_names = {};
         for p = 1:length(participant)
            folders = dir([root_dir_bs '/data/' participant{p} '/@intra/']);
            results =  contains({folders.name},['Normal_' modality_data{mode} '_' wave{w}]) & ...
                endsWith({folders.name},[condition{c} '.mat']);
            infolder = find(results);
            if isempty(infolder) % In case, for instance, MLR files are not there yet
                error(['No ' modality_data{mode} '_' wave{w} '_' condition{c} ' file for participant' participant{p}]);
            elseif size(infolder,2)> 1 % in case of more than one coincidence (should not happen), tell me
               error(['More than one ' modality_data{mode} '_' wave{w} '_' condition{c} ' file for participant' participant{p}]);
            end

            file_names{p} = [participant{p} '/@intra/' folders(infolder).name];
         end
            sFiles = file_names;

            if isempty(sFiles)
                error('for whatever reason, sFiles is empty, my friend...')
            end

            disp(' ');
            disp('-------------------------');  
            disp(['GAVR sensor data for ' modality_data{mode} '_' wave{w} '_' condition{c}]);
            disp(datetime)
            disp(' ');

            % If stated, find and delete any previous GAVR SENSOR data
            if delete_previous_file == 1
                % check if there is already GAVR source in Group analysis folder
                folders_delete = dir([root_dir_bs '/data/Group_analysis/@intra/']);
                % results_average_200520_2029_GAVR_Source_MEG_Regul_MLR_Quietest
                results_delete = contains({folders_delete.name},'GAVR_SENSOR_')...
                & endsWith({folders_delete.name}, [modality_data{mode} '_' wave{w} '_' condition{c} '.mat']); 
                infolder_delete = find(results_delete);
                if ~isempty(infolder_delete) % file exists, therefore delete it
                   delete([root_dir_bs '/data/Group_analysis/@intra/' folders_delete(infolder_delete).name]);
                end
            end

            % Average using default function (because we are using vertices, not channels)
            sFiles = bst_process('CallProcess', 'process_average', sFiles, [], ...
                'avgtype',         1, ...  % Everything
                'avg_func',        1, ...  % Arithmetic average:  mean(x)
                'weighted',        1, ...
                'scalenormalized', 0);

            % USE OPTION TO NORMALIZE HERE???

            % Process: Add tag
            sFiles = bst_process('CallProcess', 'process_add_tag', sFiles, [], ...
                'tag',           ['GAVR_SENSOR_' modality_data{mode} '_' wave{w} '_' condition{c}], ...
                'output',        2);  % Add to file name (1 to add a tag)

            % Process: Set name
            sFiles = bst_process('CallProcess', 'process_set_comment', sFiles, [], ...
                'tag',           ['GAVR_SENSOR_' modality_data{mode} '_' wave{w} '_' condition{c}], ...
                'isindex',       1);
          
    end
end
end

% Average source categories (dB and ISI) across subjects
for mode = 1:length(modality_data)
for w = 1:length(wave)
    for ex = 1:length(Exp_cond)
         file_names = {};
         for p = 1:length(participant)
            folders = dir([root_dir_bs '/data/' participant{p} '/@intra/']);
            results =  contains({folders.name},['Normal_' modality_data{mode} '_' wave{w}]) & ...
                endsWith({folders.name},[Exp_cond{ex} '.mat']);
            infolder = find(results);
            if isempty(infolder) % In case, for instance, MLR files are not there yet
                error(['No ' modality_data{mode} '_' wave{w} '_' Exp_cond{ex} ' file for participant' participant{p}]);
            elseif size(infolder,2)> 1 % in case of more than one coincidence (should not happen), tell me
               error(['More than one ' modality_data{mode} '_' wave{w} '_' Exp_cond{ex} ' file for participant' participant{p}]);
            end

            file_names{p} = [participant{p} '/@intra/' folders(infolder).name];
         end
            sFiles = file_names;

            if isempty(sFiles)
                error('for whatever reason, sFiles is empty, my friend...')
            end

            disp(' ');
            disp('-------------------------');  
            disp(['GAVR sensor data for ' modality_data{mode} '_' wave{w} '_' Exp_cond{ex}]);
            disp(datetime)
            disp(' ');

            % If stated, find and delete any previous GAVR SENSOR data
            if delete_previous_file == 1
                % check if there is already GAVR source in Group analysis folder
                folders_delete = dir([root_dir_bs '/data/Group_analysis/@intra/']);
                % results_average_200520_2029_GAVR_Source_MEG_Regul_MLR_Quietest
                results_delete = contains({folders_delete.name},'GAVR_SENSOR_')...
                & endsWith({folders_delete.name}, [modality_data{mode} '_' wave{w} '_' Exp_cond{ex} '.mat']); 
                infolder_delete = find(results_delete);
                if ~isempty(infolder_delete) % file exists, therefore delete it
                   delete([root_dir_bs '/data/Group_analysis/@intra/' folders_delete(infolder_delete).name]);
                end
            end

            % Average using default function (because we are using vertices, not channels)
            sFiles = bst_process('CallProcess', 'process_average', sFiles, [], ...
                'avgtype',         1, ...  % Everything
                'avg_func',        1, ...  % Arithmetic average:  mean(x)
                'weighted',        1, ...
                'scalenormalized', 0);

            % USE OPTION TO NORMALIZE HERE???

            % Process: Add tag
            sFiles = bst_process('CallProcess', 'process_add_tag', sFiles, [], ...
                'tag',           ['GAVR_SENSOR_' modality_data{mode} '_' wave{w} '_' Exp_cond{ex}], ...
                'output',        2);  % Add to file name (1 to add a tag)

            % Process: Set name
            sFiles = bst_process('CallProcess', 'process_set_comment', sFiles, [], ...
                'tag',           ['GAVR_SENSOR_' modality_data{mode} '_' wave{w} '_' Exp_cond{ex}], ...
                'isindex',       1);
          
    end
end
end

% Average collapsed waveform ("ALL") sources across subjects
for mode = 1:length(modality_data)
for w = 1:length(wave)
     file_names = {};
     for p = 1:length(participant)
        folders = dir([root_dir_bs '/data/' participant{p} '/@intra/']);
        results =  contains({folders.name},['Normal_' modality_data{mode} '_' wave{w}]) & ...
            endsWith({folders.name},'ALL.mat');
        infolder = find(results);
        if isempty(infolder) % In case, for instance, MLR files are not there yet
            error(['No ' modality_data{mode} '_' wave{w} '_ALL file for participant' participant{p}]);
        elseif size(infolder,2)> 1 % in case of more than one coincidence (should not happen), tell me
           error(['More than one ' modality_data{mode} '_' wave{w} '_ALL  file for participant' participant{p}]);
        end

        file_names{p} = [participant{p} '/@intra/' folders(infolder).name];
     end
        sFiles = file_names;

        if isempty(sFiles)
            error('for whatever reason, sFiles is empty, my friend...')
        end

        disp(' ');
        disp('-------------------------');  
        disp(['GAVR sensor data for ' modality_data{mode} '_' wave{w} '_ALL']);
        disp(datetime)
        disp(' ');

        % If stated, find and delete any previous GAVR SENSOR data
        if delete_previous_file == 1
            % check if there is already GAVR source in Group analysis folder
            folders_delete = dir([root_dir_bs '/data/Group_analysis/@intra/']);
            % results_average_200520_2029_GAVR_Source_MEG_Regul_MLR_Quietest
            results_delete = contains({folders_delete.name},'GAVR_SENSOR_')...
            & endsWith({folders_delete.name}, [modality_data{mode} '_' wave{w} '_ALL.mat']); 
            infolder_delete = find(results_delete);
            if ~isempty(infolder_delete) % file exists, therefore delete it
               delete([root_dir_bs '/data/Group_analysis/@intra/' folders_delete(infolder_delete).name]);
            end
        end

        % Average using default function (because we are using vertices, not channels)
        sFiles = bst_process('CallProcess', 'process_average', sFiles, [], ...
            'avgtype',         1, ...  % Everything
            'avg_func',        1, ...  % Arithmetic average:  mean(x)
            'weighted',        1, ...
            'scalenormalized', 0);

        % USE OPTION TO NORMALIZE HERE???

        % Process: Add tag
        sFiles = bst_process('CallProcess', 'process_add_tag', sFiles, [], ...
            'tag',           ['GAVR_SENSOR_' modality_data{mode} '_' wave{w} '_ALL'], ...
            'output',        2);  % Add to file name (1 to add a tag)

        % Process: Set name
        sFiles = bst_process('CallProcess', 'process_set_comment', sFiles, [], ...
            'tag',           ['GAVR_SENSOR_' modality_data{mode} '_' wave{w} '_ALL'], ...
            'isindex',       1);     
end
end

% reset modality_data variable to its original
modality_data = {'EEG','MEG','both_mod'};

clearvars('-except', initialVars{:});
disp 'DONE WITH GAVR SENSOR DATA (LLR)!!!'
disp(datetime)
toc

%% CORREGISTRATION MUST BE DONE BEFORE GOING ANY FURTHER
% Import the anatomy (HCP import list), set the fiducials and corregister MRI with sensors
% Corregistration should be done once per session only, not per block
% The reason is we are corregistering MRI with the HPI coils and EEG, which
% do not change their position within the session
% HOWEVER, A SCRIPT WILL BE NEEDED TO COPY FROM ONE BLOCK TO OTHERS WITHIN
% SESSION... include here

% Use information here if there are problems of fitting electrodes in the
% head shape: https://neuroimage.usc.edu/brainstorm/Tutorials/LabelFreeSurfer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Create a security copy of channel data as it is to 'reset' if needed

% FOR EEG OR BOTH MOD SOURCE SOLUTIONS, WE SHOULD PROJECT ELECTRODES ON THE
% SURFACE OF THE HEAD, WHICH CAN BE DONE EASILY WITH THE GUI (RESPONSE FROM
% BRAINSTORM DEVELOPERS). THAT WON'T AFFECT THE POSITION OF THE MEG
% SENSORS AND IF CORREGISTRATION IS REPEATED THE POSITIONS OF THE
% ELECTRODES WHEN DOING MEG EDIT IS THE SAME. HOWEVER, JUST TO BE SURE,
% BETTER NOT TO PROJECT ELECTRODES TO SURFACE UNTIL ALL CORREGISTRATIONS
% ARE CORRECT OR, IF NOT, PROJECT TO SURFACE AGAIN EVERYTIME YOU CHANGE
% SOMETHING IN CORREGISTRATION.

% JUST TO BE CAUTIONS, WE WILL CREATE A SECURITY COPY
% AND OVERWRITE EVERYTIME IS NECESSARY ALL THE CHANNEL FILES IN THE PROPER
% FOLDERS (THIS FILES CAN BE THEN COPIED OVER TO THE MLR SEG WITH THE OTHER
% SCRIPT THERE)

tic
disp(' ');      
disp('-------------------------');  
disp('CREATING A SECURITY COPY OF CHANNEL FILES');  
disp(datetime)
disp('-------------------------');     
disp(' '); 

for p = 1:length(participant)
      
    % Define sessions within this participant
    pos_par = find(strcmp(session_block_array(:,1), participant{p}(1:4)));
    session = session_block_array{pos_par,2}(1,:);

    for s = 1:length(session)
        
        % Define blocks within this session
        pos_ses = find(strcmp(session_block_array{pos_par,2}(1,:), session{s}));
        block = session_block_array{pos_par,2}{2,pos_ses}(1,:);
        
        for b = 1:length(block)
            for c =1:length(condition)
                for w = 1:length(wave)
                    % If that channel file exist in LLR seg
                    if exist([root_dir_bs '/data/' participant{p} '/' wave{w} '_' condition{c} session{s} '_' block{b} '_normal/channel_vectorview306_acc1.mat'],'file')
                        % Create a folder to save the channel file if needed
                        if ~exist([root_dir '/Channel_backup/' participant{p} '/' condition{c} session{s} '_' block{b} '_normal'], 'dir')
                            mkdir([root_dir '/Channel_backup/'], [participant{p} '/' condition{c} session{s} '_' block{b} '_normal']);
                        end
                        % copy the channel file in that folder
                        current_channel_file = [root_dir_bs '/data/' participant{p} '/' wave{w} '_' condition{c} session{s} '_' block{b} '_normal/channel_vectorview306_acc1.mat'];
                        copyfile(current_channel_file,[root_dir '/Channel_backup/' participant{p} '/' condition{c} session{s} '_' block{b} '_normal/channel_vectorview306_acc1.mat']);
                    end
                end
            end
        end
    end
end

clearvars('-except', initialVars{:});
disp 'DONE WITH CREATING A SECURITY COPY OF CHANNEL FILES!!!'
disp(datetime)
toc

%% 'Reset' channel files with backup copy if needed

% Option to overwrite current channel files with previous channel info 
% ('reset' to old version saved)

if reset_channel_files == 1
    
tic
disp(' ');      
disp('-------------------------');  
disp('RESSETING CHANNEL FILES');  
disp(datetime)
disp('-------------------------');     
disp(' '); 

for p = 1:length(participant)
      
    % Define sessions within this participant
    pos_par = find(strcmp(session_block_array(:,1), participant{p}(1:4)));
    session = session_block_array{pos_par,2}(1,:);

    for s = 1:length(session)
        
        % Define blocks within this session
        pos_ses = find(strcmp(session_block_array{pos_par,2}(1,:), session{s}));
        block = session_block_array{pos_par,2}{2,pos_ses}(1,:);
        
        for b = 1:length(block)
            for c =1:length(condition)
                for w = 1:length(wave)
                    % If that channel file exist in LLR seg
                    if exist([root_dir_bs '/data/' participant{p} '/' wave{w} '_' condition{c} session{s} '_' block{b} '_normal/channel_vectorview306_acc1.mat'],'file')
                        % If that file exists in the backup folder
                        if exist([root_dir '/Channel_backup/' participant{p} '/' condition{c} session{s} '_' block{b} '_normal/channel_vectorview306_acc1.mat'], 'file')
                            current_channel_file = [root_dir_bs '/data/' participant{p} '/' wave{w} '_' condition{c} session{s} '_' block{b} '_normal/channel_vectorview306_acc1.mat'];
                            backup_channel_file = [root_dir '/Channel_backup/' participant{p} '/' condition{c} session{s} '_' block{b} '_normal/channel_vectorview306_acc1.mat'];
                            % delete current_channel_file before copying backup
                            delete(current_channel_file);
                            % copy backup channel file into equivalent current channel folder
                            copyfile(backup_channel_file,current_channel_file);
                        end
                    end
                end
            end
        end
    end
end

clearvars('-except', initialVars{:});
disp 'DONE WITH RESETING CHANNEL FILES!!!'
disp(datetime)
toc

end

%% Apply tranformation matrices to other blocks and project EEG to surface in ALL blocks

% IMPORTANT: THIS 'NEW' VERSION IS TO BE USED WHEN MANUAL TRANSFORMATION
% WAS USED ONLY IN FIRST BLOCK FIRST CONDITION (THAT IS, WE SAID 'NO' WHEN
% BRAINSTORM ASKS TO COPY THE TRANSFORMATION TO OTHER CONDITIONS OF THE
% SAME BLOCK, BECAUSE IT FAILS SOMETIMES OR COPIES IT TO OTHER BLOCKS
% SOMETIMES AND SOMETIMES DON'T). SO, FOR CONSISTENCY, ONLY FIRST BLOCK
% FIRST CONDITION IS MODIFIED AND COPIED FORM THERE ANYWHERE ELSE

% IMPORTANT: MODIFIED SO THAT IT WILL ONLY COPY TRANSFORMATIONS IF "MANUAL
% TRANSFORMATION" IS FOUND IN LABELS OF MEG. ALSO, IT WILL COPY ALL (OR THE FINAL ONE) 
% MANUAL(NOT OTHERS) TRANSFORMATIONS FOUND

% IMPORTANT: IT WILL ALSO PROJECT TO SURFACE ALL BLOCKS (SO IF FIRST BLOCK
% HAS A TRANSFORMATION, IT WILL COPY THAT TO OTHER BLOCKS AND PROJECT EEG TO
% SURFACE IN ALL BLOCKS; IF FIRST BLOCK HAS NO TRANSFORMATION, THEN IT
% WILL ONLY PROJECT EEG SURFACE TO ALL BLOCKS AND FINISH)

% IMPORTANT! IF FILES ALREADY HAVE EEG moved to the surface it's not a
% problem to move them to surface again if a modification was made in its
% positions (and if no extra modification was made, they will be already in
% surface so when moving again positions should not change). 

tic
disp(' ');      
disp('-------------------------');  
disp('APLYING TRANSFORMATION MATRICES TO OTHER BLOCKS WITHIN SESSION (LLR)');  
disp(datetime)
disp('-------------------------');     
disp(' '); 

no_manual_transf = {}; no_transf_pos = 1;

for p = 1:length(participant)
    
    disp(' ');      
    disp('-------------------------');
    disp(['loading participant ' participant{p}]);
    disp(datetime)
    disp(' '); 
    
    prot_subs = bst_get('ProtocolSubjects');
    current_sub = find(strcmp({prot_subs.Subject.Name}, participant{p}));
    db_reload_conditions(current_sub);
      
    % Define sessions within this participant
    pos_par = find(strcmp(session_block_array(:,1), participant{p}(1:4)));
    session = session_block_array{pos_par,2}(1,:);

    for s = 1:length(session)
        
        % Define blocks within this session
        pos_ses = find(strcmp(session_block_array{pos_par,2}(1,:), session{s}));
        block = session_block_array{pos_par,2}{2,pos_ses}(1,:);
        
        % Find first block of the session and check if it exist
        if ~exist([root_dir_bs '/data/' participant{p} '/LLR_11' session{s} '_' block{1} '_normal'], 'dir') % If folder for condition 11 does not exists
            % Crucially, this is the ONLY folder where you are gonna apply
            % a transformation (exceptions would have been noted)
            error([participant{p} '/LLR_11' session{s} '_' block{1} '_normal does not exist'])
        end

        % Open file from first block to obtain transformation matrix
        % If the folder exist this would exist, if not, an error will pop
        % out so no worries
        load([root_dir_bs '/data/' participant{p} '/LLR_11' session{s} '_' block{1} '_normal/channel_vectorview306_acc1.mat'])

        % Find which columns have manual correction (if any)
        [x,y]= find(contains(TransfMegLabels,'manual correction'));
    
        if isempty(x)
            % If no manual transformation is found, do not copy anything to other folders
            % But keep record of it
            no_manual_transf{no_transf_pos,1} = [participant{p} '/LLR_11' session{s} '_' block{1}]; %#ok<*SAGROW>
            no_transf_pos = no_transf_pos +1;

            % And anyway I need to project EEG electrodes to surface!!

            % Search for all other blocks within the session that are not the first
            % (12, 13 and so on from first block are already changed with the gui)
            % Only needed in 'normal' folders: others don't need to be adjusted
            folders = dir([root_dir_bs '/data/' participant{p} '/']);
            for w = 1:length(wave) % always gonna be LLR now
                for c = 1:length(condition) % it's going to have to be applied to all 15 condition folders
                    for b = 1:length(block) % ALSO FIRST BLOCK NEEDS EEG PROJECTION TO SURFACE
                        results = contains({folders.name},[wave{w} '_' condition{c} session{s} '_' block{b} '_normal']);
                        infolder = find(results);
                        if isempty(infolder) % if no coincidence exists
                            continue 
                        end
                        if size(infolder,2)> 1 % in case of more than one coincidence (should not happen)
                           error(['more than one folder for ' wave{w} '_' condition{c} session{s} '_' block{b} '_normal']);
                        end

                        % Get name of any brainstorm file in the folder (from which
                        % channel file will be retrieved by function)
                        dir_sweeps = dir([root_dir_bs '/data/' participant{p} '/' folders(infolder).name]);
                        sweeps_norm = contains({dir_sweeps.name},'_trial');
                        sweeps_list = find(sweeps_norm);
                        position = sweeps_list(1);  % first trial found
                        file_name = [participant{p} '/' folders(infolder).name '/' dir_sweeps(position).name];
                        sFiles = file_name;

                        % Project electrodes to surface ONLY
                        sFiles = bst_process('CallProcess', 'process_channel_project', sFiles, []);
                    end
                end
            end

        else % there is at least one manual transformation

            % First of all, check if the manual corregistration meets the
            % requirements stablished (e.g. was made in November 2020).
            % Only if this requirement is made will the corregistration be
            % copied from that first block to the others
            % Find position of rows where Align channels manually appears
            pos_manual = find(contains(History(:,3),'Align channels manually:'));
            % Find, whithin these rows, if date coincides with our requirements
            data_manual = find(contains(History(pos_manual,1),string_to_look));
            if ~isempty(data_manual) 
            % if not empty, manual corregistration was applied in selected
            % date, therefore, apply transformation from block 1 to others
           
            
            % First, since original block has a transformation, copy to all
            % other blocks within the session
            
            % Search for all other blocks within the session that are not the first
            % (12, 13 and so on from first block are already changed with the gui)
            % Only needed in 'normal' folders: others don't need to be adjusted
            folders = dir([root_dir_bs '/data/' participant{p} '/']);
            for w = 1:length(wave) % always gonna be LLR now
                for c = 1:length(condition) % it's going to have to be applied to all 15 condition folders
                    for b = 1:length(block) 
                        % All blocks now (because when using the gui we
                        % won't apply to other files, given that freaking
                        % brainstorm sometimes copies them to other blocks
                        % aside the first and sometimes don't!!
                        if c == 1 && b == 1
                            % So first block of each session, condition 11,
                            % are the ones that were modified manually
                            % (none others with new configuration). So if
                            % it's an 11_b1 (or whichever block is first in
                            % session_block_array), then continue to next
                            % block
                            continue
                        end
                        results = contains({folders.name},[wave{w} '_' condition{c} session{s} '_' block{b} '_normal']);
                        infolder = find(results);
                        if isempty(infolder) % if no coincidence exists
                            continue 
                        end
                        if size(infolder,2)> 1 % in case of more than one coincidence (should not happen)
                           error(['more than one folder for ' wave{w} '_' condition{c} session{s} '_' block{b} '_normal']);
                        end

                        % Get target channel file name
                        target_channel_name = [root_dir_bs '/data/' participant{p} '/' folders(infolder).name '/channel_vectorview306_acc1.mat'];

                        if copy_all_transf == 1 % if choosen to copy all
                            % Apply each of the transformations found in the original file
                            for lt = 1:length(x) % for each of the transformations found
                                % Select column of transformation
                                column_transf = y(lt); % Meg always
                                % Meg always, because EEG can have the 0 and 1 matrix when projecting electrodes to surface
                                Transf = TransfMeg{1,column_transf};
                                % Apply transformation (channel_apply_transf is a brainstorm function)
                                ChannelMats = channel_apply_transf(target_channel_name,  Transf);
                                % It's automatically saved
                            end
                        else % if choosen to copy last one
                            % Select last column of transformation
                            column_transf = y(end); % Meg always
                            % Meg always, because EEG can have the 0 and 1 matrix when projecting electrodes to surface
                            Transf = TransfMeg{1,column_transf};
                            % Apply transformation (channel_apply_transf is a brainstorm function)
                            ChannelMats = channel_apply_transf(target_channel_name,  Transf);
                            % It's automatically saved
                        end   
                    end
                end
            end      
            
            % Second, IN ALL blocks (including first), project EEG electrodes
            % to surface
            
            % Search for all other blocks within the session that are not the first
            % (12, 13 and so on from first block are already changed with the gui)
            % Only needed in 'normal' folders: others don't need to be adjusted
            folders = dir([root_dir_bs '/data/' participant{p} '/']);
            for w = 1:length(wave) % always gonna be LLR now
                for c = 1:length(condition) % it's going to have to be applied to all 15 condition folders
                    for b = 1:length(block) % ALSO FIRST BLOCK NEEDS EEG PROJECTION TO SURFACE
                        results = contains({folders.name},[wave{w} '_' condition{c} session{s} '_' block{b} '_normal']);
                        infolder = find(results);
                        if isempty(infolder) % if no coincidence exists
                            continue 
                        end
                        if size(infolder,2)> 1 % in case of more than one coincidence (should not happen)
                           error(['more than one folder for ' wave{w} '_' condition{c} session{s} '_' block{b} '_normal']);
                        end

                        % Get name of any brainstorm file in the folder (from which
                        % channel file will be retrieved by function)
                        dir_sweeps = dir([root_dir_bs '/data/' participant{p} '/' folders(infolder).name]);
                        sweeps_norm = contains({dir_sweeps.name},'_trial');
                        sweeps_list = find(sweeps_norm);
                        position = sweeps_list(1);  % first trial found
                        file_name = [participant{p} '/' folders(infolder).name '/' dir_sweeps(position).name];
                        sFiles = file_name;

                        % Project electrodes to surface
                        sFiles = bst_process('CallProcess', 'process_channel_project', sFiles, []);
                    end
                end
            end
            end
        end
    end
end

save([root_dir '/no_manual_transf.mat'],'no_manual_transf');

clearvars('-except', initialVars{:});
disp 'DONE WITH APLYING TRANSFORMATION MATRICES TO OTHER BLOCKS WITHIN SESSION (LLR)!!!'
disp(datetime)
toc

%% Copy channel files with corregistration from Backup folder

tic
disp(' ');      
disp('-------------------------');  
disp('COPYING TRANSFORMATION MATRICES FROM BACKUP FOLDER');  
disp(datetime)
disp('-------------------------');     
disp(' '); 

root_dir_destiny = '~/brainstorm_db/Project_LLRseg';
root_dir_origin = backup_corregistration_folder;
% backup_corregistration_folder = '/private/path/project/Channel_backup/EEG_equal_projected_February_2021';
channel_file_copy_log = {};
log_channel = 1;

for p = 1:length(participant)
      
    % Define sessions within this participant
    pos_par = find(strcmp(session_block_array(:,1), participant{p}(1:4)));
    session = session_block_array{pos_par,2}(1,:);

    for s = 1:length(session)
        
        % Define blocks within this session
        pos_ses = find(strcmp(session_block_array{pos_par,2}(1,:), session{s}));
        block = session_block_array{pos_par,2}{2,pos_ses}(1,:);
        
        for b = 1:length(block)
            for c =1:length(condition)
                % If that channel file exist in LLR seg
                if exist([root_dir_destiny '/data/' participant{p} '/LLR_' condition{c} session{s} '_' block{b} '_normal/channel_vectorview306_acc1.mat'],'file')
                    % AND that very same channel file exist in backup
                    if exist([root_dir_origin '/' participant{p} '/' condition{c} session{s} '_' block{b} '_normal/channel_vectorview306_acc1.mat'],'file')
                        channel_file_destiny = [root_dir_destiny '/data/' participant{p} '/LLR_' condition{c} session{s} '_' block{b} '_normal/channel_vectorview306_acc1.mat'];
                        channel_file_origin = [root_dir_origin '/' participant{p} '/' condition{c} session{s} '_' block{b} '_normal/channel_vectorview306_acc1.mat'];
                        % delete original LLR channel file
                        delete(channel_file_destiny);
                        % copy LLR backup channel file into equivalent LLR folder
                        copyfile(channel_file_origin,channel_file_destiny);
                        % it will contain all transformations
                    else
                        % Destiny but not origin exists
                        channel_file_copy_log{log_channel,1} = ['Destiny but not origin in ' condition{c} session{s} '_' block{b} '_normal'];
                        log_channel = log_channel + 1; 
                    end
                else
                    if exist([root_dir_origin '/' participant{p} '/' condition{c} session{s} '_' block{b} '_normal/channel_vectorview306_acc1.mat'],'file')
                        channel_file_copy_log{log_channel,1} = ['Origin but not destiny in ' condition{c} session{s} '_' block{b} '_normal'];
                        log_channel = log_channel + 1; 
                    else
                        channel_file_copy_log{log_channel,1} = ['not origin or destiny in ' condition{c} session{s} '_' block{b} '_normal'];
                        log_channel = log_channel + 1; 
                    end
                end
            end
        end
    end
end
save([root_dir '/Events/channel_file_copy_log.mat'],'channel_file_copy_log');

clearvars('-except', initialVars{:});
disp 'DONE WITH COPYING TRANSFORMATION MATRICES FROM LLR TO MLR!!!'
disp(datetime)
toc

%% Copy the bad channels selected from LLR new protocol to LLR old protocol before EEG sources

tic
disp(' ');      
disp('-------------------------');  
disp('COPYING BAD CHANNELS FROM HUMON_LLR TO LLRseg before EEG sources');  
disp(datetime)
disp('-------------------------');     
disp(' '); 


root_dir_Origin = [root_dir '/brainstorm_db/Project_LLR']; 
root_dir_Destiny = [root_dir '/brainstorm_db/LLRseg']; 
Copy_failures = {}; copy_failure_count = 0;

% For all LLR files
for p = 1:length(participant)
      
    % Define sessions within this participant
    pos_par = find(strcmp(session_block_array(:,1), participant{p}(1:4)));
    session = session_block_array{pos_par,2}(1,:);

    for s = 1:length(session)
        
        % Define blocks within this session
        pos_ses = find(strcmp(session_block_array{pos_par,2}(1,:), session{s}));
        block = session_block_array{pos_par,2}{2,pos_ses}(1,:);
        
        for c = 1:length(condition)
            for b = 1:length(block) 
            
            % If that file does not exist in HUMON_LLR
            Origin_file_base = dir([root_dir_Origin '/data/' participant{p} '/@raw' session{s} '_' block{b} '/data_0raw*']);
            if length(Origin_file_base)>1;continue;end
            Origin_file = [Origin_file_base.folder '/' Origin_file_base.name];
            
            %%%%%%%%%%%%%% Define destiny files
            folders = dir([root_dir_Destiny '/data/' participant{p} '/']);
            results = contains({folders.name},['LLR_' condition{c} session{s} '_' block{b} '_normal']) | ...
                contains({folders.name},['LLR_' condition{c} session{s} '_' block{b} '_covariance']);
            infolder = find(results);
                    
            if isempty(infolder)
                continue
            end
            
            for l = 1:length(infolder) % for each folder found
                Destiny_files = {};
                line = infolder(l); % first position index and so on
                dir_sweeps = dir([root_dir_Destiny '/data/' participant{p} '/' folders(line).name]);
                sweeps_norm = startsWith({dir_sweeps.name},'data_');
                sweeps_list = find(sweeps_norm);
                for j= 1:length(sweeps_list) % anything that has data on it
                    position = sweeps_list(j);
                    Destiny_files{j} = [root_dir_Destiny '/data/' participant{p} '/' folders(line).name '/' dir_sweeps(position).name];
                end
            if ~exist(Origin_file,'file')
                Copy_failures{copy_failure_count,1} = [participant{p} '_' session{s} '_' block{b}]; %#ok<*SAGROW>
                Copy_failures{copy_failure_count,2} = 'No Origin file';
                copy_failure_count = copy_failure_count +1;
                continue
            elseif isempty(Destiny_files)
                Copy_failures{copy_failure_count,1} = [participant{p} '_' session{s} '_' block{b}];
                Copy_failures{copy_failure_count,2} = 'No Destiny file';
                copy_failure_count = copy_failure_count +1;
                continue
            else
            disp(' ');      
            disp('-------------------------');  
            disp(['Copying bad chans from New to old protocol for ' participant{p} '_' session{s} '_' block{b} '_' condition{c}]);
            disp(datetime)
            disp(' '); 
            
            for df = 1:length(Destiny_files)
                current_destiny_file = Destiny_files{df};
                eval(['load ' Origin_file]) % Load origin to retrieve ChannelFlag to copy
                ChannelFlag_copy = ChannelFlag; 
                eval(['load ' current_destiny_file]) % Load destiny file
                ChannelFlag = ChannelFlag_copy; % Copy channel flag from origin file
                delete(current_destiny_file) % Delete destiny file before copy
                save(current_destiny_file, 'ChannelFlag', 'ColormapType', 'Comment', 'DataType', 'Device', 'DisplayUnits', 'Events', 'F', 'History', 'Leff', 'nAvg', 'Std', 'Time')
            end
            end            
            end     
            end
        end
    end
end

% % For Intra files
% for p = 1:length(participant)
%       
%     % Define sessions within this participant
%     pos_par = find(strcmp(session_block_array(:,1), participant{p}(1:4)));
%     session = session_block_array{pos_par,2}(1,:);
% 
%     for s = 1:length(session)
%         
%         % Define blocks within this session
%         pos_ses = find(strcmp(session_block_array{pos_par,2}(1,:), session{s}));
%         block = session_block_array{pos_par,2}{2,pos_ses}(1,:);
%         
%         for c = 1:length(condition)
%             for b = 1:length(block) 
%             
%             % If that file does not exist in HUMON_LLR
%             Origin_file_base = dir([root_dir_Origin '/data/' participant{p} '/@raw' session{s} '_' block{b} '/data_0raw*']);
%             if length(Origin_file_base)>1;continue;end
%             Origin_file = [Origin_file_base.folder '\' Origin_file_base.name];
%             
%             %%%%%%%%%%%%%% Define destiny files
%             folders = dir([root_dir_Destiny '/data/' participant{p} '/']);
%             results = strcmp({folders.name},'@intra');
%             infolder = find(results);
%                     
%             if isempty(infolder)
%                 continue
%             end
%             
%             for l = 1:length(infolder) % for each folder found
%                 Destiny_files = {};
%                 line = infolder(l); % first position index and so on
%                 dir_sweeps = dir([root_dir_Destiny '/data/' participant{p} '/' folders(line).name]);
%                 sweeps_norm = startsWith({dir_sweeps.name},'data_');
%                 sweeps_list = find(sweeps_norm);
%                 for j= 1:length(sweeps_list) % anything that has data on it
%                     position = sweeps_list(j);
%                     Destiny_files{j} = [root_dir_Destiny '/data/' participant{p} '/' folders(line).name '/' dir_sweeps(position).name];
%                 end
%             if ~exist(Origin_file,'file')
%                 Copy_failures{copy_failure_count,1} = [participant{p} '_' session{s} '_' block{b}]; %#ok<*SAGROW>
%                 Copy_failures{copy_failure_count,2} = 'No Origin file';
%                 copy_failure_count = copy_failure_count +1;
%                 continue
%             elseif isempty(Destiny_files)
%                 Copy_failures{copy_failure_count,1} = [participant{p} '_' session{s} '_' block{b}];
%                 Copy_failures{copy_failure_count,2} = 'No Destiny file';
%                 copy_failure_count = copy_failure_count +1;
%                 continue
%             else
%             disp(' ');      
%             disp('-------------------------');  
%             disp(['Copying bad chans from New to old protocol for ' participant{p} '_' session{s} '_' block{b} '_' condition{c}]);
%             disp(datetime)
%             disp(' '); 
%             
%             for df = 1:length(Destiny_files)
%                 current_destiny_file = Destiny_files{df};
%                 eval(['load ' Origin_file]) % Load origin to retrieve ChannelFlag to copy
%                 ChannelFlag_copy = ChannelFlag; 
%                 eval(['load ' current_destiny_file]) % Load destiny file
%                 ChannelFlag = ChannelFlag_copy; % Copy channel flag from origin file
%                 delete(current_destiny_file) % Delete destiny file before copy
%                 save(current_destiny_file, 'ChannelFlag', 'ColormapType', 'Comment', 'DataType', 'Device', 'DisplayUnits', 'Events', 'F', 'History', 'Leff', 'nAvg', 'Std', 'Time')
%             end
%             end            
%             end     
%             end
%         end
%     end
% end


Copy_failures_badchans = Copy_failures;
save([root_dir '/Channel_backup/Copy_failures_badchans.mat'],'Copy_failures_badchans');
disp 'DONE WITH COPYING BAD CHANNELS FROM HUMON_LLR TO LLRseg!!!'
disp(datetime)
toc

%% Calculate noise covariance

% If wanting to do sources with EEG, MEG or both we should repeat this step
% which whatever sensor_analysis variable, since there cannot be more than
% one noise covariance in the folder. It will overwrite previous automatically

if covar_epochs == 1

tic
disp(' ');      
disp('-------------------------');  
disp('CALCULATING NOISE COVARIANCE'); 
disp(datetime)
disp('-------------------------');     
disp(' ');

% Not necessary, since it's better to have noise covariance with both
% EEG and MEG, then use whichever needed in the inverse solution
% if source_analysis == 1
%     modality_data = {'EEG'};
%     sensortype_cov = 'EEG';
% elseif source_analysis == 2
%     modality_data = {'MEG'};
%     sensortype_cov = 'MEG';
% else % Because if it's 4 we still want to use MEG and EEG for the matrix
%     modality_data = {'both_mod'};
%     sensortype_cov = 'MEG, EEG';
% end
    
% Compute matrices
for p = 1:length(participant)
    eval(['Insuff_cov_trials_LLR_' participant{p} '= {};'])
    count_insuff_cov = 1;
    
    disp(' ');      
    disp('-------------------------');
    disp(['loading participant ' participant{p}]);
    disp(datetime)
    disp(' '); 
    
    prot_subs = bst_get('ProtocolSubjects');
    current_sub = find(strcmp({prot_subs.Subject.Name}, participant{p}));
    db_reload_conditions(current_sub);
    
    % Define sessions within this participant
    pos_par = find(strcmp(session_block_array(:,1), participant{p}(1:4)));
    session = session_block_array{pos_par,2}(1,:);
    
    for s = 1:length(session)
        
        % Define blocks within this session
        pos_ses = find(strcmp(session_block_array{pos_par,2}(1,:), session{s}));
        block = session_block_array{pos_par,2}{2,pos_ses}(1,:);
        
        for b = 1:length(block)
            for w = 1:length(wave)
                for mode = 1:length(modality_data) % EEG, MEG and both mod    
                    folders = dir([root_dir_bs '/data/' participant{p} '/']);
                    results = contains({folders.name},[wave{w} '_91' session{s} '_' block{b} '_covariance'])...
                        | contains({folders.name},[wave{w} '_92' session{s} '_' block{b} '_covariance'])...
                        | contains({folders.name},[wave{w} '_93' session{s} '_' block{b} '_covariance']);
                    infolder = find(results);
                        if isempty(infolder)
                            continue
                        end
                    file_names = {}; rej_names = {};
                    count = 0; count_rej = 0; % dirty way to put files from different conditions in the same file_names structure
                    for l= 1:length(infolder) % For each of the folders (91, 92 and 93 sweeps togheter)
                        line = infolder(l);
                        dir_sweeps = dir([root_dir_bs '/data/' participant{p} '/' folders(line).name]);
                        sweeps_cov = contains({dir_sweeps.name},'_trial');
                        sweeps_count = find(sweeps_cov);
                        if isempty(sweeps_count) % if the folder exists but there are no trials
                            continue
                        end

                        % FIND REJECTED TRIALS
                        load([root_dir '/Events/BadTrials_Cov/' ProtocolName(7:9) '/' modality_data{w} '/' participant{p} '/' folders(line).name '/BadTrials.mat']); 
                        if ~isempty(BadTrials)
                            for rej=1:length(BadTrials)
                                rej_names{count_rej + rej} = BadTrials{rej};
                            end
                            count_rej = count_rej + rej;
                        end
                        % FIND ALL TRIALS
                        for j= 1:length(sweeps_count)
                            position = sweeps_count(j);
                            file_names{count + j} = [participant{p} '/' folders(line).name '/' dir_sweeps(position).name];
                        end
                        count = count + j;
                    end
                    % REMOVE REJECTED TRIALS FROM ALL TRIALS LIST (IF THEY EXIST)
                    if ~isempty(rej_names)
                        rejected_idx = [];
                        for rt=1:length(rej_names)
                        rejected = contains(file_names, rej_names{rt});
                        rejected_idx(rt) = find(rejected);
                        end
                        file_names(rejected_idx) = [];
                    end

                    sFiles = file_names; % ALL 91, 92 AND 93 from a particular session, block mod and wave
                    
                    % Minimum to calculate noise cov, or else we don´t
                    % get enough samples according to (N*(N+1)/2 samples)
                    % considering sampling rate and modality
                    if source_analysis == 1 % EEG
                        if size(sFiles,2) < 1
                            eval(['Insuff_cov_trials_LLR_' participant{p} '{count_insuff_cov,1} = [wave{w} ''_'' session{s} ''_'' modality_data{mode} ''_'' block{b}];'])
                            count_insuff_cov = count_insuff_cov + 1;
                            continue
                        end
                    elseif source_analysis == 2 % MEG
                        if size(sFiles,2) < 10
                            eval(['Insuff_cov_trials_LLR_' participant{p} '{count_insuff_cov,1} = [wave{w} ''_'' session{s} ''_'' modality_data{mode} ''_'' block{b}];'])
                            count_insuff_cov = count_insuff_cov + 1;
                            continue
                        end
                    elseif source_analysis == 3 % Both mod
                        if size(sFiles,2) < 14
                            eval(['Insuff_cov_trials_LLR_' participant{p} '{count_insuff_cov,1} = [wave{w} ''_'' session{s} ''_'' modality_data{mode} ''_'' block{b}];'])
                            count_insuff_cov = count_insuff_cov + 1;
                            continue
                        end
                    end

                    disp(' ');      
                    disp('-------------------------');  
                    disp(['Computing noise covariance for ' participant{p} '_' session{s} '_' modality_data{mode} '_' block{b} '_' wave{w}]);
                    disp(datetime)
                    disp(' ');

                    % Process: Compute noise covariance
                    sFiles = bst_process('CallProcess', 'process_noisecov', sFiles, [], ...
                        'baseline',       time_noise_cov, ...
                        'datatimewindow', [0, 0], ...
                        'sensortypes',    'MEG, EEG', ... % always both
                        'target',         1, ...  % Noise covariance (covariance over baseline time window)
                        'dcoffset',       1, ...  % Block by block, to avoid effects of slow shifts in data
                        'identity',       0, ...
                        'copycond',       0, ...
                        'copysubj',       0, ...
                        'copymatch',      0, ...
                        'replacefile',    1);  % Replace

                    % Process: Add tag
                    % NOT NECESSARY, BECAUSE THERE CAN BE ONLY NOISE
                    % COVARIANCE AT THE SAME TIME AND IT HAS TO BE EITHER
                    % EEG OR MEG OR BOTH, its always gonna be called noisecov_full
                end
            end
        end
    end
    save([root_dir '/Sources/Insuff_cov_trials_LLR_' participant{p} '.mat'],['Insuff_cov_trials_LLR_' participant{p}]);
end

% Copy matrices among conditions
for p = 1:length(participant)

    % Define sessions within this participant
    pos_par = find(strcmp(session_block_array(:,1), participant{p}(1:4)));
    session = session_block_array{pos_par,2}(1,:);
    
    for s = 1:length(session)
        
        % Define blocks within this session
        pos_ses = find(strcmp(session_block_array{pos_par,2}(1,:), session{s}));
        block = session_block_array{pos_par,2}{2,pos_ses}(1,:);
        for b = 1:length(block) 
            for w = 1:length(wave)           
            % Copy noise cov from 91 (any from 91, 92 and 93 works) to the rest of folders
            if ~exist([root_dir_bs '/data/' participant{p} '/' wave{w} '_91' session{s} '_' block{b} '_covariance/noisecov_full.mat'], 'file')
               continue
            end

            disp(' ');      
            disp('-------------------------');  
            disp(['Copying noise cov among conditions for ' participant{p} '_' session{s} '_' block{b} '_' wave{w}]);
            disp(datetime)
            disp(' ');

            % Theoretically, if it's not in 91 shouldn´t be in 92 or 93,
            % but just in case there were enough 92s trials, but not 91s,
            % for instance. Also, if no file is there, don´t try to copy
            if exist([root_dir_bs '/data/' participant{p} '/' wave{w} '_91' session{s} '_' block{b} '_covariance/noisecov_full.mat'], 'file')
                file_name_noise = [root_dir_bs '/data/' participant{p} '/' wave{w} '_91' session{s} '_' block{b} '_covariance/noisecov_full.mat'];
            elseif exist([root_dir_bs '/data/' participant{p} '/' wave{w} '_92' session{s} '_' block{b} '_covariance/noisecov_full.mat'], 'file')
                file_name_noise = [root_dir_bs '/data/' participant{p} '/' wave{w} '_92' session{s} '_' block{b} '_covariance/noisecov_full.mat'];
            elseif exist([root_dir_bs '/data/' participant{p} '/' wave{w} '_93' session{s} '_' block{b} '_covariance/noisecov_full.mat'], 'file')
                file_name_noise = [root_dir_bs '/data/' participant{p} '/' wave{w} '_93' session{s} '_' block{b} '_covariance/noisecov_full.mat'];
            else
                continue
            end
            for c = 1:length(condition)
                % Check if destiny folder exist
                if exist([root_dir_bs '/data/' participant{p} '/' wave{w} '_' condition{c} session{s} '_' block{b} '_normal/'], 'dir')
                    copyfile(file_name_noise,[root_dir_bs '/data/' participant{p} '/' wave{w} '_' condition{c} session{s} '_' block{b} '_normal/noisecov_full.mat']);
                end
            end
            end
        end 
    end
end

% reset modality_data variable to its original
modality_data = {'EEG','MEG','both_mod'};

clearvars('-except', initialVars{:});
disp 'DONE WITH CALCULATING NOISE COVARIANCE!!!'
disp(datetime)
toc

end

%% If needed before next: copy combined EEG/MEG forward models from new protocol

% Since forward models for some subjects are already ran there, and the
% positions of the sensors are always the same, you can copy leadfields

% WARNING: if corregistration has been done for channel files in one but
% not the other protocol, you should make them equal first


tic
disp(' ');      
disp('-------------------------');  
disp('COPYING COMBINED EEG-MEG LEADFIELD FROM HUMON_LLR TO LLRseg');  
disp(datetime)
disp('-------------------------');     
disp(' '); 


root_dir_Origin = [root_dir '/brainstorm_db/Project_LLR']; 
root_dir_Destiny = [root_dir '/brainstorm_db/LLRseg']; 

% For all LLR files
for p = 1:length(participant)
      
    % Define sessions within this participant
    pos_par = find(strcmp(session_block_array(:,1), participant{p}(1:4)));
    session = session_block_array{pos_par,2}(1,:);

    for s = 1:length(session)
        
        % Define blocks within this session
        pos_ses = find(strcmp(session_block_array{pos_par,2}(1,:), session{s}));
        block = session_block_array{pos_par,2}{2,pos_ses}(1,:);
        
        for b = 1:length(block) 
            
        % If destiny directory or origin file do not exist, abort    
        if ~exist([root_dir_Destiny '/data/' participant{p} '/LLR_11' session{s} '_' block{b} '_normal'], 'dir'); continue;end    
        if ~exist([root_dir_Origin '/data/' participant{p} '/' session{s} '_' block{b} '/headmodel_surf_duneuro_os_meg.mat'],'file'); continue; end
        
        disp(' ');      
        disp('-------------------------');  
        disp(['Copying leadfields from New to old protocol for ' participant{p} '_' session{s} '_' block{b}]);
        disp(datetime)
        disp(' '); 
        
        % Get the combined EEG/MEG forward model from HUMON_LLR
        Origin_file = [root_dir_Origin '/data/' participant{p} '/' session{s} '_' block{b} '/headmodel_surf_duneuro_os_meg.mat'];
        
        % Before copy, delete any previous leadfields there
        if exist([root_dir_Destiny '/data/' participant{p} '/LLR_11' session{s} '_' block{b} '_normal/headmodel_surf_duneuro_os_meg.mat'],'file')
            delete([root_dir_Destiny '/data/' participant{p} '/LLR_11' session{s} '_' block{b} '_normal/headmodel_surf_duneuro_os_meg.mat'])
        end    
        if exist([root_dir_Destiny '/data/' participant{p} '/LLR_11' session{s} '_' block{b} '_normal/headmodel_surf_duneuro.mat'],'file')
            delete([root_dir_Destiny '/data/' participant{p} '/LLR_11' session{s} '_' block{b} '_normal/headmodel_surf_duneuro.mat'])
        end  
        if exist([root_dir_Destiny '/data/' participant{p} '/LLR_11' session{s} '_' block{b} '_normal/headmodel_surf_os_meg.mat'],'file')
            delete([root_dir_Destiny '/data/' participant{p} '/LLR_11' session{s} '_' block{b} '_normal/headmodel_surf_os_meg.mat'])
        end  
        
        copyfile(Origin_file,[root_dir_Destiny '/data/' participant{p} '/LLR_11' session{s} '_' block{b} '_normal/headmodel_surf_duneuro_os_meg.mat']);  
        end
    end
end

disp 'DONE WITH COPYING LEADFIELDS FROM HUMON_LLR TO LLRseg!!!'
disp(datetime)
toc

%% Calculate forward model (HAS TO BE IN 98)
% Compute leadfields only once per each participant, session and block
tic
disp(' ');      
disp('-------------------------');
disp('COMPUTING FORWARD MODEL (EEG and MEG) FOR EACH PARTICIPANT, SESSION AND BLOCK (LLR)');
disp(datetime)
disp('-------------------------');
disp(' ');

% This is no longer needed as we are always going to compute a forward
% model that contains both
% if source_analysis == 1
%     modality_data = {'EEG'};
% elseif source_analysis == 2
%     modality_data = {'MEG'};
% elseif source_analysis == 3 % if it's combined EEG/MEG we need the forward model for EEG
%     modality_data = {'EEG'};
% end

No_head_model_created = {};
count_no_cov = 1;

% Compute one EEG leadfield per session and one MEG leadfield per block
for p = 1:length(participant)
    
    disp(' ');      
    disp('-------------------------');
    disp(['loading participant ' participant{p}]);
    disp(datetime)
    disp(' '); 
    
    prot_subs = bst_get('ProtocolSubjects');
    current_sub = find(strcmp({prot_subs.Subject.Name}, participant{p}));
    db_reload_conditions(current_sub);
    % Reload anatomy too
    db_reload_subjects(current_sub);
    
    
    % Define sessions within this participant
    pos_par = find(strcmp(session_block_array(:,1), participant{p}(1:4)));
    session = session_block_array{pos_par,2}(1,:);
    
    for s = 1:length(session)
        
        % Define blocks within this session
        pos_ses = find(strcmp(session_block_array{pos_par,2}(1,:), session{s}));
        block = session_block_array{pos_par,2}{2,pos_ses}(1,:);
        
        for b = 1:length(block)
            
            % We could chose any other, but let's take the LLR_11
            if ~exist([root_dir_bs '/data/' participant{p} '/LLR_11' session{s} '_' block{b} '_normal'], 'dir') 
                continue
            end

            folders = dir([root_dir_bs '/data/' participant{p} '/LLR_11' session{s} '_' block{b} '_normal']);            
            % Before it had EEG or MEG in the second contains, but that's
            % not needed, any file will have the same EEG and MEG sensor
            % positions
            results =  contains({folders.name},['LLR_11' session{s} '_' block{b} '_normal']) & contains({folders.name}, '_average_source');
            infolder = find(results);
            if isempty(infolder)
                % In case folder exists but does not contain anything (channel vector or average)
                No_head_model_created{count_no_cov,1} = [participant{p} '_LLR_11' session{s} '_' block{b} '_normal_source'];
                count_no_cov = count_no_cov +1;
                continue
            end
            if length(infolder) > 1 % Most likely, since it will find EEG and MEG
                % This is just in case it happens that EEG file is present
                % but not the MEG or viceversa, for the purpose of
                % computing forward model it does not matter which, but it
                % has to be only one
                infolder = infolder(1);
            end
            
            file_name = [participant{p} '/LLR_11' session{s} '_' block{b} '_normal/' folders(infolder).name];
            sFiles_orig = file_name;
            
            if isempty(sFiles_orig)
                error('for whatever reason, sFiles is empty, my friend...')
            end
            
            % HERE WE MAY COMPUTE A VOLUME (VS CORTEX) SOURCE SPACE for MLR
                        
            if b == 1 % If it's first block of session, compute both EEG and MEG forward models
            
            % OVERLAPPING SPHERES   
            disp(' ');
            disp('-------------------------');  
            disp(['Computing forward model (Overlapping spheres & Cortical surface) for ' participant{p} '_' session{s} '_' block{b}]);
            disp(datetime)
            disp(' ');

            % If stated, find and delete any previous head model (overwrite it)  
            if delete_previous_file == 1
                % check if there is already a headmodel OF THIS TYPE in the folder
                if isfile([root_dir_bs '/data/' participant{p} '/LLR_11' session{s} '_' block{b} '_normal/headmodel_surf_os_meg.mat'])
                    delete([root_dir_bs '/data/' participant{p} '/LLR_11' session{s} '_' block{b} '_normal/headmodel_surf_os_meg.mat'])
                end
                % head model = headmodel_surf_os_meg.mat
                % source = contains KERNEL endswith _Source_Regul_MEG.mat                 
            end

            % Process: Compute head model
            sFiles_MEG = bst_process('CallProcess', 'process_headmodel', sFiles_orig, [], ...
                'Comment',     '', ...
                'sourcespace', 1, ...  % Cortex surface
                'volumegrid',  struct(...
                     'Method',        'isotropic', ...
                     'nLayers',       17, ...
                     'Reduction',     3, ...
                     'nVerticesInit', 4000, ...
                     'Resolution',    0.005, ...
                     'FileName',      []), ...
                'meg',         3, ...  % Overlapping spheres
                'eeg',         1, ...  % 
                'ecog',        1, ...  % 
                'seeg',        1, ...  % 
                'openmeeg',    struct(...
                     'BemFiles',     {{}}, ...
                     'BemNames',     {{'Scalp', 'Skull', 'Brain'}}, ...
                     'BemCond',      [1, 0.0125, 1], ...
                     'BemSelect',    [1, 1, 1], ...
                     'isAdjoint',    0, ...
                     'isAdaptative', 1, ...
                     'isSplit',      0, ...
                     'SplitLength',  4000));
                 
            % FEM DUNNEURO SIMNIBS     
            disp(' ');
            disp('-------------------------');  
            disp(['Computing forward model (DUNNEURO SIMNIBS) for ' participant{p} '_' session{s} '_' block{b}]);
            disp(datetime)
            disp(' ');   
            
            % If stated, find and delete any previous head model (overwrite it)  
            if delete_previous_file == 1
                % check if there is already a headmodel OF THIS TYPE in the folder
                if isfile([root_dir_bs '/data/' participant{p} '/LLR_11' session{s} '_' block{b} '_normal/headmodel_surf_duneuro.mat'])
                    delete([root_dir_bs '/data/' participant{p} '/LLR_11' session{s} '_' block{b} '_normal/headmodel_surf_duneuro.mat'])
                end
                % head model = headmodel_surf_os_meg.mat
                % source = contains KERNEL endswith _Source_Regul_MEG.mat                 
            end
                 
            % FEM FROM DUNNEURO SIMNIBS
            sFiles_EEG = bst_process('CallProcess', 'process_headmodel', sFiles_orig, [], ...
                'Comment',     '', ...
                'sourcespace', 1, ...  % Cortex surface
                'meg',         1, ...  % 
                'eeg',         4, ...  % DUNEuro FEM
                'ecog',        1, ...  % 
                'seeg',        1, ...  % 
                'duneuro',     struct(...
                     'FemCond',             [0.14, 0.33, 1.79, 0.008, 0.43], ...
                     'FemSelect',           [1, 1, 1, 1, 1], ...
                     'UseTensor',           0, ...
                     'Isotropic',           1, ...
                     'SrcShrink',           0, ...
                     'SrcForceInGM',        1, ...
                     'FemType',             'fitted', ...
                     'SolverType',          'cg', ...
                     'GeometryAdapted',     0, ...
                     'Tolerance',           1e-08, ...
                     'ElecType',            'normal', ...
                     'MegIntorderadd',      0, ...
                     'MegType',             'physical', ...
                     'SolvSolverType',      'cg', ...
                     'SolvPrecond',         'amg', ...
                     'SolvSmootherType',    'ssor', ...
                     'SolvIntorderadd',     0, ...
                     'DgSmootherType',      'ssor', ...
                     'DgScheme',            'sipg', ...
                     'DgPenalty',           20, ...
                     'DgEdgeNormType',      'houston', ...
                     'DgWeights',           1, ...
                     'DgReduction',         1, ...
                     'SolPostProcess',      1, ...
                     'SolSubstractMean',    0, ...
                     'SolSolverReduction',  1e-10, ...
                     'SrcModel',            'venant', ...
                     'SrcIntorderadd',      0, ...
                     'SrcIntorderadd_lb',   2, ...
                     'SrcNbMoments',        3, ...
                     'SrcRefLen',           20, ...
                     'SrcWeightExp',        1, ...
                     'SrcRelaxFactor',      6, ...
                     'SrcMixedMoments',     1, ...
                     'SrcRestrict',         1, ...
                     'SrcInit',             'closest_vertex', ...
                     'BstSaveTransfer',     0, ...
                     'BstEegTransferFile',  'eeg_transfer.dat', ...
                     'BstMegTransferFile',  'meg_transfer.dat', ...
                     'BstEegLfFile',        'eeg_lf.dat', ...
                     'BstMegLfFile',        'meg_lf.dat', ...
                     'UseIntegrationPoint', 1, ...
                     'EnableCacheMemory',   0, ...
                     'MegPerBlockOfSensor', 0), ...
                'channelfile', '');

            else % It's not the first block, so compute only the overlapping spheres for MEG
            
            % OVERLAPPING SPHERES   
            disp(' ');
            disp('-------------------------');  
            disp(['Computing forward model (Overlapping spheres & Cortical surface) for ' participant{p} '_' session{s} '_' block{b}]);
            disp(datetime)
            disp(' ');

            % If stated, find and delete any previous head model (overwrite it)  
            if delete_previous_file == 1
                % check if there is already a headmodel OF THIS TYPE in the folder
                if isfile([root_dir_bs '/data/' participant{p} '/LLR_11' session{s} '_' block{b} '_normal/headmodel_surf_os_meg.mat'])
                    delete([root_dir_bs '/data/' participant{p} '/LLR_11' session{s} '_' block{b} '_normal/headmodel_surf_os_meg.mat'])
                end
                % head model = headmodel_surf_os_meg.mat
                % source = contains KERNEL endswith _Source_Regul_MEG.mat                 
            end

            % Process: Compute head model
            sFiles_MEG = bst_process('CallProcess', 'process_headmodel', sFiles_orig, [], ...
                'Comment',     '', ...
                'sourcespace', 1, ...  % Cortex surface
                'volumegrid',  struct(...
                     'Method',        'isotropic', ...
                     'nLayers',       17, ...
                     'Reduction',     3, ...
                     'nVerticesInit', 4000, ...
                     'Resolution',    0.005, ...
                     'FileName',      []), ...
                'meg',         3, ...  % Overlapping spheres
                'eeg',         1, ...  % 
                'ecog',        1, ...  % 
                'seeg',        1, ...  % 
                'openmeeg',    struct(...
                     'BemFiles',     {{}}, ...
                     'BemNames',     {{'Scalp', 'Skull', 'Brain'}}, ...
                     'BemCond',      [1, 0.0125, 1], ...
                     'BemSelect',    [1, 1, 1], ...
                     'isAdjoint',    0, ...
                     'isAdaptative', 1, ...
                     'isSplit',      0, ...
                     'SplitLength',  4000));
            end
        end
    end
end

% Copy EEG forward model from block 1 to next
for p = 1:length(participant)
        
    % Define sessions within this participant
    pos_par = find(strcmp(session_block_array(:,1), participant{p}(1:4)));
    session = session_block_array{pos_par,2}(1,:);
    
    for s = 1:length(session)
        
        % Define blocks within this session
        pos_ses = find(strcmp(session_block_array{pos_par,2}(1,:), session{s}));
        block = session_block_array{pos_par,2}{2,pos_ses}(1,:);
       
        % LLR_11 first block is where EEG forward model is
        % If folder does not exists or does not contain the head model, don't even try
        if ~exist([root_dir_bs '/data/' participant{p} '/LLR_11' session{s} '_' block{1} '_normal'], 'dir');continue;end
        if ~exist([root_dir_bs '/data/' participant{p} '/LLR_11' session{s} '_' block{1} '_normal/headmodel_surf_duneuro.mat'],'file')
            error([' no EEG forward model in ' participant{p} '_LLR_11' session{s} '_' block{1} '_normal']);
        end


        disp(' ');      
        disp('-------------------------');  
        disp(['Copying EEG forward model from first block to others for ' participant{p} '_' session{s}]);
        disp(datetime)
        disp(' ');

        % The EEG forward model from first block
        file_name_leadfield = [root_dir_bs '/data/' participant{p} '/LLR_11' session{s} '_' block{1} '_normal/headmodel_surf_duneuro.mat'];
        
        % ONLY SECOND BLOCK ONWARDS (FIRST ALREADY HAS IT)
        for b = 2:length(block)

            % If stated, find and delete any previous head model (overwrite it) before copying 
            if delete_previous_file == 1
                % check if there is already a headmodel OF THIS
                % TYPE in the folder where we are going to copy the file
                if isfile([root_dir_bs '/data/' participant{p} '/LLR_11' session{s} '_' block{b} '_normal/headmodel_surf_duneuro.mat'])
                    delete([root_dir_bs '/data/' participant{p} '/LLR_11' session{s} '_' block{b} '_normal/headmodel_surf_duneuro.mat'])
                end
                % head model = headmodel_surf_os_meg.mat
                % source = contains KERNEL endswith _Source_Regul_MEG.mat                 
            end
            if exist([root_dir_bs '/data/' participant{p} '/LLR_11' session{s} '_' block{b} '_normal/'],'dir')
                copyfile(file_name_leadfield,[root_dir_bs '/data/' participant{p} '/LLR_11' session{s} '_' block{b} '_normal/headmodel_surf_duneuro.mat']);
            end
        end
    end
end

% Merge EEG/MEG forward models
for p = 1:length(participant)
        
    % Define sessions within this participant
    pos_par = find(strcmp(session_block_array(:,1), participant{p}(1:4)));
    session = session_block_array{pos_par,2}(1,:);
    
    for s = 1:length(session)
        
        % Define blocks within this session
        pos_ses = find(strcmp(session_block_array{pos_par,2}(1,:), session{s}));
        block = session_block_array{pos_par,2}{2,pos_ses}(1,:);
        
        for b = 1:length(block)
            
            % If folder does not exists or 
            % one of the forward models to merge is not present, go to next block
            if ~exist([root_dir_bs '/data/' participant{p} '/LLR_11' session{s} '_' block{b} '_normal'], 'dir'); continue; end
            if ~exist([root_dir_bs '/data/' participant{p} '/LLR_11' session{s} '_' block{b} '_normal/headmodel_surf_duneuro.mat'],'file'); continue; end
            if ~exist([root_dir_bs '/data/' participant{p} '/LLR_11' session{s} '_' block{b} '_normal/headmodel_surf_os_meg.mat'],'file'); continue; end

            EEG_leadfield = [root_dir_bs '/data/' participant{p} '/LLR_11' session{s} '_' block{b} '_normal/headmodel_surf_duneuro.mat'];
            MEG_leadfield = [root_dir_bs '/data/' participant{p} '/LLR_11' session{s} '_' block{b} '_normal/headmodel_surf_os_meg.mat'];
    
            disp(' ');
            disp('-------------------------');  
            disp(['Merging EEG/MEG leadfields for ' participant{p} '_' session{s} '_' block{b}]);
            disp(datetime)
            disp(' ');   
            
            % If stated, find and delete any previous combined leadfield 
            if delete_previous_file == 1
                if isfile([root_dir_bs '/data/' participant{p} '/LLR_11' session{s} '_' block{b} '_normal/headmodel_surf_duneuro_os_meg.mat'])
                    delete([root_dir_bs '/data/' participant{p} '/LLR_11' session{s} '_' block{b} '_normal/headmodel_surf_duneuro_os_meg.mat'])
                end              
            end
           
            % FUNCTION TO COMBINE LEADFIELDS
            merge_leadfields(EEG_leadfield, MEG_leadfield)
            % Will combine and delete the original ones
            % Outcome will be saved with same name as if it was
            % computed manually: headmodel_surf_duneuro_os_meg.mat
            
        end
    end
end

% Copy combined leadfield to other conditions (up to now it was all in '11')
for p = 1:length(participant)
    
    % Define sessions within this participant
    pos_par = find(strcmp(session_block_array(:,1), participant{p}(1:4)));
    session = session_block_array{pos_par,2}(1,:);
    
    for s = 1:length(session)
        
        % Define blocks within this session
        pos_ses = find(strcmp(session_block_array{pos_par,2}(1,:), session{s}));
        block = session_block_array{pos_par,2}{2,pos_ses}(1,:);
        
        for b = 1:length(block)
            
            % LLR_11 is the folder in which they were created in the step before
            if ~exist([root_dir_bs '/data/' participant{p} '/LLR_11' session{s} '_' block{b} '_normal'], 'dir'); continue;end
           
            if ~exist([root_dir_bs '/data/' participant{p} '/LLR_11' session{s} '_' block{b} '_normal/headmodel_surf_duneuro_os_meg.mat'],'file') 
                error(['no combined forward model for ' participant{p} '_LLR_11' session{s} '_' block{b} '_normal'])
            end
            % If folder exists but does not contain the combined head model, don't even try
            
            disp(' ');      
            disp('-------------------------');  
            disp(['Copying between conditions forward model (EEG + MEG) for ' participant{p} '_' session{s} '_' block{b}]);
            disp(datetime)
            disp(' ');

            % This name is the one we obtain for the MEG overlapping spheres
            file_name_leadfield = [root_dir_bs '/data/' participant{p} '/LLR_11' session{s} '_' block{b} '_normal/headmodel_surf_duneuro_os_meg.mat'];
            
            for c = 2:length(condition)
                if  strcmp(condition{c}, '11') % IF ITS LLR_11 do not copy, just in case
                    continue
                end

                % If stated, find and delete any previous head model (overwrite it) before copying 
                if delete_previous_file == 1
                    % check if there is already a headmodel OF THIS
                    % TYPE in the folder where we are going to copy the file
                    if isfile([root_dir_bs '/data/' participant{p} '/LLR_' condition{c} session{s} '_' block{b} '_normal/headmodel_surf_duneuro_os_meg.mat'])
                        delete([root_dir_bs '/data/' participant{p} '/LLR_' condition{c} session{s} '_' block{b} '_normal/headmodel_surf_duneuro_os_meg.mat'])
                    end
                    % head model = headmodel_surf_os_meg.mat
                    % source = contains KERNEL endswith _Source_Regul_MEG.mat                 
                end
                if exist([root_dir_bs '/data/' participant{p} '/LLR_' condition{c} session{s} '_' block{b} '_normal/'],'dir')
                    copyfile(file_name_leadfield,[root_dir_bs '/data/' participant{p} '/LLR_' condition{c} session{s} '_' block{b} '_normal/headmodel_surf_duneuro_os_meg.mat']);
                end
            end
        end
    end
end

% Reset modality_data variable to its original self
modality_data = {'EEG','MEG','both_mod'};

save([root_dir '/Sources/No_head_model_created.mat'],'No_head_model_created');
clearvars('-except', initialVars{:});
disp 'DONE WITH COMPUTING FORWARD MODEL FOR EACH PARTICIPANT, SESSION AND BLOCK (LLR)!!!'
disp(datetime)
toc

%% Compute inverse solutions 
% Separately for every condition, wave, session block and participant

% IDEA: If data (or condition) really matters to perform min norm,
% why not averaging across all conditions for each block to have a
% waveform that is representative and feed min-norm with that

% IDEA 2: For the EEG, that representative waveform would include data
% averaged across the three blocks

tic
disp(' ');      
disp('-------------------------');
disp('COMPUTING INVERSE SOLUTIONS (LLR)');
disp(datetime)
disp('-------------------------');
disp(' ');

% Here this matters
if source_analysis == 1
    modality_data = {'EEG'};
elseif source_analysis == 2
    modality_data = {'MEG'};
elseif source_analysis == 3
    modality_data = {'both_mod'};
elseif source_analysis == 4
    modality_data = {'EEG','MEG','both_mod'};
end

No_inverse_solution_created = {};
count_no_sour = 1;

for p = 1:length(participant)
    
    disp(' ');      
    disp('-------------------------');
    disp(['loading participant ' participant{p}]);
    disp(datetime)
    disp(' '); 
    
    prot_subs = bst_get('ProtocolSubjects');
    current_sub = find(strcmp({prot_subs.Subject.Name}, participant{p}));
    db_reload_conditions(current_sub);
    
    % Define sessions within this participant
    pos_par = find(strcmp(session_block_array(:,1), participant{p}(1:4)));
    session = session_block_array{pos_par,2}(1,:);
    
    for s = 1:length(session)
        
        % Define blocks within this session
        pos_ses = find(strcmp(session_block_array{pos_par,2}(1,:), session{s}));
        block = session_block_array{pos_par,2}{2,pos_ses}(1,:);
        
        for b = 1:length(block)
            for w = 1:length(wave) %  A WAY TO CHANGE IT EASILY TO MLR
                for mode = 1:length(modality_data)
                    for c = 1:length(condition)
                        if ~exist([root_dir_bs '/data/' participant{p} '/' wave{w} '_' condition{c} session{s} '_' block{b} '_normal'], 'dir')
                        continue
                        end
                    
                        folders = dir([root_dir_bs '/data/' participant{p} '/' wave{w} '_' condition{c} session{s} '_' block{b} '_normal']);            
                        results =  contains({folders.name},[wave{w} '_' condition{c} session{s} '_' block{b} '_normal']) & contains({folders.name}, [modality_data{mode} '_average_source']);
                        infolder = find(results);
                    
                        % If both or either the average or the noise cov file
                        % are not in the folder, keep track and go to next condition
                        if (~exist([root_dir_bs '/data/' participant{p} '/' wave{w} '_' condition{c} session{s} '_' block{b} '_normal/noisecov_full.mat'], 'file')) && (isempty(infolder))
                            No_inverse_solution_created{count_no_sour,1} = [participant{p} '/' wave{w} '_' condition{c} session{s} '_' block{b} '_normal']; %#ok<*SAGROW>
                            No_inverse_solution_created{count_no_sour,2} = 'nothing, which is sad';
                            count_no_sour = count_no_sour +1;
                            continue
                        elseif isempty(infolder)            
                            No_inverse_solution_created{count_no_sour,1} = [participant{p} '/' wave{w} '_' condition{c} session{s} '_' block{b} '_normal']; %#ok<*SAGROW>
                            No_inverse_solution_created{count_no_sour,2} = 'no average file';
                            count_no_sour = count_no_sour +1;
                            continue
                        elseif ~exist([root_dir_bs '/data/' participant{p} '/' wave{w} '_' condition{c} session{s} '_' block{b} '_normal/noisecov_full.mat'], 'file')
                            No_inverse_solution_created{count_no_sour,1} = [participant{p} '/' wave{w} '_' condition{c} session{s} '_' block{b} '_normal']; %#ok<*SAGROW>
                            No_inverse_solution_created{count_no_sour,2} = 'no covariance matrix';
                            count_no_sour = count_no_sour +1;
                            continue
                        end
                    
                        % Rare case in which more than one average file meeting
                        % the criteria exists
                        if length(infolder) > 1; infolder = infolder(1); end
                    
                        file_name = [participant{p} '/' wave{w} '_' condition{c} session{s} '_' block{b} '_normal/' folders(infolder).name];
                        sFiles = file_name;
                    
                        if isempty(sFiles)
                            error('for whatever reason, sFiles is empty, my friend...')
                        end
                                       
                        % TAKE INTO ACCOUNT
                        % 1) THIS WILL BE DIFFERENT FOR MLR, USING VOLUME
                        % along with forward models using volume too
                        % 2) We discussed with Brian at one point that
                        % dSPM option was not to be used, but current density option
                        % intead, then average, then do normalization like in baseline (z
                        % score) at the subject average. But we saw how results
                        % did not make any sense and decided to switch back to
                        % dSPM
                    
                        if strcmp(modality_data{mode},'EEG') % if EEG source solution alone was requested
                        
                            disp(' ');
                            disp('-------------------------');  
                            disp(['Computing inverse solutions (EEG) for ' participant{p} '_' wave{w} '_' condition{c} session{s} '_' block{b}]);
                            disp(datetime)
                            disp(' ');
                        
                            % IMPORTANT: keep the original sFiles variable alive within this loop
                            for nopt = 1:length(source_noise_options)

                                % If stated, find and delete any previous source kernel (overwrite it) before copying 
                                if delete_previous_file == 1
                                    % check if there is already a source of THIS
                                    % TYPE in the folder where we are gonna calculate it
                                    folders_delete = dir([root_dir_bs '/data/' participant{p} '/' wave{w} '_' condition{c} session{s} '_' block{b} '_normal']);
                                    results_delete =  contains({folders_delete.name},'KERNEL') & endsWith({folders_delete.name}, [source_noise_tags{nopt} '_' modality_data{mode} '.mat']); 
                                    infolder_delete = find(results_delete);
                                    if ~isempty(infolder_delete) % file exists, therefore delete it
                                       delete([root_dir_bs '/data/' participant{p} '/' wave{w} '_' condition{c} session{s} '_' block{b} '_normal/' folders_delete(infolder_delete).name]);
                                    end
                                    % source = contains KERNEL endswith _Source_Regul_MEG.mat                 
                                end

                                % Process: Compute sources [2018]
                                sFiles_source = bst_process('CallProcess', 'process_inverse_2018', sFiles, [], ...
                                    'output',  1, ...  % Kernel only: shared
                                    'inverse', struct(...
                                         'Comment',        ['dSPM-unscaled: ' modality_data{mode} ' ALL'], ...
                                         'InverseMethod',  'minnorm', ...
                                         'InverseMeasure', 'dspm2018', ...
                                         'SourceOrient',   {{'loose'}}, ...
                                         'Loose',          0.4, ...
                                         'UseDepth',       1, ...
                                         'WeightExp',      0.5, ...
                                         'WeightLimit',    10, ...
                                         'NoiseMethod',    source_noise_options{nopt}, ...
                                         'NoiseReg',       0.1, ...
                                         'SnrMethod',      'fixed', ...
                                         'SnrRms',         1e-06, ...
                                         'SnrFixed',       3, ...
                                         'ComputeKernel',  1, ...
                                         'DataTypes',      {{'EEG'}}));

                                        % Process: Change name in folder (for code)
                                        sFiles_name = bst_process('CallProcess', 'process_add_tag', sFiles_source, [], ...
                                            'tag',           [source_noise_tags{nopt} '_' modality_data{mode}], ...
                                            'output',        2);  % Add to file path

                                        % Process: Set name to see in GUI
                                        sFiles_name = bst_process('CallProcess', 'process_set_comment', sFiles_name, [], ...
                                            'tag',           [source_noise_tags{nopt} '_' modality_data{mode}], ...
                                            'isindex',       0);
                                        clear sFiles_source sFiles_name

                            end
                        end
                    
                        if strcmp(modality_data{mode}, 'MEG')
                        
                            disp(' ');
                            disp('-------------------------');  
                            disp(['Computing inverse solutions (MEG) for ' participant{p} '_' wave{w} '_' condition{c} session{s} '_' block{b}]);
                            disp(datetime)
                            disp(' ');

                            % IMPORTANT: keep the original sFiles variable alive within this loop
                            for nopt = 1:length(source_noise_options)

                                % If stated, find and delete any previous source kernel (overwrite it) before copying 
                                if delete_previous_file == 1
                                    % check if there is already a source of THIS
                                    % TYPE in the folder where we are gonna calculate it
                                    folders_delete = dir([root_dir_bs '/data/' participant{p} '/' wave{w} '_' condition{c} session{s} '_' block{b} '_normal']);
                                    results_delete =  contains({folders_delete.name},'KERNEL') & endsWith({folders_delete.name}, [source_noise_tags{nopt} '_' modality_data{mode} '.mat']); 
                                    infolder_delete = find(results_delete);
                                    if ~isempty(infolder_delete) % file exists, therefore delete it
                                       delete([root_dir_bs '/data/' participant{p} '/' wave{w} '_' condition{c} session{s} '_' block{b} '_normal/' folders_delete(infolder_delete).name]);
                                    end
                                    % source = contains KERNEL endswith _Source_Regul_MEG.mat                 
                                end

                                % Process: Compute sources [2018]
                                sFiles_source = bst_process('CallProcess', 'process_inverse_2018', sFiles, [], ...
                                    'output',  1, ...  % Kernel only: shared
                                    'inverse', struct(...
                                         'Comment',        ['dSPM-unscaled: ' modality_data{mode} ' ALL'], ...
                                         'InverseMethod',  'minnorm', ...
                                         'InverseMeasure', 'dspm2018', ...
                                         'SourceOrient',   {{'loose'}}, ...
                                         'Loose',          0.4, ...
                                         'UseDepth',       1, ...
                                         'WeightExp',      0.5, ...
                                         'WeightLimit',    10, ...
                                         'NoiseMethod',    source_noise_options{nopt}, ...
                                         'NoiseReg',       0.1, ...
                                         'SnrMethod',      'fixed', ...
                                         'SnrRms',         1e-06, ...
                                         'SnrFixed',       3, ...
                                         'ComputeKernel',  1, ...
                                         'DataTypes',      {{'MEG GRAD', 'MEG MAG'}}));

                                        % Process: Change name in folder (for code)
                                        sFiles_name = bst_process('CallProcess', 'process_add_tag', sFiles_source, [], ...
                                            'tag',           [source_noise_tags{nopt} '_' modality_data{mode}], ...
                                            'output',        2);  % Add to file path

                                        % Process: Set name to see in GUI
                                        sFiles_name = bst_process('CallProcess', 'process_set_comment', sFiles_name, [], ...
                                            'tag',           [source_noise_tags{nopt} '_' modality_data{mode}], ...
                                            'isindex',       0);
                                        clear sFiles_source sFiles_name

                            end
                        
                        % {{'MEG GRAD', 'MEG MAG', 'EEG'}}
                    
                        end
                    
                        if strcmp(modality_data{mode}, 'both_mod')
                        
                            disp(' ');
                            disp('-------------------------');  
                            disp(['Computing inverse solutions (MEG) for ' participant{p} '_' wave{w} '_' condition{c} session{s} '_' block{b}]);
                            disp(datetime)
                            disp(' ');

                            % IMPORTANT: keep the original sFiles variable alive within this loop
                            for nopt = 1:length(source_noise_options)

                                % If stated, find and delete any previous source kernel (overwrite it) before copying 
                                if delete_previous_file == 1
                                    % check if there is already a source of THIS
                                    % TYPE in the folder where we are gonna calculate it
                                    folders_delete = dir([root_dir_bs '/data/' participant{p} '/' wave{w} '_' condition{c} session{s} '_' block{b} '_normal']);
                                    results_delete =  contains({folders_delete.name},'KERNEL') & endsWith({folders_delete.name}, [source_noise_tags{nopt} '_' modality_data{mode} '.mat']); 
                                    infolder_delete = find(results_delete);
                                    if ~isempty(infolder_delete) % file exists, therefore delete it
                                       delete([root_dir_bs '/data/' participant{p} '/' wave{w} '_' condition{c} session{s} '_' block{b} '_normal/' folders_delete(infolder_delete).name]);
                                    end
                                    % source = contains KERNEL endswith _Source_Regul_MEG.mat                 
                                end

                                % Process: Compute sources [2018]
                                sFiles_source = bst_process('CallProcess', 'process_inverse_2018', sFiles, [], ...
                                    'output',  1, ...  % Kernel only: shared
                                    'inverse', struct(...
                                         'Comment',        ['dSPM-unscaled: ' modality_data{mode} ' ALL'], ...
                                         'InverseMethod',  'minnorm', ...
                                         'InverseMeasure', 'dspm2018', ...
                                         'SourceOrient',   {{'loose'}}, ...
                                         'Loose',          0.4, ...
                                         'UseDepth',       1, ...
                                         'WeightExp',      0.5, ...
                                         'WeightLimit',    10, ...
                                         'NoiseMethod',    source_noise_options{nopt}, ...
                                         'NoiseReg',       0.1, ...
                                         'SnrMethod',      'fixed', ...
                                         'SnrRms',         1e-06, ...
                                         'SnrFixed',       3, ...
                                         'ComputeKernel',  1, ...
                                         'DataTypes',      {{'MEG GRAD', 'MEG MAG', 'EEG'}}));

                                        % Process: Change name in folder (for code)
                                        sFiles_name = bst_process('CallProcess', 'process_add_tag', sFiles_source, [], ...
                                            'tag',           [source_noise_tags{nopt} '_' modality_data{mode}], ...
                                            'output',        2);  % Add to file path

                                        % Process: Set name to see in GUI
                                        sFiles_name = bst_process('CallProcess', 'process_set_comment', sFiles_name, [], ...
                                            'tag',           [source_noise_tags{nopt} '_' modality_data{mode}], ...
                                            'isindex',       0);
                                        clear sFiles_source sFiles_name

                            end
                        end
                    end
                end
            end
        end
    end
end

% Reset modality_data variable to its original self
modality_data = {'EEG','MEG','both_mod'};

save([root_dir '/Sources/No_inverse_solution_created.mat'],'No_inverse_solution_created');
clearvars('-except', initialVars{:});
disp 'DONE COMPUTING INVERSE SOLUTIONS (LLR)!!!'
disp(datetime)
toc

%% Average inverse solutions across blocks from all sessions

tic
disp(' ');
disp('-------------------------');  
disp('AVERAGING INVERSE SOLUTIONS ACROSS BLOCKS FROM ALL SESSIONS (LLR)'); 
disp(datetime)
disp('-------------------------');
disp(' ');

% Which ones will be averaged
if source_analysis == 1
    modality_data = {'EEG'};
elseif source_analysis == 2
    modality_data = {'MEG'};
elseif source_analysis == 3
    modality_data = {'both_mod'};
elseif source_analysis == 4
    modality_data = {'EEG','MEG','both_mod'};
end

for p = 1:length(participant)

    disp(' ');      
    disp('-------------------------');
    disp(['loading participant ' participant{p}]);
    disp(datetime)
    disp(' '); 
    
    prot_subs = bst_get('ProtocolSubjects');
    current_sub = find(strcmp({prot_subs.Subject.Name}, participant{p}));
    db_reload_conditions(current_sub);
    
    % Define sessions within this participant
    pos_par = find(strcmp(session_block_array(:,1), participant{p}(1:4)));
    session = session_block_array{pos_par,2}(1,:);
    
    for w = 1:length(wave)
    for mode = 1:length(modality_data)
        for nopt = 1:length(source_noise_tags)
        for c = 1:length(condition)
            file_names = {}; % create one for each p, c and w, so that it includes all sessions and blocks
            list_order = 1;
            for s = 1:length(session)
                
                % Define blocks within this session (you are not gonna loop
                % through them, but you will discard any not present in the list for this session)
                pos_ses = find(strcmp(session_block_array{pos_par,2}(1,:), session{s}));
                block = session_block_array{pos_par,2}{2,pos_ses}(1,:);
                
                folders = dir([root_dir_bs '/data/' participant{p} '/']); 
                results =  contains({folders.name},[wave{w} '_' condition{c} session{s}]) & contains({folders.name},'_normal');
                infolder = find(results); % So get all blocks you find for that session
                if isempty(infolder)
                    error(['no ' wave{w} '_' condition{c} session{s} '_normal blocks for ' participant{p}]);
                end
                
                    for l = 1:length(infolder) % for each block found
                        line = infolder(l); % first position index and so on
                        % Check if that block is in the block list, if not, go to next block
                        init_block_name = regexp(folders(line).name, '_b', 'once');
                        fol_nam = contains(block, folders(line).name(init_block_name+1:init_block_name+2));
                        in_fol_nam = find(fol_nam); % index of position in logical
                        if isempty(in_fol_nam); continue; end % the block in the folder name is not present in the list, so go to next
                        % End checking if block exist in list
                        sub_folders = dir([root_dir_bs '/data/' participant{p} '/' folders(line).name '/']);
                        sub_results_kernel =  contains({sub_folders.name},[source_noise_tags{nopt} '_' modality_data{mode}]);
                        sub_infolder_kernel = find(sub_results_kernel);
                        sub_results_average =  contains({sub_folders.name},[modality_data{mode} '_average_source']);
                        sub_infolder_average = find(sub_results_average);
                        % If the folder exists but not the Kernel or the average
                        if isempty(sub_infolder_kernel) || isempty(sub_infolder_average); continue; end
                        % If there are more than one kernels or averages (shouln't be) just pick the first
                        if length(sub_infolder_average) > 1; sub_infolder_average = sub_infolder_average(1); end
                        if length(sub_infolder_kernel) > 1; sub_infolder_kernel = sub_infolder_kernel(1); end

                        % Name the links, not the kernel or the avg file separately
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        file_names{list_order} = ['link|' participant{p} '/' folders(line).name '/' ...
                            sub_folders(sub_infolder_kernel).name '|' ...
                            participant{p} '/' folders(line).name '/' sub_folders(sub_infolder_average).name];
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                        list_order = list_order + 1;
                    end
                
            end
            file_names = file_names(~cellfun('isempty', file_names')); % to avoid empty cells
            sFiles = file_names; % should contain all blocks from all sessions
                
            if isempty(sFiles)
            error('for whatever reason, sFiles is empty, my friend...')
            end
            
            if ~isempty(sFiles) % in case 'continue' option was selected for all
                
               % If stated, find and delete any previous source in intra folder (overwrite it)
                if delete_previous_file == 1
                    % check if there is already a source of THIS
                    % TYPE in the folder where we are gonna calculate it
                    folders_delete = dir([root_dir_bs '/data/' participant{p} '/@intra/']);
                    results_delete =  endsWith({folders_delete.name}, [modality_data{mode} '_' source_noise_tags{nopt} '_' wave{w} '_' condition{c} '.mat']); 
                    infolder_delete = find(results_delete);
                    if ~isempty(infolder_delete) % file exists, therefore delete it
                       delete([root_dir_bs '/data/' participant{p} '/@intra/' folders_delete(infolder_delete).name]);
                    end
                    % source = contains KERNEL endswith _Source_Regul_MEG.mat                 
                end
            
                disp(' ');
                disp('-------------------------');  
                disp(['Average blocks from all sessions (' modality_data{mode} '_' source_noise_tags{nopt} ') for ' participant{p} '_' wave{w} '_' condition{c}]);
                disp(datetime)
                disp(' ');

                % WILL HAVE TO AVERAGE ALL OTHER SOURCE SOLUTIONS THAT MAY
                % HAVE BEEN CREATED. FOR NOW, ONLY Overlapping spheres & Cortical surface
                
                % Average using default function (because we are using vertices, not channels)
                sFiles = bst_process('CallProcess', 'process_average', sFiles, [], ...
                    'avgtype',         1, ...  % Everything
                    'avg_func',        1, ...  % Arithmetic average:  mean(x)
                    'weighted',        1, ...
                    'scalenormalized', 0);
                
                % USE OPTION TO NORMALIZE HERE???

                % Process: Add tag:
                sFiles = bst_process('CallProcess', 'process_add_tag', sFiles, [], ...
                    'tag',           [modality_data{mode} '_' source_noise_tags{nopt} '_' wave{w} '_' condition{c}], ...
                    'output',        2);  % Add to file name

                % Process: Set name
                sFiles = bst_process('CallProcess', 'process_set_comment', sFiles, [], ...
                    'tag',           [modality_data{mode} '_' source_noise_tags{nopt} '_' wave{w} '_' condition{c}], ...
                    'isindex',       1); 
            end
        end
        end
    end
    end
end

% Reset modality_data variable to its original self
modality_data = {'EEG','MEG','both_mod'};

clearvars('-except', initialVars{:});
disp 'DONE AVERAGING INVERSE SOLUTIONS ACROSS BLOCKS FROM ALL SESSIONS (LLR)!!!'
disp(datetime)
toc

%% (If chosen) Average inverse solutions across blocks (source representation of each session)

if sources_single_session == 1

tic
disp(' ');      
disp('-------------------------');  
disp('AVERAGING INVERSE SOLUTIONS ACROSS BLOCKS (LLR)'); 
disp(datetime)
disp('-------------------------');
disp(' ');

% Which ones will be averaged
if source_analysis == 1
    modality_data = {'EEG'};
elseif source_analysis == 2
    modality_data = {'MEG'};
elseif source_analysis == 3
    modality_data = {'both_mod'};
elseif source_analysis == 4
    modality_data = {'EEG','MEG','both_mod'};
end

for p = 1:length(participant)
    
    disp(' ');      
    disp('-------------------------');
    disp(['loading participant ' participant{p}]);
    disp(datetime)
    disp(' '); 
    
    prot_subs = bst_get('ProtocolSubjects');
    current_sub = find(strcmp({prot_subs.Subject.Name}, participant{p}));
    db_reload_conditions(current_sub);
    
    % Define sessions within this participant
    pos_par = find(strcmp(session_block_array(:,1), participant{p}(1:4)));
    session = session_block_array{pos_par,2}(1,:);
    
    for s = 1:length(session)
        
        % Define blocks within this session
        pos_ses = find(strcmp(session_block_array{pos_par,2}(1,:), session{s}));
        block = session_block_array{pos_par,2}{2,pos_ses}(1,:);
        
        for w = 1:length(wave)
        for mode = 1:length(modality_data)
            for c = 1:length(condition)
                folders = dir([root_dir_bs '/data/' participant{p} '/']);
                results =  contains({folders.name},[wave{w} '_' condition{c} session{s}]) & contains({folders.name},'_normal');
                infolder = find(results); % so get all the blocks form that session
                if isempty(infolder) % 
                    error(['no ' wave{w} '_' condition{c} session{s} '_normal blocks for ' participant{p}]);
                end
                for nopt = 1:length(source_noise_tags)
                    for l = 1:length(infolder) % for every block found, find kernel and averages names
                        line = infolder(l); % first position index and so on
                        % Check if that block is in the block list, if not, go to next block
                        init_block_name = regexp(folders(line).name, '_b', 'once');
                        fol_nam = contains(block, folders(line).name(init_block_name+1:init_block_name+2));
                        in_fol_nam = find(fol_nam); % index of position in logical
                        if isempty(in_fol_nam); continue; end % the block in the folder name is not present in the list, so go to next
                        % End checking if block exist in list
                        sub_folders = dir([root_dir_bs '/data/' participant{p} '/' folders(line).name '/']);
                        sub_results_kernel =  contains({sub_folders.name},[source_noise_tags{nopt} '_' modality_data{mode}]);
                        sub_infolder_kernel = find(sub_results_kernel);
                        sub_results_average =  contains({sub_folders.name},[modality_data{mode} '_average_source']);
                        sub_infolder_average = find(sub_results_average);
                        % If the folder exists but not the Kernel or the average
                        if isempty(sub_infolder_kernel) || isempty(sub_infolder_average); continue; end
                        % If there are more than one kernels or averages (shouln't be) just pick the first
                        if length(sub_infolder_average) > 1; sub_infolder_average = sub_infolder_average(1); end
                        if length(sub_infolder_kernel) > 1; sub_infolder_kernel = sub_infolder_kernel(1); end

                        % Name the links, not the kernel or the avg file separately
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        file_names{l} = ['link|' participant{p} '/' folders(line).name '/' ...
                            sub_folders(sub_infolder_kernel).name '|' ...
                            participant{p} '/' folders(line).name '/' sub_folders(sub_infolder_average).name];
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                    end
                    % if all three blocks are empty, file_names variable
                    % won't exist, so go to next condition/wave/session
                    if ~exist('file_names','var'); continue;end
                    file_names = file_names(~cellfun('isempty', file_names')); % to avoid empty cells
                    sFiles = file_names;
                    
                    if isempty(sFiles)
                    error('for whatever reason, sFiles is empty, my friend...')
                    end

                    disp(' ');
                    disp('-------------------------');  
                    disp(['Average sources (' source_noise_tags{nopt} ') across blocks for ' participant{p} '_' modality_data{mode} '_' wave{w} '_' condition{c} session{s}]);
                    disp(datetime)
                    disp(' ');

                    % WILL HAVE TO AVERAGE ALL OTHER SOURCE SOLUTIONS THAT MAY
                    % HAVE BEEN CREATED. FOR NOW, ONLY Overlapping spheres & Cortical surface

                    % If stated, find and delete any previous source in intra folder (overwrite it)
                    if delete_previous_file == 1
                        % check if there is already a source of THIS
                        % TYPE in the folder where we are gonna calculate it
                        folders_delete = dir([root_dir_bs '/data/' participant{p} '/@intra/']);
                        results_delete =  endsWith({folders_delete.name}, [modality_data{mode} '_' source_noise_tags{nopt} '_' wave{w} '_' condition{c} session{s} '.mat']); 
                        infolder_delete = find(results_delete);
                        if ~isempty(infolder_delete) % file exists, therefore delete it
                           delete([root_dir_bs '/data/' participant{p} '/@intra/' folders_delete(infolder_delete).name]);
                        end
                        % source = contains KERNEL endswith _Source_Regul_MEG.mat                 
                    end
                    
                    % Process: Weighted Average: Everything (Because we are still at block level)
                    sFiles = bst_process('CallProcess', 'process_average', sFiles, [], ...
                        'avgtype',         1, ...  % Everything
                        'avg_func',        1, ...  % Arithmetic average:  mean(x)
                        'weighted',        1, ...
                        'scalenormalized', 0);
                    
                    % USE OPTION TO NORMALIZE HERE???

                    % Process: Add tag:
                    sFiles = bst_process('CallProcess', 'process_add_tag', sFiles, [], ...
                        'tag',           [modality_data{mode} '_' source_noise_tags{nopt} '_' wave{w} '_' condition{c} session{s}], ...
                        'output',        2);  % Add to file name

                    % Process: Set name
                    sFiles = bst_process('CallProcess', 'process_set_comment', sFiles, [], ...
                        'tag',           [modality_data{mode} '_' source_noise_tags{nopt} '_' wave{w} '_' condition{c} session{s}], ...
                        'isindex',       1);
                end
            end
        end
        end
    end
end

% Reset modality_data variable to its original self
modality_data = {'EEG','MEG','both_mod'};

clearvars('-except', initialVars{:});
disp 'DONE WITH AVERAGING INVERSE SOLUTIONS ACROSS BLOCKS (LLR)!!!'
disp(datetime)
toc

end

%% (If chosen) Average inverse solutions across ALL conditions per each session

if sources_single_session == 1

tic
disp(' ');      
disp('-------------------------');  
disp('OBTAINING ALL_COND SOURCE AVERAGES PER SESSION'); 
disp(datetime)
disp('-------------------------');
disp(' ');

% Which ones will be averaged
if source_analysis == 1
    modality_data = {'EEG'};
elseif source_analysis == 2
    modality_data = {'MEG'};
elseif source_analysis == 3
    modality_data = {'both_mod'};
elseif source_analysis == 4
    modality_data = {'EEG','MEG','both_mod'};
end

for p = 1:length(participant)
    
    disp(' ');      
    disp('-------------------------');
    disp(['loading participant ' participant{p}]);
    disp(datetime)
    disp(' '); 

    % Reload subject normally
    prot_subs = bst_get('ProtocolSubjects');
    current_sub = find(strcmp({prot_subs.Subject.Name}, participant{p}));
    db_reload_conditions(current_sub);
    
    % Define sessions within this participant
    pos_par = find(strcmp(session_block_array(:,1), participant{p}(1:4)));
    session = session_block_array{pos_par,2}(1,:);
    for mode = 1:length(modality_data)
        for w = 1:length(wave)        
            for s = 1:length(session)
                
                % Define blocks within this session
                pos_ses = find(strcmp(session_block_array{pos_par,2}(1,:), session{s}));
                block = session_block_array{pos_par,2}{2,pos_ses}(1,:);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ALL %%%%%%%%%%%%%%%%%%%%%%%%
                folders = dir([root_dir_bs '/data/' participant{p} '/']);
                results =  (contains({folders.name},[wave{w} '_11' session{s}])... 
                | contains({folders.name},[wave{w} '_12' session{s}]) ...
                | contains({folders.name},[wave{w} '_13' session{s}]) ...
                | contains({folders.name},[wave{w} '_31' session{s}]) ...
                | contains({folders.name},[wave{w} '_32' session{s}]) ...
                | contains({folders.name},[wave{w} '_33' session{s}]) ...
                | contains({folders.name},[wave{w} '_51' session{s}]) ...
                | contains({folders.name},[wave{w} '_52' session{s}]) ...
                | contains({folders.name},[wave{w} '_53' session{s}]) ...
                | contains({folders.name},[wave{w} '_71' session{s}]) ...
                | contains({folders.name},[wave{w} '_72' session{s}]) ...
                | contains({folders.name},[wave{w} '_73' session{s}]) ...
                | contains({folders.name},[wave{w} '_91' session{s}]) ...
                | contains({folders.name},[wave{w} '_92' session{s}]) ...
                | contains({folders.name},[wave{w} '_93' session{s}])) ...
                & contains({folders.name},'_normal');   
                infolder = find(results); % So get all blocks you find for that session and gather all conditions
                if isempty(infolder)
                    error('no coincidences');
                end            

                for nopt = 1:length(source_noise_tags)
                    for l = 1:length(infolder) % for every block found, find kernel and averages names
                        line = infolder(l); % first position index and so on
                        % Check if that block is in the block list, if not, go to next block
                        init_block_name = regexp(folders(line).name, '_b', 'once');
                        fol_nam = contains(block, folders(line).name(init_block_name+1:init_block_name+2));
                        in_fol_nam = find(fol_nam); % index of position in logical
                        if isempty(in_fol_nam); continue; end % the block in the folder name is not present in the list, so go to next
                        % End checking if block exist in list
                        sub_folders = dir([root_dir_bs '/data/' participant{p} '/' folders(line).name '/']);
                        sub_results_kernel =  contains({sub_folders.name},[source_noise_tags{nopt} '_' modality_data{mode}]);
                        sub_infolder_kernel = find(sub_results_kernel);
                        sub_results_average =  contains({sub_folders.name},[modality_data{mode} '_average_source']);
                        sub_infolder_average = find(sub_results_average);
                        % If the folder exists but not the Kernel or the average
                        if isempty(sub_infolder_kernel) || isempty(sub_infolder_average); continue; end
                        % If there are more than one kernels or averages (shouln't be) just pick the first
                        if length(sub_infolder_average) > 1; sub_infolder_average = sub_infolder_average(1); end
                        if length(sub_infolder_kernel) > 1; sub_infolder_kernel = sub_infolder_kernel(1); end

                        % Name the links, not the kernel or the avg file separately
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        file_names{l} = ['link|' participant{p} '/' folders(line).name '/' ...
                            sub_folders(sub_infolder_kernel).name '|' ...
                            participant{p} '/' folders(line).name '/' sub_folders(sub_infolder_average).name];
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                    end
                    % if all three blocks are empty, file_names variable
                    % won't exist, so go to next condition/wave/session
                    if ~exist('file_names','var'); continue;end
                    file_names = file_names(~cellfun('isempty', file_names')); % to avoid empty cells
                    sFiles = file_names;
                    
                    if isempty(sFiles)
                    error('for whatever reason, sFiles is empty, my friend...')
                    end

                    disp(' ');
                    disp('-------------------------');  
                    disp(['Average sources (' source_noise_tags{nopt} ') across all blocks and conditions for ' participant{p} '_' modality_data{mode} '_' wave{w} '_' session{s}]);
                    disp(datetime)
                    disp(' ');

                    % WILL HAVE TO AVERAGE ALL OTHER SOURCE SOLUTIONS THAT MAY
                    % HAVE BEEN CREATED. FOR NOW, ONLY Overlapping spheres & Cortical surface

                    % If stated, find and delete any previous source in intra folder (overwrite it)
                    if delete_previous_file == 1
                        % Check if there is already a source of THIS
                        % TYPE in the folder where we are gonna calculate it
                        folders_delete = dir([root_dir_bs '/data/' participant{p} '/@intra/']);
                        results_delete = endsWith({folders_delete.name}, [modality_data{mode} '_' source_noise_tags{nopt} '_' wave{w} '_ALL_' session{s} '.mat']); 
                        infolder_delete = find(results_delete);
                        if ~isempty(infolder_delete) % file exists, therefore delete it
                           delete([root_dir_bs '/data/' participant{p} '/@intra/' folders_delete(infolder_delete).name]);
                        end
                        % source = contains KERNEL endswith _Source_Regul_MEG.mat                 
                    end
                    
                    % Process: Weighted Average: Everything (Because we are still at block level)
                    sFiles = bst_process('CallProcess', 'process_average', sFiles, [], ...
                        'avgtype',         1, ...  % Everything
                        'avg_func',        1, ...  % Arithmetic average:  mean(x)
                        'weighted',        1, ...
                        'scalenormalized', 0);
                    
                    % USE OPTION TO NORMALIZE HERE???

                    % Process: Add tag:
                    sFiles = bst_process('CallProcess', 'process_add_tag', sFiles, [], ...
                        'tag',           [modality_data{mode} '_' source_noise_tags{nopt} '_' wave{w} '_ALL_' session{s}], ...
                        'output',        2);  % Add to file name

                    % Process: Set name
                    sFiles = bst_process('CallProcess', 'process_set_comment', sFiles, [], ...
                        'tag',           [modality_data{mode} '_' source_noise_tags{nopt} '_' wave{w} '_ALL_' session{s}], ...
                        'isindex',       1);
                end
            end
        end
    end
end

% reset modality_data variable to its original
modality_data = {'EEG','MEG','both_mod'};

clearvars('-except', initialVars{:});
disp 'DONE OBTAINING ALL_COND SOURCE AVERAGES PER SESSION!!!'
disp(datetime)
toc

end

%% Average SOURCES across dB and ISI within subjects

tic
disp(' ');      
disp('-------------------------');  
disp('AVERAGING dB AND ISI SOURCES WITHIN SUBJECTS (LLR)'); 
disp(datetime)
disp('-------------------------');     
disp(' '); 

% Which ones will be averaged
if source_analysis == 1
    modality_data = {'EEG'};
elseif source_analysis == 2
    modality_data = {'MEG'};
elseif source_analysis == 3
    modality_data = {'both_mod'};
elseif source_analysis == 4
    modality_data = {'EEG','MEG','both_mod'};
end

for p = 1:length(participant)

    disp(' ');      
    disp('-------------------------');
    disp(['loading participant ' participant{p}]);
    disp(datetime)
    disp(' '); 
    
    % Reload only intra folder from subjects, as this is the one used
    [~, iStudies] = bst_get('StudyWithCondition', [participant{p} filesep '@intra']);
    db_reload_studies(iStudies);
    
    for w = 1:length(wave)
    for mode = 1:length(modality_data)
         for nopt = 1:length(source_noise_tags)
            folders = dir([root_dir_bs '/data/' participant{p} '/@intra/']);
            % Quietest
            results =  contains({folders.name},[modality_data{mode} '_' source_noise_tags{nopt} '_' wave{w}]) & ... 
                (endsWith({folders.name},'11.mat') | endsWith({folders.name},'31.mat') | endsWith({folders.name},'51.mat') | ...
                endsWith({folders.name},'71.mat') | endsWith({folders.name},'91.mat'));
            infolder = find(results);
            Quietest_names = {};
            for l= 1:length(infolder)
                line = infolder(l);
                Quietest_names{l} = [participant{p} '/@intra/' folders(line).name];
            end
            % Medium intensity
            results =  contains({folders.name},[modality_data{mode} '_' source_noise_tags{nopt} '_' wave{w}]) & ... 
                (endsWith({folders.name},'12.mat') | endsWith({folders.name},'32.mat') | endsWith({folders.name},'52.mat') | ...
                endsWith({folders.name},'72.mat') | endsWith({folders.name},'92.mat'));
            infolder = find(results);
            Medium_dB_names = {};
            for l= 1:length(infolder)
                line = infolder(l);
                Medium_dB_names{l} = [participant{p} '/@intra/' folders(line).name];
            end
            % Loudest
            results =  contains({folders.name},[modality_data{mode} '_' source_noise_tags{nopt} '_' wave{w}]) & ... 
                (endsWith({folders.name},'13.mat') | endsWith({folders.name},'33.mat') | endsWith({folders.name},'53.mat') | ...
                endsWith({folders.name},'73.mat') | endsWith({folders.name},'93.mat'));
            infolder = find(results);
            Loudest_names = {};
            for l= 1:length(infolder)
                line = infolder(l);
                Loudest_names{l} = [participant{p} '/@intra/' folders(line).name];
            end

            % Fastest
            results =  contains({folders.name},[modality_data{mode} '_' source_noise_tags{nopt} '_' wave{w}]) & ... 
                (endsWith({folders.name},'11.mat') | endsWith({folders.name},'12.mat') | ...
                endsWith({folders.name},'13.mat'));
            infolder = find(results);
            Fastest_names = {};
            for l= 1:length(infolder)
                line = infolder(l);
                Fastest_names{l} = [participant{p} '/@intra/' folders(line).name];
            end
            % Fast
            results =  contains({folders.name},[modality_data{mode} '_' source_noise_tags{nopt} '_' wave{w}]) & ... 
                (endsWith({folders.name},'31.mat') | endsWith({folders.name},'32.mat') | ...
                endsWith({folders.name},'33.mat'));
            infolder = find(results);
            Fast_names = {};
            for l= 1:length(infolder)
                line = infolder(l);
                Fast_names{l} = [participant{p} '/@intra/' folders(line).name];
            end
            % Medium ISI
            results =  contains({folders.name},[modality_data{mode} '_' source_noise_tags{nopt} '_' wave{w}]) & ... 
                (endsWith({folders.name},'51.mat') | endsWith({folders.name},'52.mat') | ...
                endsWith({folders.name},'53.mat'));
            infolder = find(results);
            Medium_ISI_names = {};
            for l= 1:length(infolder)
                line = infolder(l);
                Medium_ISI_names{l} = [participant{p} '/@intra/' folders(line).name];
            end
            % Slow
            results =  contains({folders.name},[modality_data{mode} '_' source_noise_tags{nopt} '_' wave{w}]) & ... 
                (endsWith({folders.name},'71.mat') | endsWith({folders.name},'72.mat') | ...
                endsWith({folders.name},'73.mat'));
            infolder = find(results);
            Slow_names = {};
            for l= 1:length(infolder)
                line = infolder(l);
                Slow_names{l} = [participant{p} '/@intra/' folders(line).name];
            end
            % Slowest
            results =  contains({folders.name},[modality_data{mode} '_' source_noise_tags{nopt} '_' wave{w}]) & ... 
                (endsWith({folders.name},'91.mat') | endsWith({folders.name},'92.mat') | ...
                endsWith({folders.name},'93.mat'));
            infolder = find(results);
            Slowest_names = {};
            for l= 1:length(infolder)
                line = infolder(l);
                Slowest_names{l} = [participant{p} '/@intra/' folders(line).name];
            end

            for ex = 1:length(Exp_cond) 
                eval(['sFiles =' Exp_cond{ex} '_names;'])          
                
                if ~isempty(sFiles)
                    
                    if delete_previous_file == 1
                        % check if there is already a source of THIS
                        % TYPE in the folder where we are gonna calculate it
                        folders_delete = dir([root_dir_bs '/data/' participant{p} '/@intra/']);
                        results_delete =  endsWith({folders_delete.name}, [modality_data{mode} '_' source_noise_tags{nopt} '_' wave{w} '_' Exp_cond{ex} '.mat']); 
                        infolder_delete = find(results_delete);
                        if ~isempty(infolder_delete) % file exists, therefore delete it
                           delete([root_dir_bs '/data/' participant{p} '/@intra/' folders_delete(infolder_delete).name]);
                        end
                        % source = contains KERNEL endswith _Source_Regul_MEG.mat                 
                    end
                    
                    disp(' ');
                    disp('-------------------------');  
                    disp(['Average sources across db and ISI (' modality_data{mode} '_' source_noise_tags{nopt} ') for ' participant{p} '_' wave{w}]);
                    disp(datetime)
                    disp(' ');

                    % WILL HAVE TO AVERAGE ALL OTHER SOURCE SOLUTIONS THAT MAY
                    % HAVE BEEN CREATED. FOR NOW, ONLY Overlapping spheres & Cortical surface

                    % Average using default function (because we are using vertices, not channels)
                    sFiles = bst_process('CallProcess', 'process_average', sFiles, [], ...
                        'avgtype',         1, ...  % Everything
                        'avg_func',        1, ...  % Arithmetic average:  mean(x)
                        'weighted',        1, ...
                        'scalenormalized', 0);

                    % USE OPTION TO NORMALIZE HERE???

                    % Process: Add tag
                    sFiles = bst_process('CallProcess', 'process_add_tag', sFiles, [], ...
                        'tag',           [modality_data{mode} '_' source_noise_tags{nopt} '_' wave{w} '_' Exp_cond{ex}], ...
                        'output',        2);  % Add to file name (1 to add a tag)

                    % Process: Set name
                    sFiles = bst_process('CallProcess', 'process_set_comment', sFiles, [], ...
                        'tag',           [modality_data{mode} '_' source_noise_tags{nopt} '_' wave{w} '_' Exp_cond{ex}], ...
                        'isindex',       1);
                end
            end
        end    
    end
    end
end

% Reset modality_data variable to its original self
modality_data = {'EEG','MEG','both_mod'};

clearvars('-except', initialVars{:});
disp 'DONE WITH AVERAGING SOURCES ACROSS dB and ISI (LLR)!!!'
disp(datetime)
toc

%% Collapse ALL SOURCES into a single one (regarless of condition)
tic
disp(' ');      
disp('-------------------------');  
disp('AVERAGING SOURCES ACROSS ALL CONDITIONS (LLR)'); 
disp(datetime)
disp('-------------------------');     
disp(' '); 

% Which ones will be averaged
if source_analysis == 1
    modality_data = {'EEG'};
elseif source_analysis == 2
    modality_data = {'MEG'};
elseif source_analysis == 3
    modality_data = {'both_mod'};
elseif source_analysis == 4
    modality_data = {'EEG','MEG','both_mod'};
end

for p = 1:length(participant)
    
    disp(' ');      
    disp('-------------------------');
    disp(['loading participant ' participant{p}]);
    disp(datetime)
    disp(' '); 
    
    % Reload only intra folder from subjects, as this is the one used
    [~, iStudies] = bst_get('StudyWithCondition', [participant{p} filesep '@intra']);
    db_reload_studies(iStudies);
    
    for w = 1:length(wave)
    for mode = 1:length(modality_data)
         for nopt = 1:length(source_noise_tags)
            folders = dir([root_dir_bs '/data/' participant{p} '/@intra/']);
            results =  contains({folders.name},[modality_data{mode} '_' source_noise_tags{nopt} '_' wave{w}]) & ... 
            (endsWith({folders.name},'11.mat') | endsWith({folders.name},'12.mat') | endsWith({folders.name},'13.mat') | ...
            endsWith({folders.name},'31.mat') | endsWith({folders.name},'32.mat') | endsWith({folders.name},'33.mat') | ...
            endsWith({folders.name},'51.mat') | endsWith({folders.name},'52.mat') | endsWith({folders.name},'53.mat') | ...
            endsWith({folders.name},'71.mat') | endsWith({folders.name},'72.mat') | endsWith({folders.name},'73.mat') | ...
            endsWith({folders.name},'91.mat') | endsWith({folders.name},'92.mat') | endsWith({folders.name},'93.mat'));
            infolder = find(results); 
            if isempty(infolder) % If there are files for MLR, for instance
                continue
            end
            file_names = {};
            for l= 1:length(infolder)
                line = infolder(l);
                file_names{l} = [participant{p} '/@intra/' folders(line).name];
            end
            file_names = file_names(~cellfun('isempty', file_names')); % to avoid empty cells (probably not necessary)
            sFiles = file_names;
            
            if isempty(sFiles)
            error('for whatever reason, sFiles is empty, my friend...')
            end
            
            if  ~isempty(sFiles)
                if delete_previous_file == 1
                    % check if there is already a source of THIS
                    % TYPE in the folder where we are gonna calculate it
                    folders_delete = dir([root_dir_bs '/data/' participant{p} '/@intra/']);
                    results_delete =  endsWith({folders_delete.name}, [modality_data{mode} '_' source_noise_tags{nopt} '_' wave{w} '_ALL.mat']); 
                    infolder_delete = find(results_delete);
                    if ~isempty(infolder_delete) % file exists, therefore delete it
                       delete([root_dir_bs '/data/' participant{p} '/@intra/' folders_delete(infolder_delete).name]);
                    end
                    % source = contains KERNEL endswith _Source_Regul_MEG.mat                 
                end

                disp(' ');
                disp('-------------------------');  
                disp(['Average sources across all conditions (' modality_data{mode} '_' source_noise_tags{nopt} ') for ' participant{p} '_' wave{w}]);
                disp(datetime)
                disp(' ');

                % WILL HAVE TO AVERAGE ALL OTHER SOURCE SOLUTIONS THAT MAY
                % HAVE BEEN CREATED. FOR NOW, ONLY Overlapping spheres & Cortical surface

                % Average using default function (because we are using vertices, not channels)
                sFiles = bst_process('CallProcess', 'process_average', sFiles, [], ...
                    'avgtype',         1, ...  % Everything
                    'avg_func',        1, ...  % Arithmetic average:  mean(x)
                    'weighted',        1, ...
                    'scalenormalized', 0);

                % USE OPTION TO NORMALIZE HERE???

                % Process: Add tag
                sFiles = bst_process('CallProcess', 'process_add_tag', sFiles, [], ...
                    'tag',           [modality_data{mode} '_' source_noise_tags{nopt} '_' wave{w} '_ALL'], ...
                    'output',        2);  % Add to file name (1 to add a tag)

                % Process: Set name
                sFiles = bst_process('CallProcess', 'process_set_comment', sFiles, [], ...
                    'tag',           [modality_data{mode} '_' source_noise_tags{nopt} '_' wave{w} '_ALL'], ...
                    'isindex',       1);
            end
        end    
    end
    end
end

% Reset modality_data variable to its original self
modality_data = {'EEG','MEG','both_mod'};

clearvars('-except', initialVars{:});
disp 'DONE WITH AVERAGING SOURCES ACROSS ALL CONDITIONS (LLR)!!!'
disp(datetime)
toc

%% Before GAVR source level, set to absolute and project to anatomy 

tic
disp(' ');      
disp('-------------------------');  
disp('SET SOURCES TO ABSOLUTE AND PROJECTING TO COMMON ANATOMY (LLR)'); 
disp(datetime)
disp('-------------------------');     
disp(' '); 

% Which ones will be averaged
if source_analysis == 1
    modality_data = {'EEG'};
elseif source_analysis == 2
    modality_data = {'MEG'};
elseif source_analysis == 3
    modality_data = {'both_mod'};
elseif source_analysis == 4
    modality_data = {'EEG','MEG','both_mod'};
end

% First, load all participants whose files are gonna be used
for p = 1:length(participant)  
    disp(' ');      
    disp('-------------------------');
    disp(['loading for GAVR participant ' participant{p}]);
    disp(datetime)
    disp(' '); 
    
    % Reload only intra folder from subjects, as this is the one used
    [~, iStudies] = bst_get('StudyWithCondition', [participant{p} filesep '@intra']);
    db_reload_studies(iStudies);
end

% Then, ensure appropiate cortex is loaded for cortical sources
load([root_dir '/anat/@default_subject/brainstormsubject.mat']);
if ~strcmp(Cortex,['@default_subject/' group_default_cortex])
    Cortex = ['@default_subject/' group_default_cortex];
    variableInfo = who('-file',[root_dir '/anat/@default_subject/brainstormsubject.mat']);
    save([root_dir '/anat/@default_subject/brainstormsubject.mat'],variableInfo{:});
end
% Finnally, reload group anatomy, as that folder is the destiny of the abs
db_reload_subjects(0); % As 0 is the default anatomy

% Abs values for individual condition (11, 12, etc) across subjects
for mode = 1:length(modality_data)
for w = 1:length(wave)
    for c = 1:length(condition)
         for nopt = 1:length(source_noise_tags)
             file_names = {};
             for p = 1:length(participant)
                folders = dir([root_dir_bs '/data/' participant{p} '/@intra/']);
                results =  contains({folders.name},[modality_data{mode} '_' source_noise_tags{nopt} '_' wave{w}]) & ...
                    endsWith({folders.name},[condition{c} '.mat']);
                infolder = find(results);
                if isempty(infolder) % In case, for instance, MLR files are not there yet
                    error(['No ' modality_data{mode} '_' source_noise_tags{nopt} '_' wave{w} '_' condition{c} ' file for participant' participant{p}]);
                elseif size(infolder,2)> 1 % in case of more than one coincidence (should not happen), pick first
                   infolder = infolder(1);
                end

                file_name = [participant{p} '/@intra/' folders(infolder).name];

                sFiles = file_name;

                if isempty(sFiles)
                error('for whatever reason, sFiles is empty, my friend...')
                end
                
                disp(' ');
                disp('-------------------------');  
                disp(['Obtain absolute source values (' modality_data{mode} '_' source_noise_tags{nopt} ') for ' participant{p} '_' wave{w} '_' condition{c}]);
                disp(datetime)
                disp(' ');

                % WILL HAVE TO DO THIS FOR ALL OTHER SOURCE SOLUTIONS THAT MAY
                % HAVE BEEN CREATED. FOR NOW, ONLY Overlapping spheres & Cortical surface
                
                % If stated, find and delete any previous ABS source in intra folder (overwrite it)
                if delete_previous_file == 1
                    % check if there is already a source ABS of this
                    % TYPE in the folder where we are gonna calculate it
                    % This is important because if recomputing sources the
                    % abs present in the folder may be the old one
                    folders_delete = dir([root_dir_bs '/data/' participant{p} '/@intra/']);
                    results_delete =  endsWith({folders_delete.name}, [modality_data{mode} '_' source_noise_tags{nopt} '_' wave{w} '_' condition{c} '_abs.mat']); 
                    infolder_delete = find(results_delete);
                    if ~isempty(infolder_delete) % file exists, therefore delete it
                       delete([root_dir_bs '/data/' participant{p} '/@intra/' folders_delete(infolder_delete).name]);
                    end
                end
                
                sFiles = bst_process('CallProcess', 'process_absolute', sFiles, [], ...
                    'overwrite', 0);
                
                % Leave default tag, which will be '_abs.mat', as in: 
                % '2250/@intra/results_average_200514_2211_Source_MEG_Regul_LLR_ALL_abs.mat'
                
                % If stated, find and delete any previous individual
                % sources in group average with projected sources
                if delete_previous_file == 1
                    % check if there is already an individual source in
                    % group folder (which would be the one with projected
                    % sources into common anatomy)
                    folders_delete = dir([root_dir_bs '/data/Group_analysis/@intra/']);
                    % results_average_200512_2239_Source_MEG_Regul_MLR_Quietest_abs_2005_2005_Source_MEG_Regul_MLR_Quietest
                    % results_average_200512_2003_Source_MEG_Regul_MLR_11_abs_2005_2005_Source_MEG_Regul_MLR_11
                    results_delete = contains({folders_delete.name},['_abs_' participant{p}])...
                    & endsWith({folders_delete.name}, [modality_data{mode} '_' source_noise_tags{nopt} '_' wave{w} '_' condition{c} '.mat']); 
                    infolder_delete = find(results_delete);
                    if ~isempty(infolder_delete) % file exists, therefore delete it
                       delete([root_dir_bs '/data/Group_analysis/@intra/' folders_delete(infolder_delete).name]);
                    end
                end

                % Process: Project on default anatomy: surface
                sFiles = bst_process('CallProcess', 'process_project_sources', sFiles, [], ...
                    'headmodeltype', 'surface');  % Cortex surface
                
                % Process: Add tag
                sFiles = bst_process('CallProcess', 'process_add_tag', sFiles, [], ...
                    'tag',           [participant{p} '_' modality_data{mode} '_' source_noise_tags{nopt} '_' wave{w} '_' condition{c}], ...
                    'output',        2);  % Add to file name (1 to add a tag)

                % Process: Set name
                sFiles = bst_process('CallProcess', 'process_set_comment', sFiles, [], ...
                    'tag',           [participant{p} '_' modality_data{mode} '_' source_noise_tags{nopt} '_' wave{w} '_' condition{c}], ...
                    'isindex',       1);
             end
        end    
    end
end
end

% Abs values for source categories (dB and ISI) across subjects
for mode = 1:length(modality_data)
for w = 1:length(wave)
    for ex = 1:length(Exp_cond)
         for nopt = 1:length(source_noise_tags)
             for p = 1:length(participant)
                folders = dir([root_dir_bs '/data/' participant{p} '/@intra/']);
                results =  contains({folders.name},[modality_data{mode} '_' source_noise_tags{nopt} '_' wave{w}]) & ...
                    endsWith({folders.name},[Exp_cond{ex} '.mat']);
                infolder = find(results);
                if isempty(infolder) % In case, for instance, MLR files are not there yet
                    error(['No ' modality_data{mode} '_' source_noise_tags{nopt} '_' wave{w} '_' Exp_cond{ex} ' file for participant' participant{p}]);
                elseif size(infolder,2)> 1 % in case of more than one coincidence (should not happen), pick first
                   infolder = infolder(1);
                end

                file_name = [participant{p} '/@intra/' folders(infolder).name];

                sFiles = file_name;

                if isempty(sFiles)
                error('for whatever reason, sFiles is empty, my friend...')
                end
                
                disp(' ');
                disp('-------------------------');  
                disp(['Obtain absolute source values (' modality_data{mode} '_' source_noise_tags{nopt} ') for ' participant{p} '_' wave{w} '_' Exp_cond{ex}]);
                disp(datetime)
                disp(' ');

                % WILL HAVE TO DO THIS FOR ALL OTHER SOURCE SOLUTIONS THAT MAY
                % HAVE BEEN CREATED. FOR NOW, ONLY Overlapping spheres & Cortical surface
                
                % If stated, find and delete any previous ABS source in intra folder (overwrite it)
                if delete_previous_file == 1
                    % check if there is already a source ABS of this
                    % TYPE in the folder where we are gonna calculate it
                    % This is important because if recomputing sources the
                    % abs present in the folder may be the old one
                    folders_delete = dir([root_dir_bs '/data/' participant{p} '/@intra/']);
                    results_delete =  endsWith({folders_delete.name}, [modality_data{mode} '_' source_noise_tags{nopt} '_' wave{w} '_' Exp_cond{ex} '_abs.mat']); 
                    infolder_delete = find(results_delete);
                    if ~isempty(infolder_delete) % file exists, therefore delete it
                       delete([root_dir_bs '/data/' participant{p} '/@intra/' folders_delete(infolder_delete).name]);
                    end
                end

                sFiles = bst_process('CallProcess', 'process_absolute', sFiles, [], ...
                    'overwrite', 0);
                
                % Leave default tag, which will be '_abs.mat', as in: 
                % '2250/@intra/results_average_200514_2211_Source_MEG_Regul_LLR_ALL_abs.mat'

                % If stated, find and delete any previous individual
                % sources in group average with projected sources
                if delete_previous_file == 1
                    % check if there is already an individual source in
                    % group folder (which would be the one with projected
                    % sources into common anatomy)
                    folders_delete = dir([root_dir_bs '/data/Group_analysis/@intra/']);
                    % results_average_200512_2239_Source_MEG_Regul_MLR_Quietest_abs_2005_2005_Source_MEG_Regul_MLR_Quietest
                    % results_average_200512_2003_Source_MEG_Regul_MLR_11_abs_2005_2005_Source_MEG_Regul_MLR_11
                    results_delete = contains({folders_delete.name},['_abs_' participant{p}])...
                    & endsWith({folders_delete.name}, [modality_data{mode} '_' source_noise_tags{nopt} '_' wave{w} '_' Exp_cond{ex} '.mat']); 
                    infolder_delete = find(results_delete);
                    if ~isempty(infolder_delete) % file exists, therefore delete it
                       delete([root_dir_bs '/data/Group_analysis/@intra/' folders_delete(infolder_delete).name]);
                    end
                end
                
                % Process: Project on default anatomy: surface
                sFiles = bst_process('CallProcess', 'process_project_sources', sFiles, [], ...
                    'headmodeltype', 'surface');  % Cortex surface
                
                % Process: Add tag
                sFiles = bst_process('CallProcess', 'process_add_tag', sFiles, [], ...
                    'tag',           [participant{p} '_' modality_data{mode} '_' source_noise_tags{nopt} '_' wave{w} '_' Exp_cond{ex}], ...
                    'output',        2);  % Add to file name (1 to add a tag)

                % Process: Set name
                sFiles = bst_process('CallProcess', 'process_set_comment', sFiles, [], ...
                    'tag',           [participant{p} '_' modality_data{mode} '_' source_noise_tags{nopt} '_' wave{w} '_' Exp_cond{ex}], ...
                    'isindex',       1);
             end
        end    
    end
end
end

% Abs values for averaged waveform ("ALL") sources across subjects
for mode = 1:length(modality_data)
for w = 1:length(wave)
     for nopt = 1:length(source_noise_tags)
         file_names = {};
         for p = 1:length(participant)
            folders = dir([root_dir_bs '/data/' participant{p} '/@intra/']);
            results =  contains({folders.name},[modality_data{mode} '_' source_noise_tags{nopt} '_' wave{w}]) & ...
                endsWith({folders.name},'ALL.mat');
            infolder = find(results);
            if isempty(infolder) % In case, for instance, MLR files are not there yet
                error(['No ' modality_data{mode} '_' source_noise_tags{nopt} '_' wave{w} '_ALL file for participant' participant{p}]);
            elseif size(infolder,2)> 1 % in case of more than one coincidence (should not happen), pick first
               infolder = infolder(1);
            end

            file_name = [participant{p} '/@intra/' folders(infolder).name];
            
            sFiles = file_name;
            
            if isempty(sFiles)
            error('for whatever reason, sFiles is empty, my friend...')
            end
            
                disp(' ');
                disp('-------------------------');  
                disp(['Obtain absolute source values (' modality_data{mode} '_' source_noise_tags{nopt} ') for ' participant{p} '_' wave{w} '_ALL']);
                disp(datetime)
                disp(' ');

                % WILL HAVE TO DO THIS FOR ALL OTHER SOURCE SOLUTIONS THAT MAY
                % HAVE BEEN CREATED. FOR NOW, ONLY Overlapping spheres & Cortical surface

                % If stated, find and delete any previous ABS source in intra folder (overwrite it)
                if delete_previous_file == 1
                    % check if there is already a source ABS of this
                    % TYPE in the folder where we are gonna calculate it
                    % This is important because if recomputing sources the
                    % abs present in the folder may be the old one
                    folders_delete = dir([root_dir_bs '/data/' participant{p} '/@intra/']);
                    results_delete =  endsWith({folders_delete.name}, [modality_data{mode} '_' source_noise_tags{nopt} '_' wave{w} '_ALL_abs.mat']); 
                    infolder_delete = find(results_delete);
                    if ~isempty(infolder_delete) % file exists, therefore delete it
                       delete([root_dir_bs '/data/' participant{p} '/@intra/' folders_delete(infolder_delete).name]);
                    end
                end
                
                sFiles = bst_process('CallProcess', 'process_absolute', sFiles, [], ...
                    'overwrite', 0);
                
                % Leave default tag, which will be '_abs.mat', as in: 
                % '2250/@intra/results_average_200514_2211_Source_MEG_Regul_LLR_ALL_abs.mat'  
                
                % If stated, find and delete any previous individual
                % sources in group average with projected sources
                if delete_previous_file == 1
                    % check if there is already an individual source in
                    % group folder (which would be the one with projected
                    % sources into common anatomy)
                    folders_delete = dir([root_dir_bs '/data/Group_analysis/@intra/']);
                    % results_average_200512_2239_Source_MEG_Regul_MLR_Quietest_abs_2005_2005_Source_MEG_Regul_MLR_Quietest
                    % results_average_200512_2003_Source_MEG_Regul_MLR_11_abs_2005_2005_Source_MEG_Regul_MLR_11
                    results_delete = contains({folders_delete.name},['_abs_' participant{p}])...
                    & endsWith({folders_delete.name}, [modality_data{mode} '_' source_noise_tags{nopt} '_' wave{w} '_ALL.mat']); 
                    infolder_delete = find(results_delete);
                    if ~isempty(infolder_delete) % file exists, therefore delete it
                       delete([root_dir_bs '/data/Group_analysis/@intra/' folders_delete(infolder_delete).name]);
                    end
                end
                
                % Process: Project on default anatomy: surface
                sFiles = bst_process('CallProcess', 'process_project_sources', sFiles, [], ...
                    'headmodeltype', 'surface');  % Cortex surface
                
                % Process: Add tag
                sFiles = bst_process('CallProcess', 'process_add_tag', sFiles, [], ...
                    'tag',           [participant{p} '_' modality_data{mode} '_' source_noise_tags{nopt} '_' wave{w} '_ALL'], ...
                    'output',        2);  % Add to file name (1 to add a tag)

                % Process: Set name
                sFiles = bst_process('CallProcess', 'process_set_comment', sFiles, [], ...
                    'tag',           [participant{p} '_' modality_data{mode} '_' source_noise_tags{nopt} '_' wave{w} '_ALL'], ...
                    'isindex',       1);
         end
    end    
end
end

% Reset modality_data variable to its original self
modality_data = {'EEG','MEG','both_mod'};

clearvars('-except', initialVars{:});
disp 'SOURCES TO ABSOLUTE AND PROJECTING TO COMMON ANATOMY (LLR)!!!'
disp(datetime)
toc

%% GAVR source level

tic
disp(' ');      
disp('-------------------------');  
disp('GAVR SOURCE DATA (LLR)'); 
disp(datetime)
disp('-------------------------');
disp(' ');

% To correct if we used a loop through sections
participant = participant_general_list;


if source_analysis == 1
    modality_data = {'EEG'};
elseif source_analysis == 2
    modality_data = {'MEG'};
elseif source_analysis == 3
    modality_data = {'both_mod'};
elseif source_analysis == 4
    modality_data = {'EEG','MEG','both_mod'};
end

% First, reload Group analysis folder
prot_subs = bst_get('ProtocolSubjects');
current_sub = find(strcmp({prot_subs.Subject.Name}, 'Group_analysis'));
db_reload_conditions(current_sub);

% Then, ensure appropiate cortex is loaded for cortical sources
load([root_dir '/anat/@default_subject/brainstormsubject.mat']);
if ~strcmp(Cortex,['@default_subject/' group_default_cortex])
    Cortex = ['@default_subject/' group_default_cortex];
    variableInfo = who('-file',[root_dir '/anat/@default_subject/brainstormsubject.mat']);
    save([root_dir '/anat/@default_subject/brainstormsubject.mat'],variableInfo{:});
end
% Finnally, reload group anatomy, as that folder is the destiny of the abs
db_reload_subjects(0); % As 0 is the default anatomy

% Average source data individual condition (11, 12, etc) across subjects
for mode = 1:length(modality_data)
for w = 1:length(wave)
    for c = 1:length(condition)
         for nopt = 1:length(source_noise_tags)
             file_names = {};
             for p = 1:length(participant)
                folders = dir([root_dir_bs '/data/Group_analysis/@intra/']);
                results =  endsWith({folders.name},[participant{p} '_' modality_data{mode} '_' source_noise_tags{nopt} '_' wave{w} '_' condition{c} '.mat']) ...
                    & ~contains({folders.name},'GAVR');
                infolder = find(results);
                if isempty(infolder) % In case, for instance, MLR files are not there yet
                    error(['No ' modality_data{mode} '_' source_noise_tags{nopt} '_' wave{w} '_' condition{c} ' file for participant' participant{p}]);
                elseif size(infolder,2)> 1 % in case of more than one coincidence (should not happen), pick first
                   error(['More than one ' modality_data{mode} '_' source_noise_tags{nopt} '_' wave{w} '_' condition{c} ' file for participant' participant{p}]);
                end

                file_names{p} = ['Group_analysis/@intra/' folders(infolder).name];
             end
                sFiles = file_names;

                if isempty(sFiles)
                error('for whatever reason, sFiles is empty, my friend...')
                end
                
                disp(' ');
                disp('-------------------------');  
                disp(['GAVR sources ' modality_data{mode} '_' source_noise_tags{nopt} ' for ' wave{w} '_' condition{c}]);
                disp(datetime)
                disp(' ');

                % WILL HAVE TO AVERAGE ALL OTHER SOURCE SOLUTIONS THAT MAY
                % HAVE BEEN CREATED. FOR NOW, ONLY Overlapping spheres & Cortical surface

                % If stated, find and delete any previous GAVR sources
                if delete_previous_file == 1
                    % check if there is already GAVR source in Group analysis folder
                    folders_delete = dir([root_dir_bs '/data/Group_analysis/@intra/']);
                    % results_average_200520_2029_GAVR_Source_MEG_Regul_MLR_Quietest
                    results_delete = contains({folders_delete.name},'_GAVR_')...
                    & endsWith({folders_delete.name}, [modality_data{mode} '_' source_noise_tags{nopt} '_' wave{w} '_' condition{c} '.mat']); 
                    infolder_delete = find(results_delete);
                    if ~isempty(infolder_delete) % file exists, therefore delete it
                       delete([root_dir_bs '/data/Group_analysis/@intra/' folders_delete(infolder_delete).name]);
                    end
                end
                
                % Average using default function (because we are using vertices, not channels)
                sFiles = bst_process('CallProcess', 'process_average', sFiles, [], ...
                    'avgtype',         1, ...  % Everything
                    'avg_func',        1, ...  % Arithmetic average:  mean(x)
                    'weighted',        0, ...
                    'scalenormalized', 0);

                % USE OPTION TO NORMALIZE HERE???

                % Process: Add tag
                sFiles = bst_process('CallProcess', 'process_add_tag', sFiles, [], ...
                    'tag',           ['GAVR_' modality_data{mode} '_' source_noise_tags{nopt} '_' wave{w} '_' condition{c}], ...
                    'output',        2);  % Add to file name (1 to add a tag)

                % Process: Set name
                sFiles = bst_process('CallProcess', 'process_set_comment', sFiles, [], ...
                    'tag',           ['GAVR_' modality_data{mode} '_' source_noise_tags{nopt} '_' wave{w} '_' condition{c}], ...
                    'isindex',       1);
        end    
    end
end
end

% Average source categories (dB and ISI) across subjects
for mode = 1:length(modality_data)
for w = 1:length(wave)
    for ex = 1:length(Exp_cond)
         for nopt = 1:length(source_noise_tags)
             file_names = {};
             for p = 1:length(participant)
                folders = dir([root_dir_bs '/data/Group_analysis/@intra/']);
                results =  endsWith({folders.name},[participant{p} '_' modality_data{mode} '_' source_noise_tags{nopt} '_' wave{w} '_' Exp_cond{ex} '.mat']) ...
                    & ~contains({folders.name},'GAVR');
                infolder = find(results);
                if isempty(infolder) % In case, for instance, MLR files are not there yet
                    error(['No ' modality_data{mode} '_' source_noise_tags{nopt} '_' wave{w} '_' Exp_cond{ex} ' file for participant' participant{p}]);
                elseif size(infolder,2)> 1 % in case of more than one coincidence (should not happen), pick first
                   error(['More than one ' modality_data{mode} '_' source_noise_tags{nopt} '_' wave{w} '_' Exp_cond{ex} ' file for participant' participant{p}]);
                end

                file_names{p} = ['Group_analysis/@intra/' folders(infolder).name];
             end
                sFiles = file_names;

                if isempty(sFiles)
                error('for whatever reason, sFiles is empty, my friend...')
                end
                
                disp(' ');
                disp('-------------------------');  
                disp(['GAVR sources ' modality_data{mode} '_' source_noise_tags{nopt} ' for ' wave{w} '_' Exp_cond{ex}]);
                disp(datetime)
                disp(' ');

                % WILL HAVE TO AVERAGE ALL OTHER SOURCE SOLUTIONS THAT MAY
                % HAVE BEEN CREATED. FOR NOW, ONLY Overlapping spheres & Cortical surface

                % If stated, find and delete any previous GAVR sources
                if delete_previous_file == 1
                    % check if there is already GAVR source in Group analysis folder
                    folders_delete = dir([root_dir_bs '/data/Group_analysis/@intra/']);
                    % results_average_200520_2029_GAVR_Source_MEG_Regul_MLR_Quietest
                    results_delete = contains({folders_delete.name},'_GAVR_')...
                    & endsWith({folders_delete.name}, [modality_data{mode} '_' source_noise_tags{nopt} '_' wave{w} '_' Exp_cond{ex} '.mat']); 
                    infolder_delete = find(results_delete);
                    if ~isempty(infolder_delete) % file exists, therefore delete it
                       delete([root_dir_bs '/data/Group_analysis/@intra/' folders_delete(infolder_delete).name]);
                    end
                end
                
                % Average using default function (because we are using vertices, not channels)
                sFiles = bst_process('CallProcess', 'process_average', sFiles, [], ...
                    'avgtype',         1, ...  % Everything
                    'avg_func',        1, ...  % Arithmetic average:  mean(x)
                    'weighted',        0, ...
                    'scalenormalized', 0);

                % USE OPTION TO NORMALIZE HERE???

                % Process: Add tag
                sFiles = bst_process('CallProcess', 'process_add_tag', sFiles, [], ...
                    'tag',           ['GAVR_' modality_data{mode} '_' source_noise_tags{nopt} '_' wave{w} '_' Exp_cond{ex}], ...
                    'output',        2);  % Add to file name (1 to add a tag)

                % Process: Set name
                sFiles = bst_process('CallProcess', 'process_set_comment', sFiles, [], ...
                    'tag',           ['GAVR_' modality_data{mode} '_' source_noise_tags{nopt} '_' wave{w} '_' Exp_cond{ex}], ...
                    'isindex',       1);
        end    
    end
end
end

% Average collapsed waveform ("ALL") sources across subjects
for mode = 1:length(modality_data)
for w = 1:length(wave)
     for nopt = 1:length(source_noise_tags)
         file_names = {};
         for p = 1:length(participant)
            folders = dir([root_dir_bs '/data/Group_analysis/@intra/']);
            results =  endsWith({folders.name},[participant{p} '_' modality_data{mode} '_' source_noise_tags{nopt} '_' wave{w} '_ALL.mat']) & ...
               ~contains({folders.name},'GAVR');
            infolder = find(results);
            if isempty(infolder) % In case, for instance, MLR files are not there yet
                error(['No ' modality_data{mode} '_' source_noise_tags{nopt} '_' wave{w} '_ALL file for participant' participant{p}]);
            elseif size(infolder,2)> 1 % in case of more than one coincidence (should not happen), pick first
               error(['More than one ' modality_data{mode} '_' source_noise_tags{nopt} '_' wave{w} '_ALL file for participant' participant{p}]);
            end

            file_names{p} = ['Group_analysis/@intra/' folders(infolder).name];
         end
            sFiles = file_names;

            if isempty(sFiles)
            error('for whatever reason, sFiles is empty, my friend...')
            end
            
            disp(' ');
            disp('-------------------------');  
            disp(['GAVR sources ' modality_data{mode} '_' source_noise_tags{nopt} ' for ' wave{w} '_ALL']);
            disp(datetime)
            disp(' ');

            % WILL HAVE TO AVERAGE ALL OTHER SOURCE SOLUTIONS THAT MAY
            % HAVE BEEN CREATED. FOR NOW, ONLY Overlapping spheres & Cortical surface

            % If stated, find and delete any previous GAVR sources
            if delete_previous_file == 1
                % check if there is already GAVR source in Group analysis folder
                folders_delete = dir([root_dir_bs '/data/Group_analysis/@intra/']);
                % results_average_200520_2029_GAVR_Source_MEG_Regul_MLR_Quietest
                results_delete = contains({folders_delete.name},'_GAVR_')...
                & endsWith({folders_delete.name}, [modality_data{mode} '_' source_noise_tags{nopt} '_' wave{w} '_ALL.mat']); 
                infolder_delete = find(results_delete);
                if ~isempty(infolder_delete) % file exists, therefore delete it
                   delete([root_dir_bs '/data/Group_analysis/@intra/' folders_delete(infolder_delete).name]);
                end
            end
            
            % Average using default function (because we are using vertices, not channels)
            sFiles = bst_process('CallProcess', 'process_average', sFiles, [], ...
                'avgtype',         1, ...  % Everything
                'avg_func',        1, ...  % Arithmetic average:  mean(x)
                'weighted',        0, ...
                'scalenormalized', 0);

            % USE OPTION TO NORMALIZE HERE???

            % Process: Add tag
            sFiles = bst_process('CallProcess', 'process_add_tag', sFiles, [], ...
                'tag',           ['GAVR_' modality_data{mode} '_' source_noise_tags{nopt} '_' wave{w} '_ALL'], ...
                'output',        2);  % Add to file name (1 to add a tag)

            % Process: Set name
            sFiles = bst_process('CallProcess', 'process_set_comment', sFiles, [], ...
                'tag',           ['GAVR_' modality_data{mode} '_' source_noise_tags{nopt} '_' wave{w} '_ALL'], ...
                'isindex',       1);
    end    
end
end

% Reset modality_data variable to its original self
modality_data = {'EEG','MEG','both_mod'};

clearvars('-except', initialVars{:});
disp 'DONE WITH GAVR SOURCES (LLR)!!!'
disp(datetime)
toc

%% Extract_source_waveforms_LLR
% Search for a separated script in Brainstorm Pipeline

%% Obtain SNR based on averaged ERP

tic
disp(' ');      
disp('-------------------------');  
disp('OBTAINING SNR BASED ON AVERAGED ERP (EEG ONLY)'); 
disp(datetime)
disp('-------------------------');
disp(' ');

% Always has to be EEG for now, as this is only relevant for FCz electrode
modality_data = {'EEG'};
SNR_average_table_MLR = {'Sub/cond' '11' '12' '13' '31' '32' '33' '51' '52' '53' '71' '72' '73' '91' '92' '93' 'AVERAGE'}; % Will be participants x conditions
SNR_average_table_LLR = {'Sub/cond' '11' '12' '13' '31' '32' '33' '51' '52' '53' '71' '72' '73' '91' '92' '93' 'AVERAGE'}; % Will be participants x conditions

for p = 1:length(participant) 
    for mode = 1:length(modality_data)
        for w = 1:length(wave)   
            % Create to store num (not chain) values to average
            eval(['Data_to_average_' wave{w} '= [];'])
            root_dir_bs = [root_dir '/brainstorm_db/' wave{w} 'seg']; 
            for c = 1:length(condition)
                
                folders = dir([root_dir_bs '/data/' participant{p} '/@intra']);
                results =  endsWith({folders.name},['Normal_' modality_data{mode} '_' wave{w} '_' condition{c} '.mat']);
                infolder = find(results); % So get all blocks you find for that condition and condition
                if isempty(infolder)
                    continue
                elseif length(infolder) > 1 % more than 1 coincidence
                    error(['more than one Normal_' modality_data{mode} '_' wave{w} '_' condition{c} ' found']);
                end            

                disp(' ');      
                disp('-------------------------');  
                disp(['Computing SNR (average) for ' participant{p} '_' modality_data{mode} '_' wave{w} '_' condition{c}]);
                disp(datetime)
                disp(' ');     
                
                current_file = [root_dir_bs '/data/' participant{p} '/@intra/' folders(infolder).name];
                
                % Set parameters to define what is baseline and what is signal in waves
                if wave{w} == 'MLR'
                    baseline = -50; % in ms
                    post = 100;
                    s_r = 1500;
                elseif wave{w} == 'LLR'
                    baseline = -150; 
                    post = 400;
                    s_r = 500;
                end
                time_samples=linspace(baseline,post,((((baseline*(-1)) + post)/1000)*s_r) +1);
                zero_time = find(time_samples == 0);
                
                % Load file
                load([root_dir_bs '/data/' participant{p} '/@intra/' folders(infolder).name]);
                
                % Compute SNR
                rms_signal = rms(F(335,zero_time+1:end)); % 335 is always FCz
                rms_noise = rms(F(335,1:zero_time-1));
                SNR = 10*log10(rms_signal/rms_noise);
                
                % Put subject name in variable
                eval(['SNR_average_table_' wave{w} '{p+1,c+1} = SNR;']);   
                
                % Store the value in num form to average later
                eval(['Data_to_average_' wave{w} '(c) = SNR;'])
            end
        end
    % Compute SNR average across conditions before moving to next subject
    eval(['SNR_average_table_' wave{w} '{p+1,17} = mean(Data_to_average_' wave{w} ');']);
    % Add row name for that subject
    eval(['SNR_average_table_' wave{w} '{p+1,1} = ''' participant{p} ''';']); 
    end
end

% Save the variables
try
save([root_dir '/Brainstorm_pipelines/SNR/SNR_average_table_MLR.mat'], 'SNR_average_table_MLR');
catch
    disp('MLR matrix not saved');
end
try
save([root_dir '/Brainstorm_pipelines/SNR/SNR_average_table_LLR.mat'], 'SNR_average_table_LLR');
catch
    disp('LLR matrix not saved');
end

% reset modality_data variable to its original
modality_data = {'EEG','MEG','both_mod'};

clearvars('-except', initialVars{:});
disp 'DONE WITH OBTAINING SNR BASED ON AVERAGED ERP (EEG ONLY)!!!'
disp(datetime)
toc

%% Obtain values for SNR single trials

tic
disp(' ');      
disp('-------------------------');  
disp('OBTAINING SIGNAL TO NOISE RATIO BASED ON INDIVIDUAL SWEEPS (EEG only)');  
disp(datetime)
disp('-------------------------');     
disp(' ');

modality_data = {'EEG'};

for p = 1:length(participant)
    for w = 1:length(wave)
        
    % Define root_dir_bs
    root_dir_bs = [root_dir '/brainstorm_db/' wave{w} 'seg']; 
    % Define variable to store single trial values for SNR
    eval(['SNR_matrix_' wave{w} '_' participant{p} ' = [];']);
        
    % Define sessions within this participant
    pos_par = find(strcmp(session_block_array(:,1), participant{p}(1:4)));
    session = session_block_array{pos_par,2}(1,:);
    
        for c = 1:length(condition)     
            
            % Re/define rows to start and end    
            rowStart = 1;
            if strcmp(wave{w},'MLR')
                rowEnd = (0.150*1500)+1; % size whole epoch times sr
            elseif strcmp(wave{w},'LLR')
                rowEnd = (0.550*500)+1; % size whole epoch times sr
            end
            
            for mode = 1:length(modality_data) % EEG and MEG
                for s = 1:length(session)
                    
                    % Define blocks within this session
                    pos_ses = find(strcmp(session_block_array{pos_par,2}(1,:), session{s}));
                    block = session_block_array{pos_par,2}{2,pos_ses}(1,:);
                    
                    for b = 1:length(block)  
                    % Search for files
                    folders = dir([root_dir_bs '/data/' participant{p} '/']);
                    results = contains({folders.name},[wave{w} '_' condition{c} session{s} '_' block{b} '_normal']);
                    infolder = find(results);
                    if isempty(infolder)
                        continue % should not be a problem unless all conditions are bad for a block, but in this case that block is not probably in the block list
                    end
                    dir_trial = dir([root_dir_bs '/data/' participant{p} '/' folders(infolder).name]);
                    trial_cov = contains({dir_trial.name},'_trial');
                    trial_count = find(trial_cov);
                    if isempty(trial_count)
                        continue % this is a rare exception, so just in case
                    end

                    % Load bad trials (EEG)
                    if ~exist([root_dir '/Events/BadTrials/' wave{w} '/' modality_data{mode} '/' participant{p} '/' folders(infolder).name '/BadTrials.mat'], 'file')
                        continue % if this folder does not have bad trials corresponding, better discard these sweeps
                    end
                    load([root_dir '/Events/BadTrials/' wave{w} '/' modality_data{mode} '/' participant{p} '/' folders(infolder).name '/BadTrials.mat']);

                    disp(' ');      
                    disp('-------------------------');  
                    disp(['Extracting ' wave{w} ' single trial matrix for SNR ' modality_data{mode} '_' participant{p} '_' session{s} '_' block{b} '_' condition{c}]);
                    disp(' ');  

                        for j= 1:length(trial_count) % for every trial
                            position = trial_count(j);
                            trial_name = [root_dir_bs '/data/' participant{p} '/' folders(infolder).name '/' dir_trial(position).name];
                            % If this trial is rejected, continue to next trial
                            if ~isempty(BadTrials)
                                if any(contains(BadTrials, [dir_trial(position).name]))
                                    disp('here')
                                    continue;
                                end
                            end      
                            eval(['load ' trial_name]) % trial loaded 
                            
                            % If for whatever reason size is not what's expected, abort!
                            if strcmp(wave{w},'MLR')
                                if size(F,2) ~= (0.150*1500)+1 % size whole epoch times sr
                                    continue
                                end
                            elseif strcmp(wave{w},'LLR')
                                if size(F,2) ~= (0.550*500)+1 % size whole epoch times sr
                                    continue
                                end
                            end
                                                                                    
                            % Store values in variable
                            eval(['SNR_matrix_' wave{w} '_' participant{p} '(rowStart:rowEnd,c) = F(335,:);']);
                            
                            % Increase size of rowStart and rowEnd as we sum trials
                            if strcmp(wave{w},'MLR')
                                rowStart = rowStart+(0.150*1500)+1;
                                rowEnd = rowEnd+(0.150*1500)+1;
                            elseif strcmp(wave{w},'LLR')
                                rowStart = rowStart+(0.550*500)+1;
                                rowEnd = rowEnd+(0.550*500)+1;
                            end
                        end
                    end
                end
            end
        end
    end
    % Save variables
    save([root_dir '/Brainstorm_pipelines/SNR/SNR_matrix_' wave{w} '_' participant{p} '.mat'],['SNR_matrix_' wave{w} '_' participant{p}]);
end

% Reset modality_data variable to its original self
modality_data = {'EEG','MEG','both_mod'};

clearvars('-except', initialVars{:});
disp 'DONE WITH OBTAINING SIGNAL TO NOISE RATIO BASED ON INDIVIDUAL SWEEPS (EEG only)!!!'
disp(datetime)
toc

%% Calculate bootstrap SNR

% Based on method by Parks, Gannon, Long and Young (2016)
% https://www.frontiersin.org/articles/10.3389/fnhum.2016.00050/full
% Available source code at: https://figshare.com/s/f6da4150953b0f9cc3bd.

% Add function to path
addpath([root_dir '/Brainstorm_pipelines/SNR']);

tic
disp(' ');      
disp('-------------------------');  
disp('OBTAINING BOOTSTRAP SNR BASED ON SINGLE TRIALS'); 
disp(datetime)
disp('-------------------------');
disp(' ');

% Always has to be EEG for now, as this is only relevant for FCz electrode
modality_data = {'EEG'};
SNR_bootstrap_table_MLR = {'Sub/cond' '11_lb' '11_ub' '11_snr' '12_lb' '12_ub' '12_snr' '13_lb' '13_ub' '13_snr' '31_lb' '31_ub' '31_snr' '32_lb' '32_ub' '32_snr' '33_lb' '33_ub' '33_snr' '51_lb' '51_ub' '51_snr' '52_lb' '52_ub' '52_snr' '53_lb' '53_ub' '53_snr' '71_lb' '71_ub' '71_snr' '72_lb' '72_ub' '72_snr' '73_lb' '73_ub' '73_snr' '91_lb' '91_ub' '91_snr' '92_lb' '92_ub' '92_snr' '93_lb' '93_ub' '93_snr' 'AVERAGE_SNR'};
SNR_bootstrap_table_LLR = {'Sub/cond' '11_lb' '11_ub' '11_snr' '12_lb' '12_ub' '12_snr' '13_lb' '13_ub' '13_snr' '31_lb' '31_ub' '31_snr' '32_lb' '32_ub' '32_snr' '33_lb' '33_ub' '33_snr' '51_lb' '51_ub' '51_snr' '52_lb' '52_ub' '52_snr' '53_lb' '53_ub' '53_snr' '71_lb' '71_ub' '71_snr' '72_lb' '72_ub' '72_snr' '73_lb' '73_ub' '73_snr' '91_lb' '91_ub' '91_snr' '92_lb' '92_ub' '92_snr' '93_lb' '93_ub' '93_snr' 'AVERAGE_SNR'};

for p = 1:length(participant) 
    for mode = 1:length(modality_data)
        for w = 1:length(wave)   
            if ~exist([root_dir '/Brainstorm_pipelines/SNR/SNR_matrix_' wave{w} '_' participant{p} '.mat'],'file')
                continue;
            end
            load([root_dir '/Brainstorm_pipelines/SNR/SNR_matrix_' wave{w} '_' participant{p} '.mat']);
            
            % Create to store num (not chain) values to average
            eval(['Data_to_average_' wave{w} '= [];'])
            % Reset position to fill confidence intervals and mean snr
            pos = 0;
            for c = 1:length(condition)
                
                disp(' ');      
                disp('-------------------------');  
                disp(['Computing SNR (bootstrap) for ' participant{p} '_' modality_data{mode} '_' wave{w} '_' condition{c}]);
                disp(datetime)
                disp(' ');     
                
                % Set parameters to define what is baseline and what is signal in waves
                if wave{w} == 'MLR'
                    baseline = -50; % in ms
                    post = 100;
                    s_r = 1500;
                elseif wave{w} == 'LLR'
                    baseline = -150; 
                    post = 400;
                    s_r = 500;
                end
                time_samples=linspace(baseline,post,((((baseline*(-1)) + post)/1000)*s_r) +1);
                zero_time = find(time_samples == 0);

                % Prepare input for function
                eval(['segment_dat = SNR_matrix_' wave{w} '_' participant{p} '(:,c);'])
                segment_points = [-zero_time length(time_samples)-zero_time];
                
                [snr_lb, snr_ub, snr_mean] = erp_snr_ci(segment_dat, segment_points, .90, 9999, [], 100, []);
                
                % Put subject name in variable
                % 1, 2, 3, 4, 5, 6, 
                eval(['SNR_bootstrap_table_' wave{w} '{p+1,c+1+pos} = snr_lb;']);
                eval(['SNR_bootstrap_table_' wave{w} '{p+1,c+1+1+pos} = snr_ub;']);  
                eval(['SNR_bootstrap_table_' wave{w} '{p+1,c+1+2+pos} = snr_mean;']);  
                
                % To properly fill table
                pos = pos+2;
                % Store the value in num form to average later
                eval(['Data_to_average_' wave{w} '(c) = snr_mean;'])
            end
        end
    % Compute SNR average across conditions before moving to next subject
    eval(['SNR_bootstrap_table_' wave{w} '{p+1,47} = mean(Data_to_average_' wave{w} ');']);
    % Add row name for that subject
    eval(['SNR_bootstrap_table_' wave{w} '{p+1,1} = ''' participant{p} ''';']); 
    % Clear single trial matrix for this subject before moving to next
    clear(['SNR_matrix_' wave{w} '_' participant{p}])
    end
end

% Save the variables
try
save([root_dir '/Brainstorm_pipelines/SNR/SNR_bootstrap_table_MLR.mat'], 'SNR_bootstrap_table_MLR');
catch
    disp('MLR matrix not saved');
end
try
save([root_dir '/Brainstorm_pipelines/SNR/SNR_bootstrap_table_LLR.mat'], 'SNR_bootstrap_table_LLR');
catch
    disp('LLR matrix not saved');
end

% reset modality_data variable to its original
modality_data = {'EEG','MEG','both_mod'};

clearvars('-except', initialVars{:});
disp 'DONE WITH OBTAINING BOOTSTRAP SNR BASED ON SINGLE TRIALS (EEG ONLY)!!!'
disp(datetime)
toc

%% Extract event times from all conditions (for SINGLE TRIAL AMPLITUDES, LLR)

tic
disp(' ');      
disp('-------------------------');  
disp('OBTAINING EVENT TIMES & ISI (FOR SINGLE TRIAL AMPLITUDES)');  
disp(datetime)
disp('-------------------------');     
disp(' ');

% First obtain original latency logs from AMICA files, which coincide with
% Brainstorm times
coincidence = 0; 
for p = 1:length(participant)
    
    % Define sessions within this participant
    pos_par = find(strcmp(session_block_array(:,1), participant{p}(1:4)));
    session = session_block_array{pos_par,2}(1,:);
    
    for s = 1:length(session)
        
        % Define runs within this block (but this time based on all present in folder)        
        run_files =dir([root_dir '/analysis/' participant{p}(1:4) session{s} '/'...
        participant{p} session{s} '_Click*_tsss_AMICA.vmrk']);

        runs_in_block = {};
        for i = 1:length(run_files)
            end_run_name = regexp(run_files(i).name, '_mc_tsss', 'once');
            runs_in_block{1,i} = run_files(i).name(7:end_run_name);
        end

        for r = 1:length(runs_in_block)

            folders = dir([root_dir '/analysis/' participant{p}(1:4) session{s} '/']);
            results =  contains({folders.name},runs_in_block{r}) & endsWith({folders.name},'AMICA.vmrk');
            infolder = find(results);
                if isempty(infolder)
                    continue
                end
                if size(infolder,2)> 1 % in case of more than one coincidence
                    infolder = infolder(1);
                end
            filename = folders(infolder).name;

            fullfilename = [root_dir '/analysis/' participant{p}(1:4) session{s} '/' filename];
            if ~exist(fullfilename, 'file')
                disp (['filename ' fullfilename ' does not exist'])
                continue
            end
            fileID = fopen(fullfilename,'r');
            textscan(fileID, '%[^\n\r]', startRow-1, 'WhiteSpace', '', 'ReturnOnError', false, 'EndOfLine', '\r\n');
            dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'ReturnOnError', false);
            fclose(fileID);
            current_run = runs_in_block{r}(10:end-1); % extract the number from the "Click_run9_", even with two numbers (10)
            eval(['eve_' participant{p}(1:4) session{s} '_r' current_run '= [dataArray{1:end-1}];'])
            eval(['events = eve_' participant{p}(1:4) session{s} '_r' current_run ';'])
            latencies_11 = []; count_11 = 1;
            latencies_12 = []; count_12 = 1;
            latencies_13 = []; count_13 = 1;
            latencies_31 = []; count_31 = 1;
            latencies_32 = []; count_32 = 1;
            latencies_33 = []; count_33 = 1;
            latencies_51 = []; count_51 = 1;
            latencies_52 = []; count_52 = 1;
            latencies_53 = []; count_53 = 1;
            latencies_71 = []; count_71 = 1;
            latencies_72 = []; count_72 = 1;
            latencies_73 = []; count_73 = 1;
            latencies_91 = []; count_91 = 1;
            latencies_92 = []; count_92 = 1;
            latencies_93 = []; count_93 = 1;
            for c = 1:length(condition)
                for e = 1:length(events)
                    if events(e,1)== condition{c}
                        % Edit: This count backwards does not make much
                        % sense, since it may consider the previous
                        % category or calculate the ISI using an stimulus
                        % ocurring before a gap, which is incorrect. For
                        % the purpose of subtraction, this was done
                        % correctly since any stimulus after a boundary was
                        % discarded for subtraction (more cautious), but
                        % for retrieving single trial data this should not
                        % be used. AND THAT IS THE CASE, we don't use these
                        % logs, but the next ones in which no boundaries
                        % are present. We only use these first ones to
                        % correct the little inconsistencies in time with
                        % the next ones (the non-AMICA vmrk extracted),
                        % which are the Logs to be actually used when
                        % retrieving single trial data
                        if e <= 5 % if less than 5 positions back are available, use the iteration size as a measure
                            for res = 1:(e-1) % check 'res' positions back until finding a coincidence with a meaningful trigger
                                if coincidence == 1; continue; end % if already found one n-back with the number, stop
                                if (events(e-res,1) == '11') || (events(e-res,1) == '12') || (events(e-res,1) == '13') ...
                                        || (events(e-res,1) == '31') || (events(e-res,1) == '32') || (events(e-res,1) == '33') ...
                                        || (events(e-res,1) == '51') || (events(e-res,1) == '52') || (events(e-res,1) == '53') ...
                                        || (events(e-res,1) == '71') || (events(e-res,1) == '72') || (events(e-res,1) == '73') ...
                                        || (events(e-res,1) == '91') || (events(e-res,1) == '92') || (events(e-res,1) == '93')
                                    coincidence = 1; % if reaching here, set this not to repeat
                                    % Absolute time (first column, in seconds)
                                    eval(['latencies_' condition{c} '(count_' condition{c} ',1) = str2double(events(e,2))/1500;']) % since this is vmrk archive, no resampling values are necessary
                                    % Latency from previous trial (second column, in seconds)
                                    eval(['latencies_' condition{c} '(count_' condition{c} ',2) = (str2double(events(e,2)) - str2double(events(e-res,2)))/1500;']) % since this is vmrk archive, no resampling values are necessary
                                    % Category of previous trial (third column)
                                    eval(['latencies_' condition{c} '(count_' condition{c} ',3) = str2double(events(e-res,1));'])
                                    eval(['count_' condition{c} ' = count_' condition{c} ' + 1;'])
                                end
                            end
                            coincidence = 0; % reset coincidence after finishing the loop
                        else
                            for res = 1:5 % number of positions back in which you are gonna look (5 for now)
                                if coincidence == 1; continue; end % if already found one n-back with the number, stop
                                if (events(e-res,1) == '11') || (events(e-res,1) == '12') || (events(e-res,1) == '13') ...
                                        || (events(e-res,1) == '31') || (events(e-res,1) == '32') || (events(e-res,1) == '33') ...
                                        || (events(e-res,1) == '51') || (events(e-res,1) == '52') || (events(e-res,1) == '53') ...
                                        || (events(e-res,1) == '71') || (events(e-res,1) == '72') || (events(e-res,1) == '73') ...
                                        || (events(e-res,1) == '91') || (events(e-res,1) == '92') || (events(e-res,1) == '93')            
                                    coincidence = 1; % if reaching here, set this not to repeat
                                    % Absolute time (first column, in seconds)
                                    eval(['latencies_' condition{c} '(count_' condition{c} ',1) = str2double(events(e,2))/1500;']) % since this is vmrk archive, no resampling values are necessary
                                    % Latency from previous trial (second column, in seconds)
                                    eval(['latencies_' condition{c} '(count_' condition{c} ',2) = (str2double(events(e,2)) - str2double(events(e-res,2)))/1500;']) % since this is vmrk archive, no resampling values are necessary
                                    % Category of previous trial (third column)
                                    eval(['latencies_' condition{c} '(count_' condition{c} ',3) = str2double(events(e-res,1));'])
                                    eval(['count_' condition{c} ' = count_' condition{c} ' + 1;'])
                                end
                            end
                            coincidence = 0; % reset coincidence after finishing the loop
                        end
                    end
                end
            end
            eval(['Ori_Log_11_' participant{p} session{s} '_r' current_run '= latencies_11;'])
            eval(['Ori_Log_12_' participant{p} session{s} '_r' current_run '= latencies_12;'])
            eval(['Ori_Log_13_' participant{p} session{s} '_r' current_run '= latencies_13;'])
            eval(['Ori_Log_31_' participant{p} session{s} '_r' current_run '= latencies_31;'])
            eval(['Ori_Log_32_' participant{p} session{s} '_r' current_run '= latencies_32;'])
            eval(['Ori_Log_33_' participant{p} session{s} '_r' current_run '= latencies_33;'])
            eval(['Ori_Log_51_' participant{p} session{s} '_r' current_run '= latencies_51;'])
            eval(['Ori_Log_52_' participant{p} session{s} '_r' current_run '= latencies_52;'])
            eval(['Ori_Log_53_' participant{p} session{s} '_r' current_run '= latencies_53;'])
            eval(['Ori_Log_71_' participant{p} session{s} '_r' current_run '= latencies_71;'])
            eval(['Ori_Log_72_' participant{p} session{s} '_r' current_run '= latencies_72;'])
            eval(['Ori_Log_73_' participant{p} session{s} '_r' current_run '= latencies_73;'])
            eval(['Ori_Log_91_' participant{p} session{s} '_r' current_run '= latencies_91;'])
            eval(['Ori_Log_92_' participant{p} session{s} '_r' current_run '= latencies_92;'])
            eval(['Ori_Log_93_' participant{p} session{s} '_r' current_run '= latencies_93;'])
            clearvars filename fileID dataArray ans;
        end
    end
    save([root_dir '/Events/Single_trials/Latencies_' participant{p} '_original.mat'],'-regexp','^Ori_Log_')
    
    % clear logs form this participant before loading next ones
    clearvars -regexp ^Log_
    clearvars -regexp ^Ori_Log_ 
    
end % Probably useless, since correction is probably not needed (>10>2)
clearvars('-except', initialVars{:});

% Next, obtain latency times from vmrk files of raw fif files (which
% contain removed sweeps at preprocessing)
coincidence = 0; 
for p = 1:length(participant)
    
    % Define sessions within this participant
    pos_par = find(strcmp(session_block_array(:,1), participant{p}(1:4)));
    session = session_block_array{pos_par,2}(1,:);
    
    for s = 1:length(session)
        
        % Define runs within this block (but this time based on all present in folder)        
        run_files =dir([root_dir '/Events/complete_event_log/' participant{p}(1:4) session{s} '/'...
        participant{p} session{s} '_Click*_tsss.vmrk']);

        runs_in_block = {};
        for i = 1:length(run_files)
            end_run_name = regexp(run_files(i).name, '_mc_tsss', 'once');
            runs_in_block{1,i} = run_files(i).name(7:end_run_name);
        end

        for r = 1:length(runs_in_block)

            folders = dir([root_dir '/Events/complete_event_log/' participant{p}(1:4) session{s} '/']);
            results =  contains({folders.name},runs_in_block{r}) & endsWith({folders.name},'vmrk');
            infolder = find(results);
                if isempty(infolder)
                    continue
                end
                if size(infolder,2)> 1 % in case of more than one coincidence, vmrk is going to be the same for all, so pick first
                    infolder = infolder(1);
                end
            filename = folders(infolder).name;

            fullfilename = [root_dir '/Events/complete_event_log/' participant{p}(1:4) session{s} '/' filename];
            if ~exist(fullfilename, 'file')
                disp (['filename ' fullfilename ' does not exist'])
                continue
            end
            fileID = fopen(fullfilename,'r');
            textscan(fileID, '%[^\n\r]', startRow-1, 'WhiteSpace', '', 'ReturnOnError', false, 'EndOfLine', '\r\n');
            dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'ReturnOnError', false);
            fclose(fileID);
            current_run = runs_in_block{r}(10:end-1); % extract the number from the "Click_run9_", even with two numbers (10)
            eval(['eve_' participant{p}(1:4) session{s} '_r' current_run '= [dataArray{1:end-1}];'])
            eval(['events = eve_' participant{p}(1:4) session{s} '_r' current_run ';'])
            latencies_11 = []; count_11 = 1;
            latencies_12 = []; count_12 = 1;
            latencies_13 = []; count_13 = 1;
            latencies_31 = []; count_31 = 1;
            latencies_32 = []; count_32 = 1;
            latencies_33 = []; count_33 = 1;
            latencies_51 = []; count_51 = 1;
            latencies_52 = []; count_52 = 1;
            latencies_53 = []; count_53 = 1;
            latencies_71 = []; count_71 = 1;
            latencies_72 = []; count_72 = 1;
            latencies_73 = []; count_73 = 1;
            latencies_91 = []; count_91 = 1;
            latencies_92 = []; count_92 = 1;
            latencies_93 = []; count_93 = 1;
            for c = 1:length(condition)
                for e = 1:length(events)
                    if events(e,1)== condition{c}    
                        if e <= 5 % if less than 5 positions back are available, use the iteration size as a measure
                            for res = 1:(e-1) % check 'res' positions back until finding a coincidence with a meaningful trigger
                                if coincidence == 1; continue; end % if already found one n-back with the number, stop
                                if (events(e-res,1) == '11') || (events(e-res,1) == '12') || (events(e-res,1) == '13') ...
                                        || (events(e-res,1) == '31') || (events(e-res,1) == '32') || (events(e-res,1) == '33') ...
                                        || (events(e-res,1) == '51') || (events(e-res,1) == '52') || (events(e-res,1) == '53') ...
                                        || (events(e-res,1) == '71') || (events(e-res,1) == '72') || (events(e-res,1) == '73') ...
                                        || (events(e-res,1) == '91') || (events(e-res,1) == '92') || (events(e-res,1) == '93')
                                    coincidence = 1; % if reaching here, set this not to repeat
                                    % Absolute time (first column, in seconds)
                                    eval(['latencies_' condition{c} '(count_' condition{c} ',1) = str2double(events(e,2))/1500;']) % since this is vmrk archive, no resampling values are necessary
                                    % Latency from previous trial (second column, in seconds)
                                    eval(['latencies_' condition{c} '(count_' condition{c} ',2) = (str2double(events(e,2)) - str2double(events(e-res,2)))/1500;']) % since this is vmrk archive, no resampling values are necessary
                                    % Category of previous trial (third column)
                                    eval(['latencies_' condition{c} '(count_' condition{c} ',3) = str2double(events(e-res,1));'])
                                    eval(['count_' condition{c} ' = count_' condition{c} ' + 1;'])
                                end
                            end
                            coincidence = 0; % reset coincidence after finishing the loop
                        else
                            for res = 1:5 % number of positions back in which you are gonna look (5 for now)
                                if coincidence == 1; continue; end % if already found one n-back with the number, stop
                                if (events(e-res,1) == '11') || (events(e-res,1) == '12') || (events(e-res,1) == '13') ...
                                        || (events(e-res,1) == '31') || (events(e-res,1) == '32') || (events(e-res,1) == '33') ...
                                        || (events(e-res,1) == '51') || (events(e-res,1) == '52') || (events(e-res,1) == '53') ...
                                        || (events(e-res,1) == '71') || (events(e-res,1) == '72') || (events(e-res,1) == '73') ...
                                        || (events(e-res,1) == '91') || (events(e-res,1) == '92') || (events(e-res,1) == '93')            
                                    coincidence = 1; % if reaching here, set this not to repeat
                                    % Absolute time (first column, in seconds)
                                    eval(['latencies_' condition{c} '(count_' condition{c} ',1) = str2double(events(e,2))/1500;']) % since this is vmrk archive, no resampling values are necessary
                                    % Latency from previous trial (second column, in seconds)
                                    eval(['latencies_' condition{c} '(count_' condition{c} ',2) = (str2double(events(e,2)) - str2double(events(e-res,2)))/1500;']) % since this is vmrk archive, no resampling values are necessary
                                    % Category of previous trial (third column)
                                    eval(['latencies_' condition{c} '(count_' condition{c} ',3) = str2double(events(e-res,1));'])
                                    eval(['count_' condition{c} ' = count_' condition{c} ' + 1;'])
                                end
                            end
                            coincidence = 0; % reset coincidence after finishing the loop
                        end
                    end
                end
            end
            eval(['Log_11_' participant{p} session{s} '_r' current_run '= latencies_11;'])
            eval(['Log_12_' participant{p} session{s} '_r' current_run '= latencies_12;'])
            eval(['Log_13_' participant{p} session{s} '_r' current_run '= latencies_13;'])
            eval(['Log_31_' participant{p} session{s} '_r' current_run '= latencies_31;'])
            eval(['Log_32_' participant{p} session{s} '_r' current_run '= latencies_32;'])
            eval(['Log_33_' participant{p} session{s} '_r' current_run '= latencies_33;'])
            eval(['Log_51_' participant{p} session{s} '_r' current_run '= latencies_51;'])
            eval(['Log_52_' participant{p} session{s} '_r' current_run '= latencies_52;'])
            eval(['Log_53_' participant{p} session{s} '_r' current_run '= latencies_53;'])
            eval(['Log_71_' participant{p} session{s} '_r' current_run '= latencies_71;'])
            eval(['Log_72_' participant{p} session{s} '_r' current_run '= latencies_72;'])
            eval(['Log_73_' participant{p} session{s} '_r' current_run '= latencies_73;'])
            eval(['Log_91_' participant{p} session{s} '_r' current_run '= latencies_91;'])
            eval(['Log_92_' participant{p} session{s} '_r' current_run '= latencies_92;'])
            eval(['Log_93_' participant{p} session{s} '_r' current_run '= latencies_93;'])
            clearvars filename fileID dataArray ans;
        end
    end
    save([root_dir '/Events/Single_trials/Latencies_' participant{p} '.mat'],'-regexp','^Log_') 
    
    % clear logs form this participant before loading next ones
    clearvars -regexp ^Log_
    clearvars -regexp ^Ori_Log_ 
    
end % Needed


clearvars('-except', initialVars{:});
disp 'DONE WITH OBTAINING EVENT TIMES & ISI (FOR SINGLE TRIAL AMPLITUDES)!!!'
disp(datetime)
toc


for p = 1:length(participant)
    
    % Load 
    load([root_dir '/Events/Single_trials/Latencies_' participant{p} '.mat'])
    load([root_dir '/Events/Single_trials/Latencies_' participant{p} '_original.mat'])
    
    % Define sessions within this participant
    pos_par = find(strcmp(session_block_array(:,1), participant{p}(1:4)));
    session = session_block_array{pos_par,2}(1,:);
    
    for s = 1:length(session)
        
        % Define runs within this block (but this time based on all present in folder)        
        run_files =dir([root_dir '/Events/complete_event_log/' participant{p}(1:4) session{s} '/'...
        participant{p} session{s} '_Click*_tsss.vmrk']);

        % This now is just to define which runs had Log files
        runs_in_block = {};
        for i = 1:length(run_files)
            end_run_name = regexp(run_files(i).name, '_mc_tsss', 'once');
            runs_in_block{1,i} = run_files(i).name(7:end_run_name); %#ok<SAGROW>
        end

        for r = 1:length(runs_in_block)
            
            current_run = runs_in_block{r}(10:end-1); % extract the number from the "Click_run9_", even with two numbers (10)
            
            for c = 1:length(condition)
                
                % Find the closest correspondant (there will be a tiny difference that
                % we aim to correct, but for that I need to find the equivalent value,
                % which will be extremely similar but with a tiny difference)
                % If I set the difference too strictly, two times that are actually the
                % same will be considered different and I won't be able to correct in
                % that log. On the other hand, if I set the difference too loose, then
                % it may match two completely different triggers (highly unlikely, as
                % the difference almost never reaches the decimals).

                % Check first values in ori_log (as any value in Ori_log will be in the normal log but not the other way around)    
                % Check first if Ori Log exists or is empty
                if ~exist(['Ori_Log_' condition{c} '_' participant{p} session{s} '_r' current_run],'var')
                    % It could be that a run exist in original fif files
                    % but not in AMICA files in analysis folder, in which case
                    % original logs won't match the current search of runs we
                    % are running here, based on the original fifs. So if a
                    % orig log is not present for comparison, leave. The sweeps
                    % there won't have been analyzed in Brainstorm them, so we
                    % will not need the number to be corrected by that tiny
                    % diffference in time to match it with brainstorm sweeps
                    continue
                end
                
                eval(['fix_bug = isempty(Ori_Log_' condition{c} '_' participant{p} session{s} '_r' current_run ');'])
                if fix_bug == 1; continue; end
                
                % For whatever reason this was not working properly
                % if isempty(['Ori_Log_' condition{c} '_' participant{p} session{s} '_r' current_run]); continue; end

                eval(['trial_time_example = Ori_Log_' condition{c} '_' participant{p} session{s} '_r' current_run '(1,1);'])

                % Find that value in complete Log
                % It's impossible that if two numbers coincide including ONE 
                % decimal (ex. 172.91 and 172.92) these are two different
                % sweeps or triggers, because at least there has to be 0.250 ms
                % of distance. In case of two repeated sweeps
                eval(['ident_pos_example = find(abs(trial_time_example  - Log_' condition{c} '_' participant{p} session{s} '_r' current_run '(:,1))<10^-1);']) 
                
                % Just in case a trigger is doubled, take first coincidence
                if length(ident_pos_example)>1
                    ident_pos_example = ident_pos_example(1);
                end
                
                % Extremely unlikely, but just in case
                if isempty(ident_pos_example)
                    % Extremely unlikely that Orig Log has a value and no
                    % correspondant in normal Log, that would mean that all values
                    % in normal Log are to be included a part, and therefore no
                    % correction is necessary
                    continue
                end

                % Now find what's the difference between the original and the
                % new obtained Log (should always be in the order of few ms)
                eval(['correction_factor = Ori_Log_' condition{c} '_' participant{p} session{s} '_r' current_run '(1,1) - Log_' condition{c} '_' participant{p} session{s} '_r' current_run '(ident_pos_example,1);'])

                % Once the correspondant has been identified, that's all I need to find
                % the difference and correct it in the rest of values from the Log
                % (including the non-presented sweeps that are not present in the
                % original Log). The difference in time with the previous sweep will
                % remain the same, as the difference within all the values of the Log
                % is consistent)

                % Now correct every value found in normal log with that correction factor
                eval(['size_log = size(Log_' condition{c} '_' participant{p} session{s} '_r' current_run ',1);'])
                for sl = 1:size_log
                    if correction_factor>=0 % so if it's positive or there are no differences
                        eval(['Log_' condition{c} '_' participant{p} session{s} '_r' current_run '(sl,1) = Log_' condition{c} '_' participant{p} session{s} '_r' current_run '(sl,1) + correction_factor;'])
                    else % rare, but if symbol is negative
                        eval(['Log_' condition{c} '_' participant{p} session{s} '_r' current_run '(sl,1) = Log_' condition{c} '_' participant{p} session{s} '_r' current_run '(sl,1)  - (correction_factor*-1);'])
                    end        
                end
            end  
        end
    end
    % now save Log files (not original ones) back to where they were, only
    % now their times are corrected
    save([root_dir '/Events/Single_trials/Latencies_' participant{p} '.mat'],'-regexp','^Log_') 
    
    % clear logs form this participant before loading next ones
    clearvars -regexp ^Log_
    clearvars -regexp ^Ori_Log_ 
end

%% Extract complete log for correction (for SINGLE TRIAL AMPLITUDES, LLR and MLR)

tic
disp(' ');      
disp('-------------------------');  
disp('OBTAINING MATLAB COMPLETE LOG (TO FILL GAPS IN FINAL SINGLE TRIAL MATRIX)');  
disp(datetime)
disp('-------------------------');     
disp(' ');

% Next, obtain latency times from vmrk files of raw fif files (which
% contain removed sweeps at preprocessing)
coincidence = 0; 
for p = 1:length(participant)
    
    % Define sessions within this participant
    pos_par = find(strcmp(session_block_array(:,1), participant{p}(1:4)));
    session = session_block_array{pos_par,2}(1,:);
    
    for s = 1:length(session)
        
        % Define runs within this block (but this time based on all present in folder)        
        run_files =dir([root_dir '/Events/complete_event_log/' participant{p}(1:4) session{s} '/'...
        participant{p} session{s} '_Click*_tsss.vmrk']);

        runs_in_block = {};
        for i = 1:length(run_files)
            end_run_name = regexp(run_files(i).name, '_mc_tsss', 'once');
            runs_in_block{1,i} = run_files(i).name(7:end_run_name);
        end

        for r = 1:length(runs_in_block)

            folders = dir([root_dir '/Events/complete_event_log/' participant{p}(1:4) session{s} '/']);
            results =  contains({folders.name},runs_in_block{r}) & endsWith({folders.name},'vmrk');
            infolder = find(results);
                if isempty(infolder)
                    continue
                end
                if size(infolder,2)> 1 % in case of more than one coincidence, vmrk is going to be the same for all, so pick first
                    infolder = infolder(1);
                end
            filename = folders(infolder).name;

            fullfilename = [root_dir '/Events/complete_event_log/' participant{p}(1:4) session{s} '/' filename];
            if ~exist(fullfilename, 'file')
                disp (['filename ' fullfilename ' does not exist'])
                continue
            end
            fileID = fopen(fullfilename,'r');
            textscan(fileID, '%[^\n\r]', startRow-1, 'WhiteSpace', '', 'ReturnOnError', false, 'EndOfLine', '\r\n');
            dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'ReturnOnError', false);
            fclose(fileID);
            dataArray_full = [];
            dataArray_full(:,1) = dataArray{1,1};
            dataArray_full(:,2) = dataArray{1,2};
%             dataArray_full = str2double(dataArray_full);
            save([root_dir '/Events/complete_event_log/' participant{p} session{s} '/' participant{p} session{s} '_' runs_in_block{r} '.mat'],'dataArray_full');
            clearvars filename fileID dataArray ans;
        end
    end    
end

clearvars('-except', initialVars{:});
disp 'OBTAINING MATLAB COMPLETE LOG (TO FILL GAPS IN FINAL SINGLE TRIAL MATRIX)!!!'
disp(datetime)
toc

%% Obtain averaged SINGLE TRIAL amplitudes (LLR)

tic
disp(' ');      
disp('-------------------------');  
disp('OBTAINING DATA TABLE FOR FOR SINGLE TRIAL AMPLITUDES (LLR)');  
disp(datetime)
disp('-------------------------');     
disp(' ');

% If LLR seg is not the protocol at this step, better stop
if ~strcmp(root_dir_bs,'~/brainstorm_db/Project_LLRseg')
    error('root_dir_bs is not LLRseg, change it!')
end

% Separate EEG and MEG files because they have different rejection criteria
modality_data = {'EEG','MEG'};

% Same for EEG and MEG (using same ones than for EEG statistics at FCz: MEG ones are similar)
components_LLR = {'N1' 'P2'};
Time_wind_N1 = [0.078 0.106]; % Has to be in seconds 
Time_wind_P2 = [0.150 0.182];

Time_wind_N1_MEG = [0.078 0.106]; % Has to be in seconds 
Time_wind_P2_MEG = [0.134 0.166];

% THESE WILL BE THE SAME FOR EEG AND MEG, SO REMOVE '_MEG' FOR MLR
% components_MLR = {'P0' 'Na' 'Pa' 'Nb' 'P50'};
% Time_wind_P0 = [0.005 0.01133]; % 8 +-3 %
% Time_wind_Na = [0.015 0.02133]; % 18 +-3
% Time_wind_Pa = [0.027 0.0333]; % 30 +-3
% Time_wind_Nb = [0.038 0.0446]; % 41 +-3
% Time_wind_P50 = [0.051 0.06133]; % 56 ms +- 5

Single_trial_amplitudes_EEG_LLR = {};
Single_trial_amplitudes_MEG_LLR = {};
% Single_trial_amplitudes_EEG_MLR = {};
% Single_trial_amplitudes_MEG_MLR = {};
% Load channel data
load ([root_dir '/Brainstorm_pipelines/Areas_channels.mat'])
% colNames_EEG_MLR = {'TrialNr','subjID','sessionNr','blockNr','run','absTime','clickIntensity','ISI',...
%     'EEG_P0','EEG_Na', 'EEG_Pa', 'EEG_Nb', 'EEG_P50'};
% colNames_MEG_MLR = {'TrialNr','subjID','sessionNr','blockNr','run','absTime','clickIntensity','ISI',...
%     'MEG_P0', 'MEG_Na', 'MEG_Pa', 'MEG_Nb', 'MEG_P50'};
colNames_EEG_LLR = {'TrialNr','subjID','sessionNr','blockNr','run','absTime','clickIntensity','ISI',...
    'EEG_N1', 'EEG_P2'};
colNames_MEG_LLR = {'TrialNr','subjID','sessionNr','blockNr','run','absTime','clickIntensity','ISI',...
    'MEG_N1', 'MEG_P2'};
% EEG/MEG channels
% In case we wanna do groups of sensors
% chan_list_EEG = find((contains({Areas.Group}, 'Frontal') | contains({Areas.Group}, 'Central') |...
%     contains({Areas.Group}, 'Fronto_Central')) & strcmp({Areas.Type}, 'EEG') ==1);
% chan_list_MEG_MAG = find((contains({Areas.Group}, 'Front') | contains({Areas.Group}, 'Pariet'))...
%     & strcmp({Areas.Type}, 'MEG MAG') ==1);
chan_list_EEG = find(strcmp({Areas.Name}, 'FCz'));
chan_list_MEG_MAG = find(strcmp({Areas.Name}, 'MEG0131')); % Left

for p = 1:length(participant)
    eval(['load ''' root_dir '/Events/Single_trials/Latencies_' participant{p} '.mat'''])
    
    % Define sessions within this participant
    pos_par = find(strcmp(session_block_array(:,1), participant{p}(1:4)));
    session = session_block_array{pos_par,2}(1,:);
    
    for s = 1:length(session)
        
        % Define blocks within this session
        pos_ses = find(strcmp(session_block_array{pos_par,2}(1,:), session{s}));
        block = session_block_array{pos_par,2}{2,pos_ses}(1,:);
        
        for b = 1:length(block)
            for mode = 1:length(modality_data) % EEG and MEG
                for c = 1:length(condition)
                eval(['ncs_' participant{p} '_' modality_data{mode} '_' condition{c} session{s} '_' block{b} '_LLR = 0;'])    
                folders = dir([root_dir_bs '/data/' participant{p} '/']); % Will be LLRseg here
                results = contains({folders.name},['LLR_' condition{c} session{s} '_' block{b} '_normal']);
                infolder = find(results);
                    if isempty(infolder)
                        continue % should not be a problem unless all conditions are bad for a block, but in this case that block is not probably in the block list
                    end
                    dir_trial = dir([root_dir_bs '/data/' participant{p} '/' folders(infolder).name]);
                    trial_cov = contains({dir_trial.name},'_trial');
                    trial_count = find(trial_cov);
                    
                    if isempty(trial_count)
                        continue % this is a rare exception, so just in case
                    end
                    
                    % Load bad trials (EEG or MEG)
                    if ~exist([root_dir '/Events/BadTrials/LLR/' modality_data{mode} '/' participant{p} '/' folders(infolder).name '/BadTrials.mat'], 'file')
                        continue % if this folder does not have bad trials corresponding, better discard these sweeps
                    end
                    load([root_dir '/Events/BadTrials/LLR/' modality_data{mode} '/' participant{p} '/' folders(infolder).name '/BadTrials.mat']);
                    
                    disp(' ');      
                    disp('-------------------------');  
                    disp(['Extracting LLR single trial amplitudes ' modality_data{mode} '_' participant{p} '_' session{s} '_' block{b} '_' condition{c}]);
                    disp(' ');  
                    
                    for j= 1:length(trial_count) % for every trial
                        position = trial_count(j);
                        trial_name = [root_dir_bs '/data/' participant{p} '/' folders(infolder).name '/' dir_trial(position).name];
                                                                            
                        eval(['load ' trial_name]) % trial loaded 
                        % contains(ChannelFlag ColormapType Comment DataType Device DisplayUnits Events F History Leff Std Time nAvg)
                        % Retrieve value and run number of each trial
                        [x,y]= find(contains(History,'import_time'));
                        if isempty(x); continue; end % Just in case
                        
                        brack = regexp(History{x,y+1},'[');com = regexp(History{x,y+1},',');      
                        % latency absolute value (for match with log data)
                        trial_time = str2double(History{x,y+1}(brack +1: com -1));
                        trial_time = trial_time + 0.150; % For LLR (+ 0.05 for MLR)               
                     
                        % run number
                        cell_import = find(contains(History,'Import from:'));
                        if isempty(cell_import); continue; end % just in case
                        init_run_name = regexp(History{cell_import}, 'run', 'once');
                        end_run_name = regexp(History{cell_import}, '_mc_tsss', 'once');
                        trial_run = History{cell_import}(init_run_name+3:end_run_name-1); % gives the two digit number of the run
                        % Check if any Log file is empty before looking into it
                        eval(['empty_logs = isempty(Log_' condition{c} '_' participant{p} session{s} '_r' trial_run ');'])
                        if empty_logs == 1 % if the log is empty
                            continue
                        end
                        % find position in Log      
                        eval(['ident_pos = find(abs(trial_time  - Log_' condition{c} '_' participant{p} session{s} '_r' trial_run '(:,1))<10^-2);']) % could be even one
                        if isempty(ident_pos)
                            eval(['ncs_' participant{p} '_' modality_data{mode} '_' condition{c} session{s} '_' block{b} '_LLR = ncs_' participant{p} '_' modality_data{mode} '_' condition{c} session{s} '_' block{b} '_LLR + 1;'])
                            continue
                        end
                        % Previous category info
                        eval(['prev_cat = Log_' condition{c} '_' participant{p} session{s} '_r' trial_run '(ident_pos,3);'])
                        % avoid rare random categories that I don't expect (like '8' or '' or 'boundary'), I rather do this
                        if (prev_cat ~= 11) && (prev_cat ~= 12) && (prev_cat ~= 13) && ...
                            (prev_cat ~= 31) && (prev_cat ~= 32) && (prev_cat ~= 33) &&...
                            (prev_cat ~= 51) && (prev_cat ~= 52) && (prev_cat ~= 53) &&...
                            (prev_cat ~= 71) && (prev_cat ~= 72) && (prev_cat ~= 73) &&...
                            (prev_cat ~= 91) && (prev_cat ~= 92) && (prev_cat ~= 93)
                            % The way the logs are programmed now this
                            % should be impossible, left this just in case
                            continue
                        end
                        eval(['prev_lat = Log_' condition{c} '_' participant{p} session{s} '_r' trial_run '(ident_pos,2);'])
                                                
                        % save trial amplitudes before overwriting F variable
                        trial_amplitudes = F; % always gonna be 392 channels
                        length_trial = size(trial_amplitudes,2);
                        clear F ChannelFlag ColormapType Comment DataType Device DisplayUnits Events F History Leff Std Time nAvg
                        
                        baseline = 0.15; 
                        % baseline = 0.05; FOR MLR!!!!
                        if strcmp(modality_data{mode}, 'EEG')
                            for comp_LLR = 1:length(components_LLR)
                                eval(['EEG_' components_LLR{comp_LLR} '= trial_amplitudes(chan_list_EEG,fix(((baseline*resample_LLR) + (Time_wind_' components_LLR{comp_LLR} '(1)*resample_LLR))):fix(((baseline*resample_LLR) + (Time_wind_' components_LLR{comp_LLR} '(2)*resample_LLR))));']) 
                                eval(['EEG_' components_LLR{comp_LLR} '= EEG_' components_LLR{comp_LLR} '.*1000000;']) % Transform from Volts to microVolts
                                % Next two columns lines in case we wanna do average across sensors
                                % Delete any empty channels before average
%                                 eval(['EEG_' components_LLR{comp_LLR}  '= EEG_' components_LLR{comp_LLR}  '(any(EEG_' components_LLR{comp_LLR} ',2),:);'])
%                                 eval(['EEG_' components_LLR{comp_LLR} '= mean(EEG_' components_LLR{comp_LLR} ',1);']) 
                                eval(['EEG_' components_LLR{comp_LLR} '= mean(EEG_' components_LLR{comp_LLR} ');'])
                            end
                            
                            if isempty(Single_trial_amplitudes_EEG_LLR)                             
                                row_pos_EEG = 1;
                            else
                                row_pos_EEG = size(Single_trial_amplitudes_EEG_LLR,1) + 1;
                            end  
                            
                            %%% IF it's a bad trial, rewrite component
                            %%% amplitudes as 'NaN', as Tobias wanted
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            if ~isempty(BadTrials)   
                                if any(contains(BadTrials, [dir_trial(position).name]))
                                    for comp_LLR = 1:length(components_LLR) % N1, P2, etc
                                        eval(['EEG_' components_LLR{comp_LLR} '= ''NaN'';'])
                                    end
                                end
                            end
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            
                            Single_trial_amplitudes_EEG_LLR{row_pos_EEG,1} = j;
                            Single_trial_amplitudes_EEG_LLR{row_pos_EEG,2} = participant{p};  
                            Single_trial_amplitudes_EEG_LLR{row_pos_EEG,3} = session{s};    
                            Single_trial_amplitudes_EEG_LLR{row_pos_EEG,4} = block{b};
                            Single_trial_amplitudes_EEG_LLR{row_pos_EEG,5} = str2double(trial_run);
                            Single_trial_amplitudes_EEG_LLR{row_pos_EEG,6} = trial_time;
                            Single_trial_amplitudes_EEG_LLR{row_pos_EEG,7} = str2double(condition{c}(2)); % 1 (65dB), 2 (75dB), 3 (85dB)
                            Single_trial_amplitudes_EEG_LLR{row_pos_EEG,8} = prev_lat;
                            Single_trial_amplitudes_EEG_LLR{row_pos_EEG,9} = EEG_N1;
                            Single_trial_amplitudes_EEG_LLR{row_pos_EEG,10} = EEG_P2;
                            
                        elseif strcmp(modality_data{mode}, 'MEG')
                            for comp_LLR = 1:length(components_LLR) 
                                eval(['Meg_Mag_' components_LLR{comp_LLR} '= trial_amplitudes(chan_list_MEG_MAG,fix(((baseline*resample_LLR) + (Time_wind_' components_LLR{comp_LLR} '_MEG(1)*resample_LLR))):fix(((baseline*resample_LLR) + (Time_wind_' components_LLR{comp_LLR} '_MEG(2)*resample_LLR))));']) 
                                eval(['Meg_Mag_' components_LLR{comp_LLR} '= Meg_Mag_' components_LLR{comp_LLR} '.*1000000000000000;']) % Tesla to fT
                                % Next two columns lines in case we wanna do average across sensors
                                % Delete any empty channels before average
%                                 eval(['Meg_Mag_' components_LLR{comp_LLR}  '= Meg_Mag_' components_LLR{comp_LLR}  '(any(Meg_Mag_' components_LLR{comp_LLR} ',2),:);'])
%                                 eval(['Meg_Mag_' components_LLR{comp_LLR} '= mean(Meg_Mag_' components_LLR{comp_LLR} ',1);']) 
                                eval(['Meg_Mag_' components_LLR{comp_LLR} '= mean(Meg_Mag_' components_LLR{comp_LLR} ');'])
                            end
                            
%                           row_pos = ((s-1)*45) + ((b-1)*15) + c; % index position in table 
                            if isempty(Single_trial_amplitudes_MEG_LLR)                             
                                row_pos_MEG = 1;
                            else
                                row_pos_MEG = size(Single_trial_amplitudes_MEG_LLR,1) + 1;
                            end
                            
                            %%% IF it's a bad trial, rewrite component
                            %%% amplitudes as 'NaN', as Tobias wanted
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            if ~isempty(BadTrials)   
                                if any(contains(BadTrials, [dir_trial(position).name]))
                                    for comp_LLR = 1:length(components_LLR) % N1, P2, etc
%                                         eval(['Meg_Grad_' components_LLR{comp_LLR} '= ''NaN'';'])
                                        eval(['Meg_Mag_' components_LLR{comp_LLR} '= ''NaN'';'])
                                    end
                                end
                            end
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            
                            Single_trial_amplitudes_MEG_LLR{row_pos_MEG,1} = j;
                            Single_trial_amplitudes_MEG_LLR{row_pos_MEG,2} = participant{p};  
                            Single_trial_amplitudes_MEG_LLR{row_pos_MEG,3} = session{s};    
                            Single_trial_amplitudes_MEG_LLR{row_pos_MEG,4} = block{b};
                            Single_trial_amplitudes_MEG_LLR{row_pos_MEG,5} = str2double(trial_run);
                            Single_trial_amplitudes_MEG_LLR{row_pos_MEG,6} = trial_time;
                            Single_trial_amplitudes_MEG_LLR{row_pos_MEG,7} = str2double(condition{c}(2)); % 1 (65dB), 2 (75dB), 3 (85dB)
                            Single_trial_amplitudes_MEG_LLR{row_pos_MEG,8} = prev_lat;
                            Single_trial_amplitudes_MEG_LLR{row_pos_MEG,9} = Meg_Mag_N1;
                            Single_trial_amplitudes_MEG_LLR{row_pos_MEG,10} = Meg_Mag_P2;
                        end   
                    end
                end
            end
            %%%%%%%%%%%%%%%% EEG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % find rows of this block (just for this session and subject)
            block_rows_EEG = find(strcmp({Single_trial_amplitudes_EEG_LLR{:,2}}, participant{p}) ...
            & strcmp({Single_trial_amplitudes_EEG_LLR{:,3}}, session{s})...
             & strcmp({Single_trial_amplitudes_EEG_LLR{:,4}}, block{b})); 
            % Delimite current block in matrix
            current_block_EEG = Single_trial_amplitudes_EEG_LLR(block_rows_EEG,:);    
            current_block_EEG = sortrows(current_block_EEG,5); % sort runs
            psbl_run_values_EEG = unique([current_block_EEG{:,5}]); % identify run numbers
            % For every run
            for k = 1:length(psbl_run_values_EEG) 
                run_rows_EEG = find(cell2mat(current_block_EEG(:,5)) == psbl_run_values_EEG(k));
                current_run_EEG = current_block_EEG(run_rows_EEG,:);
                current_run_EEG = sortrows(current_run_EEG,6); % sort sweep times
                current_run_EEG(:,1) = num2cell(1:size(current_run_EEG)); % correct sweep number list
                current_block_EEG(run_rows_EEG,:) = current_run_EEG; % take run back to block section it belongs    
            end
            Single_trial_amplitudes_EEG_LLR(block_rows_EEG,:) = current_block_EEG;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%% MEG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % find rows of this block (just for this session and subject)
            block_rows_MEG = find(strcmp({Single_trial_amplitudes_MEG_LLR{:,2}}, participant{p}) ...
            & strcmp({Single_trial_amplitudes_MEG_LLR{:,3}}, session{s})...
             & strcmp({Single_trial_amplitudes_MEG_LLR{:,4}}, block{b})); 
            % Delimite current block in matrix
            current_block_MEG = Single_trial_amplitudes_MEG_LLR(block_rows_MEG,:);    
            current_block_MEG = sortrows(current_block_MEG,5); % sort runs
            psbl_run_values_MEG = unique([current_block_MEG{:,5}]); % identify run numbers
            % For every run
            for k = 1:length(psbl_run_values_MEG) 
                run_rows_MEG = find(cell2mat(current_block_MEG(:,5)) == psbl_run_values_MEG(k));
                current_run_MEG = current_block_MEG(run_rows_MEG,:);
                current_run_MEG = sortrows(current_run_MEG,6); % sort sweep times
                current_run_MEG(:,1) = num2cell(1:size(current_run_MEG)); % correct sweep number list
                current_block_MEG(run_rows_MEG,:) = current_run_MEG; % take run back to block section it belongs    
            end
            Single_trial_amplitudes_MEG_LLR(block_rows_MEG,:) = current_block_MEG;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
    end
    save([root_dir '/Events/Single_trials/Non_coincident_sweeps_' participant{p} '.mat'],'-regexp','^ncs_');
end

Amp_singl_trials_EEG_LLR = array2table(Single_trial_amplitudes_EEG_LLR,'VariableNames',colNames_EEG_LLR);
Amp_singl_trials_MEG_LLR = array2table(Single_trial_amplitudes_MEG_LLR,'VariableNames',colNames_MEG_LLR);
writetable(Amp_singl_trials_EEG_LLR, [root_dir '/Events/Single_trials/Amp_singl_trials_EEG_LLR.txt'], 'Delimiter','\t')
writetable(Amp_singl_trials_MEG_LLR, [root_dir '/Events/Single_trials/Amp_singl_trials_MEG_LLR.txt'], 'Delimiter','\t')

% Reset modality_data variable to its original self
modality_data = {'EEG','MEG','both_mod'};

clearvars('-except', initialVars{:});
disp 'DONE WITH OBTAINING DATA TABLE FOR FOR SINGLE TRIAL AMPLITUDES (LLR)!!!'
disp(datetime)
toc

%% Fill table with missing trials for LLR (rejected at preprocessing or lost)

tic
disp(' ');      
disp('-------------------------');  
disp('FILL SINGLE TRIAL TABLE WITH MISSING TRIALS (LLR)');  
disp(datetime)
disp('-------------------------');     
disp(' ');

If LLR seg is not the protocol at this step, better stop
if ~strcmp(root_dir_bs,'~/brainstorm_db/Project_LLRseg')
    error('root_dir_bs is not LLRseg, change it!')
end

add_runs = 1; % For empty runs
Add_runs = {};

% Separate EEG and MEG files because they have different rejection criteria
modality_data = {'EEG','MEG'};
% modality_data = {'EEG'};
components_LLR = {'N1' 'P2'};
% components_MLR = {'P0' 'Na' 'Pa' 'Nb' 'P50'};

for mode = 1:length(modality_data) % EEG and MEG
    
Amp_singl_trials_LLR = readtable([root_dir '/Events/Single_trials/Amp_singl_trials_' modality_data{mode} '_LLR.txt']);
for p = 1:length(participant)
    
    % Define sessions within this participant
    pos_par = find(strcmp(session_block_array(:,1), participant{p}(1:4)));
    session = session_block_array{pos_par,2}(1,:);
    for s = 1:length(session)

        % Define runs within this block (but this time based on all present in folder)        
        run_files =dir([root_dir '/Events/complete_event_log/' participant{p}(1:4) session{s} '/'...
        participant{p} session{s} '_Click*.mat']);

        runs_in_block = {};
        for i = 1:length(run_files)
            end_run_name = regexp(run_files(i).name, '_.mat', 'once');
            runs_in_block{1,i} = run_files(i).name(7:end_run_name);
        end
        runs_10 = find(contains(runs_in_block,'run10_'));
        runs_11 = find(contains(runs_in_block,'run11_'));
        runs_12 = find(contains(runs_in_block,'run12_'));
        if ~isempty(runs_12); error('more than 11 blocks');end
        runs_13 = find(contains(runs_in_block,'run13_'));
        if ~isempty(runs_13); error('more than 12 blocks');end
        runs_14 = find(contains(runs_in_block,'run14_'));
        if ~isempty(runs_14); error('more than 13 blocks');end
        
        % Order blocks in case there are 11
        if ~isempty(runs_11) % there are just 11 runs, so order them
            name_block_11 = runs_in_block(runs_11);
            name_block_10 = runs_in_block(runs_10);
            run_matrix_size = length(runs_in_block);
            runs_in_blockCopy = runs_in_block;
            runs_in_blockCopy(1:run_matrix_size-2) = runs_in_blockCopy(3:run_matrix_size);
            runs_in_blockCopy(run_matrix_size) = name_block_11;
            runs_in_blockCopy(run_matrix_size-1) = name_block_10;
            clear('runs_in_block');
            runs_in_block = runs_in_blockCopy;
        % Order blocks in case there are 10
        elseif ~isempty(runs_10) % there are just 10 runs, so order them
            name_block_10 = runs_in_block(runs_10);
            run_matrix_size = length(runs_in_block);
            runs_in_blockCopy = runs_in_block;
            runs_in_blockCopy(1:run_matrix_size-1) = runs_in_blockCopy(2:run_matrix_size);
            runs_in_blockCopy(run_matrix_size) = name_block_10;
            clear('runs_in_block');
            runs_in_block = runs_in_blockCopy;
        end
        
        for r = 1:length(runs_in_block)
            
            disp(' ');      
            disp('-------------------------');  
            disp(['Filling single trial matrix  ' modality_data{mode} '_' participant{p} '_' session{s} '_' runs_in_block{r} 'LLR']);
            disp(' '); 
            
            % Load dataArray_full
            load ([root_dir '/Events/complete_event_log/' participant{p}(1:4) session{s} '/'...
            participant{p} session{s} '_' runs_in_block{r} '.mat']);
            % Remove rows with NaN 
            dataArray_full = rmmissing(dataArray_full);
            
            % Retrieve number of run:
            % should be the same than r, but could not, so better be sure
            init_run_name = regexp(runs_in_block{r}, 'run', 'once');
            number_of_run = str2double(runs_in_block{r}(init_run_name+3:end-1)); % gives the two digit number of the run
        
            % Rows with current participant
            rp = find(contains(Amp_singl_trials_LLR.subjID,participant{p}));
            % Rows with current session WITHIN current subject coincidences
            rs = find(contains(Amp_singl_trials_LLR.sessionNr(rp),session{s})); %#ok<*FNDSB>
            rs = rp(rs); 
            if isempty(rs); continue; end % almost impossible, but just in case
            % Rows with current run WITHIN current session coincidences
            rr = find(Amp_singl_trials_LLR.run(rs) == number_of_run);
            rr = rs(rr);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%% Modified from original %%%%%%%%%%%%%%%
            % If the run file existed in complete_event_log but not in data
            % matrix, likely it's a run that was deleted for whatever
            % reason. In that case, we need to add the sweeps to preserve
            % the integrity of the 20m blocks, and then go to next block
            % Since it's extremely rare and hard to do with script (we do not have the block)
            % we will do it manually at the end, for which we need to know which ones to add
            if isempty(rr)
                Add_runs{add_runs,1} = [participant{p} '_' session{s} '_run_' num2str(number_of_run)];
                add_runs = add_runs+1;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if number_of_run == 1 && s == 1 % first run and session of this subject
                    % Then it goes after the last subject (we know it does
                    % not happen in the first subject first run)
                    rp_1 = find(contains(Amp_singl_trials_LLR.subjID,participant{p-1}));
                    row_add_run = rp_1(end);
                    block = 'b1'; % exception because if the whole block is
                    % not in session_block array it's gonna be empty
                elseif number_of_run == 1 && s ~= 1 % so there is a previous session to check
                    % Then it goes after the last session of present subject
                    rs_1 = find(contains(Amp_singl_trials_LLR.sessionNr(rp),session{s-1})); %#ok<*FNDSB>
                    rs_1 = rp(rs_1); 
                    row_add_run = rs_1(end);
                    block = 'b1'; % exception because if the whole block is
                    % not in session_block array it's gonna 
                else % This is not the first run, so use last as reference
                    rr_1 = find(Amp_singl_trials_LLR.run(rs) == number_of_run-1);
                    rr_1 = rs(rr_1);
                    row_add_run = rr_1(end);
                    % Find out in which block is this run for this participant and session
                    pos_ses = find(strcmp(session_block_array{pos_par,2}(1,:), session{s}));
                    array_search_run = session_block_array{pos_par,2}{2,pos_ses};
                    matrix_coinc = [];
                    for sar = 1:size(array_search_run,2)
                        current_run_pos = find(contains(array_search_run{2,sar}, ['_run' num2str(number_of_run) '_']));
                        if ~isempty(current_run_pos)
                            current_block = array_search_run(1,sar);
                            current_block = current_block{1};
                        end
                    end
                    if exist('current_block','var')
                        % so a coincidence was found in session_block_array
                        % and current_block variable was created
                        block = current_block;
                        clear('current_block'); % for next time
                    else % Find it using standard rules or manual specifications
                        if number_of_run == 1 || number_of_run == 2 || number_of_run == 3
                            block = 'b1';
                        elseif number_of_run == 4 || number_of_run == 5 || number_of_run == 6
                            block = 'b2';
                        else % 7 to 10
                            block = 'b3';
                        end
                        if strcmp(participant{p},'2193') && strcmp(session{s},'C')
                            if number_of_run == 5 || number_of_run == 6 || number_of_run == 7
                                block = 'b2';
                            elseif number_of_run == 8 || number_of_run == 9 || number_of_run == 10
                                block = 'b3';
                            end
                        end
                        % Exceptions:
%                         2193C 5 to 7 are b2
%                         2193C 8 to 10 are b3
%                         H001_D_run_10 should be block 3
%                         previous_block = Amp_singl_trials_LLR.blockNr{row_add_run};
%                         previous_block = str2double(previous_block(2));
%                         previous_block = previous_block +1;
%                         block = ['b' num2str(previous_block)];
                    end   
                end
                % To ensure it gets only the number of rows (trials)
                rows_dataArray = size(dataArray_full,1);
                for i = 1:rows_dataArray
                    % First, duplicate row as many times as trials in run
                    Amp_singl_trials_LLR = Amp_singl_trials_LLR([1:row_add_run,row_add_run:end],:);
                end
                for i = 1:rows_dataArray
                    % Second, add the info in each of created rows
                    Amp_singl_trials_LLR.TrialNr(row_add_run+i) = i;
                    Amp_singl_trials_LLR.subjID{row_add_run+i} = participant{p}; % not needed
                    Amp_singl_trials_LLR.sessionNr{row_add_run+i} = session{s};
                    Amp_singl_trials_LLR.blockNr{row_add_run+i} = block;
                    Amp_singl_trials_LLR.run(row_add_run+i) = number_of_run;
                    curr_abs_time = dataArray_full(i,2)/1500; % from samples to s
                    Amp_singl_trials_LLR.absTime(row_add_run+i) = curr_abs_time;

                    % Determine dB based on trigger code of current value:
                    current_dB = num2str(dataArray_full(i,1));
                    try
                        current_dB = str2double(current_dB(2));
                    catch
                        current_dB = NaN;
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    Amp_singl_trials_LLR.clickIntensity(row_add_run+i) = current_dB;

                    % Determine ISI using current and previous abs time
                    % provided it's not first trial of run (in which case we
                    % cannot compare with abs times from previous run)
                    if i == 1 % it's the first trial from this run
                        Amp_singl_trials_LLR.ISI(row_add_run+i)= NaN;
                    else
                        current_ISI = curr_abs_time - Amp_singl_trials_LLR.absTime(row_add_run+i-1);
                        % Add value into the matrix
                        Amp_singl_trials_LLR.ISI(row_add_run+i) = current_ISI; % in s
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    % Add NaN in amplitude values
                    for comp_LLR = 1:length(components_LLR) % Will do for every column of components
                        eval(['Amp_singl_trials_LLR.' modality_data{mode} '_' components_LLR{comp_LLR} '(row_add_run+i)= NaN;'])
                    end
                end
                continue
            end
            % Happens with H001Cb1 and H002Cb1 for instance
            
            % For every trial found in run file
            for i = 1:length(dataArray_full)
                curr_abs_time = dataArray_full(i,2)/1500; % from samples to s
                % Now find that value from complete log in the single trial matrix:
                rv = find(abs(curr_abs_time - Amp_singl_trials_LLR.absTime(rr))<10^-2); % could be one if found many inconsistencies
                % If that value IS in matrix, there is nothing to add
                if ~isempty(rv); continue;end
                % If it's not, let's add it:
                % Find the row position with the closest timing to the trial to add
                [minValue,closestIndex] = min(abs(curr_abs_time - Amp_singl_trials_LLR.absTime(rr)));
                % Now is the value to add bigger or smaller than the closest coincidence in the matrix
                % It is important to do so now before modifying the size of the matrix
                % Adjust for the closest index IN THE ROWS of this run
                closestIndex = rr(closestIndex);
                if  curr_abs_time > Amp_singl_trials_LLR.absTime(closestIndex)
                    % Our value goes after the closest index, so, upon
                    % duplicating that row, the second row of these two is
                    % the one that should have the added (bigger) value
                    position = closestIndex + 1;
                else
                    % Our value goes before the closest index, so, upon
                    % duplicating that row, the first row of these two is
                    % the one that should have the added (smaller) value
                    position = closestIndex;
                end
                
                % Duplicate the row of that closet Index
                Amp_singl_trials_LLR = Amp_singl_trials_LLR([1:closestIndex,closestIndex:end],:); 
                % Now in the right position 
                
                Amp_singl_trials_LLR.TrialNr(position) = 999; % Since it will be corrected later
                % Not needed to update subject, session, block or run since they should be the same
                Amp_singl_trials_LLR.absTime(position) = curr_abs_time; % in s
                % Determine dB based on trigger code of current value:
                current_dB = num2str(dataArray_full(i,1));
                try
                    current_dB = str2double(current_dB(2));
                catch
                    current_dB = NaN;
                end
                % Add value into the matrix
                Amp_singl_trials_LLR.clickIntensity(position) = current_dB; 
                % Determine ISI using current and previous abs time
                % provided it's not first trial of run (in which case we
                % cannot compare with abs times from previous run)
                if i == 1 % it's the first trial from this run
                    Amp_singl_trials_LLR.ISI(position)= NaN;
                else
                    current_ISI = curr_abs_time - Amp_singl_trials_LLR.absTime(position-1);
                    % Add value into the matrix
                    Amp_singl_trials_LLR.ISI(position) = current_ISI; % in s
                end
                % Add NaN in amplitude values
                for comp_LLR = 1:length(components_LLR) % Will do for every column of components
                    eval(['Amp_singl_trials_LLR.' modality_data{mode} '_' components_LLR{comp_LLR} '(position)= NaN;'])
                end
                
                % After new row has been added with all new info:
                %%%%%% Gotta update these every time a row is included %%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Rows with current participant
                rp = find(contains(Amp_singl_trials_LLR.subjID,participant{p}));
                % Rows with current session WITHIN current subject coincidences
                rs = find(contains(Amp_singl_trials_LLR.sessionNr(rp),session{s})); %#ok<*FNDSB>
                rs = rp(rs); 
                if isempty(rs); continue; end % almost impossible, but just in case
                % Rows with current run WITHIN current session coincidences
                rr = find(Amp_singl_trials_LLR.run(rs) == number_of_run);
                rr = rs(rr);
                if isempty(rr); continue; end % incredibly rare, possibly only happens with H001Cb1 
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            end
            % Correct TrialNr after all the additions (1:end in each run):
            % REDEFINE rows with subject/session/run (NOW THEY MAY BE DIFFERENT IF WE ADDED ROWS)
            % Rows with current participant
            rp = find(contains(Amp_singl_trials_LLR.subjID,participant{p}));
            % Rows with current session WITHIN current subject coincidences
            rs = find(contains(Amp_singl_trials_LLR.sessionNr(rp),session{s})); %#ok<*FNDSB>
            rs = rp(rs);
            if isempty(rs); continue; end % almost impossible, but just in case
            % Rows with current run WITHIN current session coincidences
            rr = find(Amp_singl_trials_LLR.run(rs) == number_of_run);
            rr = rs(rr);
            if isempty(rr); continue; end % incredibly rare, possibly only happens with H001Cb1 
            % Now correct TrialNr column numbers      
            Amp_singl_trials_LLR.TrialNr(rr) = 1:length(Amp_singl_trials_LLR.TrialNr(rr));   
        end    
    end 
end

if strcmp(modality_data{mode},'EEG')
    Amp_singl_trials_EEG_LLR = Amp_singl_trials_LLR;
    writetable(Amp_singl_trials_EEG_LLR, [root_dir '/Events/Single_trials/Amp_singl_trials_EEG_LLR.txt'])
    clear('Amp_singl_trials_EEG_LLR'); % to free space
elseif strcmp(modality_data{mode},'MEG')
    Amp_singl_trials_MEG_LLR = Amp_singl_trials_LLR;
    writetable(Amp_singl_trials_MEG_LLR, [root_dir '/Events/Single_trials/Amp_singl_trials_MEG_LLR.txt'])
    clear('Amp_singl_trials_MEG_LLR'); % to free space
end

save([root_dir '/Events/Single_trials/Add_runs.mat'],'Add_runs');

end

% Reset modality_data variable to its original self
modality_data = {'EEG','MEG','both_mod'};

clearvars('-except', initialVars{:});
disp 'DONE WITH FILLING SINGLE TRIAL TABLE WITH MISSING TRIALS (LLR)!!!'
disp(datetime)
toc

%% Fix three sessions with inverted click intensities LLR

tic
disp(' ');      
disp('-------------------------');  
disp('FIX THREE SESSIONS WITH INVERTED CLICK INTENSITIES (LLR)');  
disp(datetime)
disp('-------------------------');     
disp(' ');

% Separate EEG and MEG files because they have different rejection criteria
modality_data = {'EEG'};

for mode = 1:length(modality_data) % EEG and MEG

    if strcmp(modality_data{mode},'EEG')
        error('This was already ran for EEG, running it again will flip these sessions to the original (wrong) state');
    end
%%%%%%%%%%%%%%%%%%%%% FIRST, IDENTIFY ROWS TO MODIFY %%%%%%%%%%%%%%%%%%
    
Amp_singl_trials_LLR = readtable([root_dir '/Events/Single_trials/Amp_singl_trials_' modality_data{mode} '_LLR.txt']);

% Find 2005D louder
A = ismember(Amp_singl_trials_LLR.subjID,'2005');
B = ismember(Amp_singl_trials_LLR.sessionNr,'D');
C = find(Amp_singl_trials_LLR.clickIntensity == 3);
count = 1;
for i = 1:size(Amp_singl_trials_LLR,1)
    if A(i)&& B(i) && sum(i == C)
    coincident_rows_2005D_louder(count) = i;
    count = count +1; 
    end
end

% Find 2005D quieter
A = ismember(Amp_singl_trials_LLR.subjID,'2005');
B = ismember(Amp_singl_trials_LLR.sessionNr,'D');
C = find(Amp_singl_trials_LLR.clickIntensity == 1);
count = 1;
for i = 1:size(Amp_singl_trials_LLR,1)
    if A(i)&& B(i) && sum(i == C)
    coincident_rows_2005D_quieter(count) = i;
    count = count +1; 
    end
end

% Find 2238G louder
A = ismember(Amp_singl_trials_LLR.subjID,'2238');
B = ismember(Amp_singl_trials_LLR.sessionNr,'G');
C = find(Amp_singl_trials_LLR.clickIntensity == 3);
count = 1;
for i = 1:size(Amp_singl_trials_LLR,1)
    if A(i)&& B(i) && sum(i == C)
    coincident_rows_2238G_louder(count) = i;
    count = count +1; 
    end
end

% Find 2238G quieter
A = ismember(Amp_singl_trials_LLR.subjID,'2238');
B = ismember(Amp_singl_trials_LLR.sessionNr,'G');
C = find(Amp_singl_trials_LLR.clickIntensity == 1);
count = 1;
for i = 1:size(Amp_singl_trials_LLR,1)
    if A(i)&& B(i) && sum(i == C)
    coincident_rows_2238G_quieter(count) = i;
    count = count +1; 
    end
end

% Find 2235H louder
A = ismember(Amp_singl_trials_LLR.subjID,'2235');
B = ismember(Amp_singl_trials_LLR.sessionNr,'H');
C = find(Amp_singl_trials_LLR.clickIntensity == 3);
count = 1;
for i = 1:size(Amp_singl_trials_LLR,1)
    if A(i)&& B(i) && sum(i == C)
    coincident_rows_2235H_louder(count) = i;
    count = count +1; 
    end
end

% Find 2235H quieter
A = ismember(Amp_singl_trials_LLR.subjID,'2235');
B = ismember(Amp_singl_trials_LLR.sessionNr,'H');
C = find(Amp_singl_trials_LLR.clickIntensity == 1);
count = 1;
for i = 1:size(Amp_singl_trials_LLR,1)
    if A(i)&& B(i) && sum(i == C)
    coincident_rows_2235H_quieter(count) = i;
    count = count +1; 
    end
end

%%%%%%%%%%%%%%%%%%%%% SECOND, MODIFY %%%%%%%%%%%%%%%%%%
% 2005D
for i = 1:length(coincident_rows_2005D_louder)
    line = coincident_rows_2005D_louder(i);
    Amp_singl_trials_LLR.clickIntensity(line) = 1;
end
for i = 1:length(coincident_rows_2005D_quieter)
    line = coincident_rows_2005D_quieter(i);
    Amp_singl_trials_LLR.clickIntensity(line) = 3;
end
% 2238G
for i = 1:length(coincident_rows_2238G_louder)
    line = coincident_rows_2238G_louder(i);
    Amp_singl_trials_LLR.clickIntensity(line) = 1;
end
for i = 1:length(coincident_rows_2238G_quieter)
    line = coincident_rows_2238G_quieter(i);
    Amp_singl_trials_LLR.clickIntensity(line) = 3;
end
% 2235H
for i = 1:length(coincident_rows_2235H_louder)
    line = coincident_rows_2235H_louder(i);
    Amp_singl_trials_LLR.clickIntensity(line) = 1;
end
for i = 1:length(coincident_rows_2235H_quieter)
    line = coincident_rows_2235H_quieter(i);
    Amp_singl_trials_LLR.clickIntensity(line) = 3;
end

if strcmp(modality_data{mode},'EEG')
    Amp_singl_trials_EEG_LLR = Amp_singl_trials_LLR;
    writetable(Amp_singl_trials_EEG_LLR, [root_dir '/Events/Single_trials/Amp_singl_trials_EEG_LLR.txt'])
    clear('Amp_singl_trials_EEG_LLR'); % to free space
elseif strcmp(modality_data{mode},'MEG')
    Amp_singl_trials_MEG_LLR = Amp_singl_trials_LLR;
    writetable(Amp_singl_trials_MEG_LLR, [root_dir '/Events/Single_trials/Amp_singl_trials_MEG_LLR.txt'])
    clear('Amp_singl_trials_MEG_LLR'); % to free space
end

end

% Reset modality_data variable to its original self
modality_data = {'EEG','MEG','both_mod'};

clearvars('-except', initialVars{:});
disp 'DONE WITH FIXING THREE SESSIONS WITH INVERTED CLICK INTENSITIES (LLR)!!!'
disp(datetime)
toc

%% Eliminate any remaining sweeps above threshold (LLR)

tic
disp(' ');      
disp('-------------------------');  
disp('ELIMINATING ANY REMAINING SWEEPS ABOVE THRESHOLD (LLR)');  
disp(datetime)
disp('-------------------------');     
disp(' ');

% Separate EEG and MEG files because they have different rejection criteria
modality_data = {'EEG'};
components_LLR = {'N1' 'P2'};
% components_MLR = {'P0' 'Na' 'Pa' 'Nb' 'P50'};

for mode = 1:length(modality_data) % EEG and MEG
    
Amp_singl_trials_LLR = readtable([root_dir '/Events/Single_trials/Amp_singl_trials_' modality_data{mode} '_LLR.txt']);

if strcmp(modality_data{mode},'EEG')
    
    A = find(Amp_singl_trials_LLR.EEG_N1 > 50);
    B = find(Amp_singl_trials_LLR.EEG_N1 < -50);
    C = find(Amp_singl_trials_LLR.EEG_P2 > 50);
    D = find(Amp_singl_trials_LLR.EEG_P2 < -50);

    % Find rows with any of these coincidences
    combined_list = [A;B;C;D];
    final_row_list = unique(combined_list);

    % Replace values in these rows by NaN
    for i = 1:length(final_row_list)
        line = final_row_list(i);
        Amp_singl_trials_LLR.EEG_N1(line) = NaN;
        Amp_singl_trials_LLR.EEG_P2(line) = NaN;
    end
    
    Amp_singl_trials_EEG_LLR = Amp_singl_trials_LLR;
    writetable(Amp_singl_trials_EEG_LLR, [root_dir '/Events/Single_trials/Amp_singl_trials_EEG_LLR.txt'])
    clear('Amp_singl_trials_EEG_LLR'); % to free space
elseif strcmp(modality_data{mode},'MEG')
    
    A = find(Amp_singl_trials_LLR.MEG_N1 > 50);
    B = find(Amp_singl_trials_LLR.MEG_N1 < -50);
    C = find(Amp_singl_trials_LLR.MEG_P2 > 50);
    D = find(Amp_singl_trials_LLR.MEG_P2 < -50);

    % Find rows with any of these coincidences
    combined_list = [A;B;C;D];
    final_row_list = unique(combined_list);

    % Replace values in these rows by NaN
    for i = 1:length(final_row_list)
        line = final_row_list(i);
        Amp_singl_trials_LLR.MEG_N1(line) = NaN;
        Amp_singl_trials_LLR.MEG_P2(line) = NaN;
    end
    
    Amp_singl_trials_MEG_LLR = Amp_singl_trials_LLR;
    writetable(Amp_singl_trials_MEG_LLR, [root_dir '/Events/Single_trials/Amp_singl_trials_MEG_LLR.txt'])
    clear('Amp_singl_trials_MEG_LLR'); % to free space
end

end

% Reset modality_data variable to its original self
modality_data = {'EEG','MEG','both_mod'};

clearvars('-except', initialVars{:});
disp 'DONE ELIMINATING ANY REMAINING SWEEPS ABOVE THRESHOLD (LLR)!!!'
disp(datetime)
toc
