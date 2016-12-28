%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TMS-EEG Single Pulse (N100 - MIN) Script %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 

% First need to create main 'N100' folder, which has 3 subfolders ('Controls',
% 'Musicians, 'MDD'), and then within each of these folders, a folder for
% each N100 stimulation paradigm (sp_lpfc, sp_lmc, sp_lpc, sp_rpfc, sp_rmc,
% sp_rpc) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Folders must have the EXACT SAME NAMES AS ABOVE %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

% Creating empty N100 structure to hold data
N100 = struct('Controls', struct('sp_lpfc', table(), 'sp_lmc', table(), 'sp_lpc', table(),...
    'sp_rpfc', table(), 'sp_rmc', table(), 'sp_rpc', table()),...
    'Musicians', struct('sp_lpfc', table(), 'sp_lmc', table(), 'sp_lpc', table(),...
    'sp_rpfc', table(), 'sp_rmc', table(), 'sp_rpc', table()),...
    'MDD', struct('sp_lpfc', table(), 'sp_lmc', table(), 'sp_lpc', table(),...
    'sp_rpfc', table(), 'sp_rmc', table(), 'sp_rpc', table())); 

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHANGE TO N100 ROOT FOLDER %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cdpath = cd;
files = dir;
directoryNames = {files([files.isdir]).name}; 
folders = directoryNames(~ismember(directoryNames,{'.','..'}));

%%

% For each group (e.g. Controls, Musicians, MDD)...
for j = 1:size(folders,2)
    
    % getting names of each subfolder 
    path_group = strcat(cdpath,'\', folders{j});
    files_groupsub = dir(path_group); 
    direcNames_groupsub = {files_groupsub([files_groupsub.isdir]).name};
    folders_groupsub = direcNames_groupsub(~ismember(direcNames_groupsub,{'.','..'}));
    
    % for each stimulation site folder (sp_lmc, sp_rmc...)...
    for stimsite = 1:size(folders_groupsub,2)
        
        % getting each .set file in folder and putting names in SP_list
        path_group_site = strcat(path_group, '\', folders_groupsub{stimsite},'\'); 
        pathin = path_group_site;
        SP_list = dir([pathin '*sp*']);
        
        % Only go through next loop (looping for each subject's site
        % specific stimulation) if there are files (a.k.a not an empty
        % folder)
        if ~isempty(SP_list) 
        
            % Loading single (first) SP to get channel information
            SP_file_chan_count = pop_loadset('filename',...
            SP_list(1).name,...
            'filepath', pathin);

            % Setting N100 epoch     
            N100_min_time = 90;
            N100_max_time = 130; 

            % Preallocating vector to hold individual N100 data (rows = subject,
            % columns = electrodes)
            N100_each_subject_min = zeros(size(SP_list,1),size(SP_file_chan_count.chanlocs(1,1:end),2));

            % for each subject at one stimulation site (e.g. N100 -> Controls -> sp_lmc --> sub1 sub2 ...)     
            for i = 1:size(SP_list,1)

                SP_file = pop_loadset('filename',...
                    SP_list(i).name,...
                    'filepath', pathin);

                SP_data = SP_file.data;

                N100_subject_min = zeros(1, numel(SP_data(:,1,1)));

                % for each electrode
                for h = 1:numel(SP_data(:,1,1))

                    SP = mean(SP_data(h,:,:), 3);
                    N100_min_time_index = find(SP_file.times == N100_min_time);
                    N100_max_time_index = find(SP_file.times == N100_max_time);
                    N100_min_value = min(SP(N100_min_time_index:N100_max_time_index));
                    N100_subject_min(h) = N100_min_value;
                    N100_each_subject_min(i,h) = N100_min_value;
                end


            % Creating dataset of N100 results for each site of stimulation
            % of each group
            N100_min_dataset = dataset({N100_each_subject_min SP_file.chanlocs(1,:).labels},...
                'obsnames',{SP_list(:,1).name});
            
            % Converting dataset to table format (recommended by MATLAB)
            N100_min_table = dataset2table(N100_min_dataset);
            
            % Resulting table is added to N100 structure 
            % folders{j} gives current group that is being looped (e.g.
            % Controls)
            % folders_groupsub(stimsite): gives for each group folder, the
            % stimulation folder being looped through current (e.g. sp_lmc)
            N100.(folders{j}).(folders_groupsub{stimsite}) = N100_min_table;
            
            % Getting specific table name (Group and site of stimulation)
            % to be used in output
            table_name = strcat(folders(j),'_',folders_groupsub{stimsite},'_', 'N100_min_table.csv'); 
            
            % Exporting dataset to .csv format to be used in excel/R
            % Will be outputted in the top directory (N100 folder)
            export(N100_min_dataset,'File',table_name{1,1},'Delimiter',',')

            end
        end
    end
end


%% 
% Uncorrected independent t-test comparison at each electrode 

t_test_results = struct('sp_lpfc',table(), 'sp_lmc', table(), 'sp_lpc', table(),...
    'sp_rpfc', table(), 'sp_rmc', table(), 'sp_rpc', table());

for stimsite = 1:size(folders_groupsub,2)
    
    N100_p_val_dataset = dataset({zeros(1,60) SP_file.chanlocs(1,:).labels});
    N100_p_val_table = dataset2table(N100_p_val_dataset);
    
    for elec = 1:size(N100.Controls.sp_rmc{1,:},2)
        
        
        % only looping through stimsites tables with data (not empty) 
        if ~isempty(N100.Musicians.(folders_groupsub{stimsite}))
        
            Musicians_elec = N100.Musicians.(folders_groupsub{stimsite})(:,elec);
            Controls_elec = N100.Controls.(folders_groupsub{stimsite})(:,elec);
        
            [H, P] = ttest2(Musicians_elec{:,:}, Controls_elec{:,:}); 
            N100_p_val_table{1,elec} = P;
        end
    end
    
    t_test_results.(folders_groupsub{stimsite}) = N100_p_val_table;
    
end



%%

% Creating new structure to make any adjustments to p-value for
% graphic/figure related reasons 
t_plot = t_test_results; 

critical_p_val = 0.05;

for stimsite = 1:size(folders_groupsub,2) 

    for i = 1:size(t_plot.(folders_groupsub{stimsite}),2)
        if t_plot.(folders_groupsub{stimsite}){1,i} > critical_p_val;
            t_plot.(folders_groupsub{stimsite}){1,i} = 0;
        else
            t_plot.(folders_groupsub{stimsite}){1,i} = -log(t_plot.(folders_groupsub{stimsite}){1,i}); 
        end
    end
end

clear criticap_p_val

%%

clear cdpath Controls_elec direcNames_groupsub directoryNames elec files files_groupsub folders...
    folders_groupsub h i j Musicians_elec N100_each_subject_min N100_max_time N100_max_time_index...
    N100_min_dataset N100_min_table N100_min_time N100_min_time_index N100_min_value N100_p_val_dataset...
    N100_p_val_table N100_subject_min P path_group path_group_site pathin SP SP_data SP_file_chan_count...
    SP_list stimsite table_name

%%
% Plot of N100 topography 

figure
subplot(2,6,1)
topoplot(mean(N100.Controls.sp_lmc{:,:}), SP_file.chanlocs); caxis([-7 7]); title('Controls:SP LMC') 

subplot(2,6,2)
topoplot(mean(N100.Controls.sp_rmc{:,:}), SP_file.chanlocs); caxis([-7 7]); title('Controls:SP RMC') 

subplot(2,6,3)
topoplot(mean(N100.Controls.sp_lpfc{:,:}), SP_file.chanlocs); caxis([-7 7]) ; title('Controls:SP LPFC')

subplot(2,6,4)
topoplot(mean(N100.Controls.sp_rpfc{:,:}), SP_file.chanlocs); caxis([-7 7]); title('Controls:SP RPFC')

subplot(2,6,5)
topoplot(mean(N100.Controls.sp_lpc{:,:}), SP_file.chanlocs); caxis([-7 7]) ; title('Controls:SP LPC')

subplot(2,6,6)
topoplot(mean(N100.Controls.sp_rpc{:,:}), SP_file.chanlocs); caxis([-7 7]); title('Controls:SP RPC')


subplot(2,6,7)
topoplot(mean(N100.Musicians.sp_lmc{:,:}), SP_file.chanlocs); caxis([-7 7]); title('Musicians:SP LMC') 

subplot(2,6,8)
topoplot(mean(N100.Musicians.sp_rmc{:,:}), SP_file.chanlocs); caxis([-7 7]); title('Musicians:SP RMC') 

subplot(2,6,9)
topoplot(mean(N100.Musicians.sp_lpfc{:,:}), SP_file.chanlocs); caxis([-7 7]) ; title('Musicians:SP LPFC') 

subplot(2,6,10)
topoplot(mean(N100.Musicians.sp_rpfc{:,:}), SP_file.chanlocs); caxis([-7 7]) ; title('Musicians:SP RPFC')

subplot(2,6,11)
topoplot(mean(N100.Musicians.sp_lpc{:,:}), SP_file.chanlocs); caxis([-7 7]) ; title('Musicians:SP LPC') 

subplot(2,6,12)
topoplot(mean(N100.Musicians.sp_rpc{:,:}), SP_file.chanlocs); caxis([-7 7]) ; title('Musicians:SP RPC')

%%

% Plot of UNCORRECTED independent t-test at each electrode

figure
subplot(2,3,1)
topoplot(-log(t_test_results.sp_lmc{:,:}), SP_file.chanlocs); caxis([0 3]); title('LMC') 

subplot(2,3,4)
topoplot(-log(t_test_results.sp_rmc{:,:}), SP_file.chanlocs); caxis([0 3]); title('RMC') 

subplot(2,3,2)
topoplot(-log(t_test_results.sp_lpfc{:,:}), SP_file.chanlocs); caxis([0 3]) ; title('LPFC')

subplot(2,3,5)
topoplot(-log(t_test_results.sp_rpfc{:,:}), SP_file.chanlocs); caxis([0 3]); title('RPFC')

subplot(2,3,3)
topoplot(-log(t_test_results.sp_lpc{:,:}), SP_file.chanlocs); caxis([0 3]) ; title('LPC')

subplot(2,3,6)
topoplot(-log(t_test_results.sp_rpc{:,:}), SP_file.chanlocs); caxis([0 3]); title('RPC')


%%

% Plotting of p-values thresholded (i.e. only significant at specific
% p-value)
figure
subplot(2,3,1)
topoplot(t_plot.sp_lmc{:,:}, SP_file.chanlocs); caxis([0 3]); title('LMC') 

subplot(2,3,4)
topoplot(t_plot.sp_rmc{:,:}, SP_file.chanlocs); caxis([0 3]); title('RMC') 

subplot(2,3,2)
topoplot(t_plot.sp_lpfc{:,:}, SP_file.chanlocs); caxis([0 3]) ; title('LPFC')

subplot(2,3,5)
topoplot(t_plot.sp_rpfc{:,:}, SP_file.chanlocs); caxis([0 3]); title('RPFC')

subplot(2,3,3)
topoplot(t_plot.sp_lpc{:,:}, SP_file.chanlocs); caxis([0 3]) ; title('LPC')

subplot(2,3,6)
topoplot(t_plot.sp_rpc{:,:}, SP_file.chanlocs); caxis([0 3]); title('RPC')
