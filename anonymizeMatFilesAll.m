fclose('all');
addpath('common');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mainPath = 'C:\path2experiment'; % Change this as required

info = dir([mainPath filesep '**' filesep 'CESTdata*.mat']); % Change name of mat files if necessary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nfiles = length(info);
for ii=1:nfiles
    fullname = fullfiles(info(ii).folder,info(ii).name);
    disp(['Processing ' num2str(ii) ' of ' num2str(nfiles) ': ' fullname]);
    anonymizeMatFile(fullname);
end

rmpath('common');


