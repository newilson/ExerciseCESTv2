fclose('all');
clear,clc
addpath('common');

filtstrength = 'weak'; % 'strong', 'moderate','weak','none'

% Change all params
%   true to pick file
%   FileName to read specific csv
changeParams = true;

% Change specific param(s)
changeIndividualParam = false;
% changeIndividualParam = {'b0tol','b1tol'};
% newIndividualParam = {0.3, 0.6};


% mainpath = 'C:\Users\CMROI\Desktop\CrCEST-reproducibility_forNeil';
% mainpath = 'D:\ZAMANI';
mainpath = 'D:\SHEIKH\823664-Sheikh';

% For Zamani
% info = dir([mainpath filesep '**' filesep 'CESTdata*2017.mat']);
% info = cat(1,info,dir([mainpath filesep '**' filesep 'CESTdata*2018.mat']));
% info = dir([mainpath filesep '**' filesep 'CESTdata*_bd_mf.mat']);
% info = dir([mainpath filesep '**' filesep 'CESTdata22-Aug-2019.mat']);
info = [dir([mainpath filesep '**' filesep 'CESTdata*_bd_mf.mat']); dir([mainpath filesep '**' filesep 'CESTdata*2021.mat'])];

% info = cat(1,info,dir([mainpath filesep '**' filesep 'CESTdata*Jan-2018_GRE.mat']));
% info = cat(1,info,dir([mainpath filesep '**' filesep 'CESTdata*Jan-2018_WASSR.mat']));
% info = cat(1,info,dir([mainpath filesep '**' filesep 'CESTdata*May-2017.mat']));
% info = cat(1,info,dir([mainpath filesep '**' filesep 'CESTdata*Jan-2018_WASSR.mat']));

% For Sheikh
% info = dir([mainpath filesep '**' filesep 'CESTdata*Apr-2018.mat']);
% info = cat(1,info,dir([mainpath filesep '**' filesep 'CESTdata*Apr-2019.mat']));


tic
nfiles = length(info)
for ii=1:nfiles
    fullname = fullfile(info(ii).folder,info(ii).name);
    [fpath,fname,fext] = fileparts(fullname);
%     newname = [fname '_bd_mf' fext];
%     newname = [fname '_new' fext];
    newname = ['CESTdata' date '.mat'];
    if ~exist(fullfile(fpath,newname),'file') && isempty(strfind(fpath,'DO_NOT_USE')) % && isempty(strfind(fpath,'ZAMANI_PAYMAN'))
        disp(['processing ' fullname]);
        clear out
        load(fullname)

        if isfolder(changeParams) || changeParams
            if isfolder(changeParams)
                out = setParams(out,changeParams);
            else
                out = setParams(out);
            end
        end

        if iscell(changeIndividualParam) || isstring(changeIndividualParam)
            if iscell(changeIndividualParam)
                for jj=1:length(changeIndividualParam)
                    out = modifyParam(out,changeIndividualParam{jj},newIndividualParam{jj});
                end
            else
                out = modifyParam(out,changeIndividualParam,newIndividualParam);
            end
        end
        
        out.mainDir = info(ii).folder;
        try
            out = reprocessExerciseCESTwithROIs(out,filtstrength);
            save(fullfile(fpath,newname),'out');
            clc
            disp(['saved #' num2str(ii) ' ' newname]);
        catch
            warndlg([fullfile(fpath,newname) ' did not save.'],'WARNING','nonmodal');
            disp(['ERROR with ' fullfile(fpath,newname)]);
        end
    end
end
toc