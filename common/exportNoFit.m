function exportNoFit(fid,app)

nrois = app.fit.nrois;
parnames = [{'ROI'},{'Pixels Post'},{'Pixels Pre'},{'Increase'},{'Recov 50'},{'Baseline'}];
nthresh = app.fit.nthr;
if isfield(app.params,'roithresh')
    nthr_abs = length(app.params.roithresh);
else
    nthr_abs = 0;
end

for kk=1:nthresh
    if nthresh>1
        if kk>nthr_abs
            fprintf(fid,'%s,%f\n','THRESHOLD',app.params.roithreshpercent(kk-nthr_abs));
        else
            fprintf(fid,'%s,%f\n','THRESHOLD',app.params.roithresh(kk));
        end
    end
    
    fprintf(fid,[repmat('%s,',1,6) '\n'],parnames{:});
    str = {'time','mean','std'};
    if isfield(app.params,'roinames') && ~isempty(app.params.roinames)
        lab = {'',app.params.roinames{1},''};
    else
        lab = {'','ROI 1',''};
    end
    if isfield(app.params,'roinames') && ~isempty(app.params.roinames)
        N = cat(1,app.fit.roipix.post(kk,:),app.fit.roipix.pre,app.fit.simplepars.increase(:,kk)',app.fit.simplepars.Tt50.spline(:,kk)',app.fit.simplepars.baseline(:,kk)');
    else
        N = cat(1,1:nrois,app.fit.roipix.post(kk,:),app.fit.roipix.pre,app.fit.simplepars.increase(:,kk)',app.fit.simplepars.Tt50.spline(:,kk)',app.fit.simplepars.baseline(:,kk)');
    end
    if nrois>1
        for ii=2:nrois
            str = [str {'mean','std'}];
            if isfield(app.params,'roinames') && ~isempty(app.params.roinames)
                lab = [lab {app.params.roinames{ii},''}];
            else
                lab = [lab {['ROI ' num2str(ii)],''}];
            end
        end
    end
    if isfield(app.params,'roinames') && ~isempty(app.params.roinames)
        for ii=1:nrois
            fprintf(fid,['%s,' repmat('%3.2f,',1,5) '\n'],app.params.roinames{ii},N(:,ii));
        end
    else
        fprintf(fid,[repmat('%3.2f,',1,6) '\n'],N);
    end
    fprintf(fid,'\n\n');
    ncols = 2*nrois + 1;
    fprintf(fid,[repmat('%s,',1,ncols) '\n'],lab{:});
    fprintf(fid,[repmat('%s,',1,ncols) '\n'],str{:});
    M = zeros(length(app.fit.timepoints),ncols);
    M(:,1) = app.fit.timepoints(:);
    M(:,2:2:end) = transpose(app.fit.means(:,:,kk));
    M(:,3:2:end) = transpose(app.fit.stds(:,:,kk));
    
    fprintf(fid,[repmat('%3.3f,',1,ncols) '\n'],M.');
    fprintf(fid,'\n\n');
end

fclose(fid);

return;