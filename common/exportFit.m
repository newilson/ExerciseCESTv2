function exportFit(fid,app)

nrois = app.fit.nrois;
npars = length(app.params.parnames);
nweights = app.fit.nweights;
parnames = [{'ROI'},repmat(app.params.parnames,[1 nweights]),{'NoFitInc'},{'NoFitRecov_sp'},{'NoFitRecov_lin'},{'NoFitBase'}];
nthresh = size(app.fit.pars,4);
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
    
    fprintf(fid,[repmat('%s,',1,npars*nweights+1+4) '\n'],parnames{:});
    str = {'time','mean','std','fit'};
    if isfield(app.params,'roinames') && ~isempty(app.params.roinames)
        lab = {'',app.params.roinames{1},'',''};
    else
        lab = {'','ROI 1','',''};
    end
    if nrois>1
        if isfield(app.params,'roinames') && ~isempty(app.params.roinames)
            N = cat(1,reshape(permute(app.fit.pars(:,:,:,kk),[1 3 2]),[],nrois),app.fit.roipix.post(kk,:),app.fit.roipix.pre,app.fit.simplepars.increase(:,kk)',app.fit.simplepars.Tt50.spline(:,kk)',app.fit.simplepars.Tt50.linear(:,kk)',app.fit.simplepars.baseline(:,kk)');
            Nlb95 = reshape(permute(squeeze(app.fit.confidence_intervals95(:,1,:,:,kk)),[1 3 2]),[],nrois);
            Nub95 = reshape(permute(squeeze(app.fit.confidence_intervals95(:,2,:,:,kk)),[1 3 2]),[],nrois);
            Nlb90 = reshape(permute(squeeze(app.fit.confidence_intervals90(:,1,:,:,kk)),[1 3 2]),[],nrois);
            Nub90 = reshape(permute(squeeze(app.fit.confidence_intervals90(:,2,:,:,kk)),[1 3 2]),[],nrois);
        else
            N = cat(1,1:nrois,reshape(permute(app.fit.pars(:,:,:,kk),[1 3 2]),[],nrois),app.fit.roipix.post(kk,:),app.fit.roipix.pre,app.fit.simplepars.increase(:,kk)',app.fit.simplepars.Tt50.spline(:,kk)',app.fit.simplepars.Tt50.linear(:,kk)',app.fit.simplepars.baseline(:,kk)');
            Nlb95 = cat(1,1:nrois,reshape(permute(squeeze(app.fit.confidence_intervals95(:,1,:,:,kk)),[1 3 2]),[],nrois));
            Nub95 = cat(1,1:nrois,reshape(permute(squeeze(app.fit.confidence_intervals95(:,2,:,:,kk)),[1 3 2]),[],nrois));
            Nlb90 = cat(1,1:nrois,reshape(permute(squeeze(app.fit.confidence_intervals90(:,1,:,:,kk)),[1 3 2]),[],nrois));
            Nub90 = cat(1,1:nrois,reshape(permute(squeeze(app.fit.confidence_intervals90(:,2,:,:,kk)),[1 3 2]),[],nrois));
        end
        for ii=2:nrois
            str = [str {'mean','std'}];
            for ll=1:nweights
                str = [str {'fit'}];
            end
            if isfield(app.params,'roinames') && ~isempty(app.params.roinames)
                lab = [lab {app.params.roinames{ii}}];
            else
                lab = [lab {['ROI ' num2str(ii)]}];
            end
            for ll=1:1+nweights
                lab = [lab {''}];
            end
        end
    end
    if isfield(app.params,'globalroi') && app.params.globalroi
        lab = [lab {'Global'}];
        for ll=1:1+nweights
            lab = [lab {''}];
        end
        str = [str {'mean','std'}];
        for ll=1:nweights
            str = [str {'fit'}];
        end
    end
    if isfield(app.params,'roinames') && ~isempty(app.params.roinames)
        for ii=1:nrois
            fprintf(fid,['%s,' repmat('%3.2f,',1,npars+4) '\n'],app.params.roinames{ii},N(:,ii));
        end
        if isfield(app.params,'globalroi') && app.params.globalroi
            fprintf(fid,['%s,' repmat('%3.2f,',1,npars+4) '\n'],'Global',app.fit.glo.pars(:,kk),app.fit.roipix.postglo(kk),app.fit.roipix.preglo,app.fit.glo.simplepars.increase(kk),app.fit.glo.simplepars.Tt50.spline(kk),app.fit.glo.simplepars.Tt50.linear(kk),app.fit.glo.simplepars.baseline(kk));
        end
        fprintf(fid,'\n Confidence Intervals');
        fprintf(fid,'\n 95%%\n');
        for ii=1:nrois
            fprintf(fid,['%s,' repmat('%3.2f,',1,npars-2) '\n'],app.params.roinames{ii},Nlb95(:,ii));
        end
        if isfield(app.params,'globalroi') && app.params.globalroi
            fprintf(fid,['%s,' repmat('%3.2f,',1,npars-2) '\n'],'Global',squeeze(app.fit.glo.confidence_intervals95(:,1,kk)));
        end
        fprintf(fid,'\n');
        for ii=1:nrois
            fprintf(fid,['%s,' repmat('%3.2f,',1,npars-2) '\n'],app.params.roinames{ii},Nub95(:,ii));
        end
        if isfield(app.params,'globalroi') && app.params.globalroi
            fprintf(fid,['%s,' repmat('%3.2f,',1,npars-2) '\n'],'Global',squeeze(app.fit.glo.confidence_intervals95(:,2,kk)));
        end
        fprintf(fid,'\n 90%%\n');
        for ii=1:nrois
            fprintf(fid,['%s,' repmat('%3.2f,',1,npars-2) '\n'],app.params.roinames{ii},Nlb90(:,ii));
        end
        if isfield(app.params,'globalroi') && app.params.globalroi
            fprintf(fid,['%s,' repmat('%3.2f,',1,npars-2) '\n'],'Global',squeeze(app.fit.glo.confidence_intervals90(:,1,kk)));
        end
        fprintf(fid,'\n');
        for ii=1:nrois
            fprintf(fid,['%s,' repmat('%3.2f,',1,npars-2) '\n'],app.params.roinames{ii},Nub90(:,ii));
        end
        if isfield(app.params,'globalroi') && app.params.globalroi
            fprintf(fid,['%s,' repmat('%3.2f,',1,npars-2) '\n'],'Global',squeeze(app.fit.glo.confidence_intervals90(:,2,kk)));
        end
    else
        fprintf(fid,[repmat('%3.2f,',1,npars+1) '\n'],N);
        if isfield(app.params,'globalroi') && app.params.globalroi
            fprintf(fid,['%s,' repmat('%3.2f,',1,npars) '\n'],'Global',app.fit.glo.pars(:,:,kk),app.fit.glo.roipix.post(kk,:),app.fit.glo.roipix.pre);
        end
        fprintf(fid,'\n Confidence Intervals');
        fprintf(fid,'\n 95%%\n');
        fprintf(fid,[repmat('%3.2f,',1,npars-1) '\n'],Nlb95);
        if isfield(app.params,'globalroi') && app.params.globalroi
            fprintf(fid,['%s,' repmat('%3.2f,',1,npars-2) '\n'],'Global',squeeze(app.fit.glo.confidence_intervals95(:,1,kk)));
        end
        fprintf(fid,'\n');
        fprintf(fid,[repmat('%3.2f,',1,npars-1) '\n'],Nub95);
        if isfield(app.params,'globalroi') && app.params.globalroi
            fprintf(fid,['%s,' repmat('%3.2f,',1,npars-2) '\n'],'Global',squeeze(app.fit.glo.confidence_intervals95(:,2,kk)));
        end
        fprintf(fid,'\n 90%%\n');
        fprintf(fid,[repmat('%3.2f,',1,npars-1) '\n'],Nlb90);
        if isfield(app.params,'globalroi') && app.params.globalroi
            fprintf(fid,['%s,' repmat('%3.2f,',1,npars-2) '\n'],'Global',squeeze(app.fit.glo.confidence_intervals90(:,1,kk)));
        end
        fprintf(fid,'\n');
        fprintf(fid,[repmat('%3.2f,',1,npars-1) '\n'],Nub90);
        if isfield(app.params,'globalroi') && app.params.globalroi
            fprintf(fid,['%s,' repmat('%3.2f,',1,npars-2) '\n'],'Global',squeeze(app.fit.glo.confidence_intervals90(:,2,kk)));
        end
    end
    fprintf(fid,'\n\n');
    ntemp = nweights + 2;
    if isfield(app.params,'globalroi') && app.params.globalroi
        ncols = ntemp*(nrois+1) + 1;
    else
        ncols = ntemp*nrois + 1;
    end
    fprintf(fid,[repmat('%s,',1,ncols) '\n'],lab{:});
    fprintf(fid,[repmat('%s,',1,ncols) '\n'],str{:});
    M = zeros(length(app.fit.timepoints),ncols);
    M(:,1) = app.fit.timepoints(:);
    if isfield(app.params,'globalroi') && app.params.globalroi
        M(:,2:ntemp:end) = transpose(cat(1,app.fit.means(:,:,kk),app.fit.glo.means(:,kk).'));
        M(:,3:ntemp:end) = transpose(cat(1,app.fit.stds(:,:,kk),app.fit.glo.stds(:,kk).'));
        M(:,4:ntemp:end) = cat(2,app.fit.fitdata(:,:,kk),app.fit.glo.fitdata(:,kk));
    else
        M(:,2:ntemp:end) = transpose(app.fit.means(:,:,kk));
        M(:,3:ntemp:end) = transpose(app.fit.stds(:,:,kk));
        M(:,4:ntemp:end) = app.fit.fitdata(:,:,kk);
    end
    fprintf(fid,[repmat('%3.3f,',1,ncols) '\n'],M.');
    fprintf(fid,'\n\n');
end

if isfield(app.params,'myfun_2') && isfield(app.params,'par0_2') && isfield(app.params,'parnames_2')
    if strcmp(app.params.parnames_2{end-3},'Voxels Post')
        app.params.parnames_2 = app.params.parnames_2{1:end-2};
    end
    npars = length(app.params.parnames_2);
    parnames = [{'ROI'},app.params.parnames_2];
    nthresh = size(app.fit.pars_2,3);
    
    for kk=1:nthresh
        if nthresh>1
            if kk>nthr_abs
                fprintf(fid,'%s,%f\n','THRESHOLD',app.params.roithreshpercent(kk-nthr_abs));
            else
                fprintf(fid,'%s,%f\n','THRESHOLD',app.params.roithresh(kk));
            end
        end
        
        fprintf(fid,[repmat('%s,',1,npars+1) '\n'],parnames{:});
        str = {'time','mean','std','fit'};
        if isfield(app.params,'roinames')
            lab = {'',app.params.roinames{1},'',''};
        else
            lab = {'','ROI 1','',''};
        end
        if nrois>1
            if isfield(app.params,'roinames') && ~isempty(app.params.roinames)
                N = cat(1,app.fit.pars_2(:,:,kk),app.fit.roipix.post(kk,:),app.fit.roipix.pre);
                Nlb95 = squeeze(app.fit.confidence_intervals95_2(:,1,:,kk));
                Nub95 = squeeze(app.fit.confidence_intervals95_2(:,2,:,kk));
                Nlb90 = squeeze(app.fit.confidence_intervals90_2(:,1,:,kk));
                Nub90 = squeeze(app.fit.confidence_intervals90_2(:,2,:,kk));
            else
                N = cat(1,1:nrois,app.fit.pars_2(:,:,kk),app.fit.roipix.post(kk,:),app.fit.roipix.pre);
                Nlb95 = cat(1,1:nrois,squeeze(app.fit.confidence_intervals95_2(:,1,:,kk)));
                Nub95 = cat(1,1:nrois,squeeze(app.fit.confidence_intervals95_2(:,2,:,kk)));
                Nlb90 = cat(1,1:nrois,squeeze (app.fit.confidence_intervals90_2(:,1,:,kk)));
                Nub90 = cat(1,1:nrois,squeeze(app.fit.confidence_intervals90_2(:,2,:,kk)));
            end
            for ii=2:nrois
                str = [str {'mean','std','fit'}];
                if isfield(app.params,'roinames')
                    lab = [lab {app.params.roinames{ii},'',''}];
                else
                    lab = [lab {['ROI ' num2str(ii)],'',''}];
                end
            end
        end
        if isfield(app.params,'globalroi') && app.params.globalroi
            lab = [lab {'Global','',''}];
            str = [str {'mean','std','fit'}];
        end
        if isfield(app.params,'roinames') && ~isempty(app.params.roinames)
            for ii=1:nrois
                fprintf(fid,['%s,' repmat('%3.2f,',1,npars) '\n'],app.params.roinames{ii},N(:,ii));
            end
            if isfield(app.params,'globalroi') && app.params.globalroi
                fprintf(fid,['%s,' repmat('%3.2f,',1,npars) '\n'],'Global',app.fit.glo.pars_2(:,kk),app.fit.roipix.postglo(kk),app.fit.roipix.preglo);
            end
            fprintf(fid,'\n Confidence Intervals');
            fprintf(fid,'\n 95%%\n');
            for ii=1:nrois
                fprintf(fid,['%s,' repmat('%3.2f,',1,npars-2) '\n'],app.params.roinames{ii},Nlb95(:,ii));
            end
            if isfield(app.params,'globalroi') && app.params.globalroi
                fprintf(fid,['%s,' repmat('%3.2f,',1,npars-2) '\n'],'Global',squeeze(app.fit.glo.confidence_intervals95_2(:,1,kk)));
            end
            fprintf(fid,'\n');
            for ii=1:nrois
                fprintf(fid,['%s,' repmat('%3.2f,',1,npars-2) '\n'],app.params.roinames{ii},Nub95(:,ii));
            end
            if isfield(app.params,'globalroi') && app.params.globalroi
                fprintf(fid,['%s,' repmat('%3.2f,',1,npars-2) '\n'],'Global',squeeze(app.fit.glo.confidence_intervals95_2(:,2,kk)));
            end
            fprintf(fid,'\n 90%%\n');
            for ii=1:nrois
                fprintf(fid,['%s,' repmat('%3.2f,',1,npars-2) '\n'],app.params.roinames{ii},Nlb90(:,ii));
            end
            if isfield(app.params,'globalroi') && app.params.globalroi
                fprintf(fid,['%s,' repmat('%3.2f,',1,npars-2) '\n'],'Global',squeeze(app.fit.glo.confidence_intervals90_2(:,1,kk)));
            end
            fprintf(fid,'\n');
            for ii=1:nrois
                fprintf(fid,['%s,' repmat('%3.2f,',1,npars-2) '\n'],app.params.roinames{ii},Nub90(:,ii));
            end
            if isfield(app.params,'globalroi') && app.params.globalroi
                fprintf(fid,['%s,' repmat('%3.2f,',1,npars-2) '\n'],'Global',squeeze(app.fit.glo.confidence_intervals90_2(:,2,kk)));
            end
        else
            fprintf(fid,[repmat('%3.2f,',1,npars+1) '\n'],N);
            if isfield(app.params,'globalroi') && app.params.globalroi
                fprintf(fid,['%s,' repmat('%3.2f,',1,npars) '\n'],'Global',app.fit.glo.pars_2(:,:,kk),app.fit.glo.roipix.post(kk,:),app.fit.glo.roipix.pre);
            end
            fprintf(fid,'\n Confidence Intervals');
            fprintf(fid,'\n 95%%\n');
            fprintf(fid,[repmat('%3.2f,',1,npars-1) '\n'],Nlb95);
            if isfield(app.params,'globalroi') && app.params.globalroi
                fprintf(fid,['%s,' repmat('%3.2f,',1,npars-2) '\n'],'Global',squeeze(app.fit.glo.confidence_intervals95_2(:,1,kk)));
            end
            fprintf(fid,'\n');
            fprintf(fid,[repmat('%3.2f,',1,npars-1) '\n'],Nub95);
            if isfield(app.params,'globalroi') && app.params.globalroi
                fprintf(fid,['%s,' repmat('%3.2f,',1,npars-2) '\n'],'Global',squeeze(app.fit.glo.confidence_intervals95_2(:,2,kk)));
            end
            fprintf(fid,'\n 90%%\n');
            fprintf(fid,[repmat('%3.2f,',1,npars-1) '\n'],Nlb90);
            if isfield(app.params,'globalroi') && app.params.globalroi
                fprintf(fid,['%s,' repmat('%3.2f,',1,npars-2) '\n'],'Global',squeeze(app.fit.glo.confidence_intervals90_2(:,1,kk)));
            end
            fprintf(fid,'\n');
            fprintf(fid,[repmat('%3.2f,',1,npars-1) '\n'],Nub90);
            if isfield(app.params,'globalroi') && app.params.globalroi
                fprintf(fid,['%s,' repmat('%3.2f,',1,npars-2) '\n'],'Global',squeeze(app.fit.glo.confidence_intervals90_2(:,2,kk)));
            end
        end
        fprintf(fid,'\n\n');
        if isfield(app.params,'globalroi') && app.params.globalroi
            ncols = 3*(nrois+1) + 1;
        else
            ncols = 3*nrois + 1;
        end
        fprintf(fid,[repmat('%s,',1,ncols) '\n'],lab{:});
        fprintf(fid,[repmat('%s,',1,ncols) '\n'],str{:});
        M = zeros(length(app.fit.timepoints),ncols);
        M(:,1) = app.fit.timepoints(:);
        if isfield(app.params,'globalroi') && app.params.globalroi
            M(:,2:3:end) = transpose(cat(1,app.fit.means(:,:,kk),app.fit.glo.means(:,kk).'));
            M(:,3:3:end) = transpose(cat(1,app.fit.stds(:,:,kk),app.fit.glo.stds(:,kk).'));
            M(:,4:3:end) = cat(2,app.fit.fitdata_2(:,:,kk),app.fit.glo.fitdata_2(:,kk));
        else
            M(:,2:3:end) = transpose(app.fit.means(:,:,kk));
            M(:,3:3:end) = transpose(app.fit.stds(:,:,kk));
            M(:,4:3:end) = app.fit.fitdata_2(:,:,kk);
        end
        fprintf(fid,[repmat('%3.3f,',1,ncols) '\n'],M.');
        fprintf(fid,'\n\n');
    end
end
    
fclose(fid);
    
return;