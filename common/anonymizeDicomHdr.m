function anonymizedHdr = anonymizeDicomHdr(fullHdr)
%
%
% hdr is a struct
%

anonymizedHdr = fullHdr;

fields = fieldnames(anonymizedHdr);

% Remove Patient info
inds = find(contains(fields,'Patient'));
for ii=1:length(inds)
    if ~(strcmp(fields{inds(ii)},'PatientPosition') || strcmp(fields{inds(ii)},'ImagePositionPatient') || strcmp(fields{inds(ii)},'ImageOrientationPatient') ) % exceptions
        anonymizedHdr = rmfield(anonymizedHdr,fields{inds(ii)});
    end
end

% Remove Instituion info
inds = find(contains(fields,'Institution'));
for ii=1:length(inds)
    if true % no exceptions yet
        anonymizedHdr = rmfield(anonymizedHdr,fields{inds(ii)});
    end
end

