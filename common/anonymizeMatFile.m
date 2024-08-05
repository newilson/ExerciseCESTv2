function anonymizeMatFile(matfile)
%
%
% only refdicomhdrs need to be anonymized
%

load(matfile);

% anonymize
out.pre.refdicomhdr = anonymizeDicomHdr(out.pre.refdicomhdr);
out.post.refdicomhdr = anonymizeDicomHdr(out.post.refdicomhdr);

save(matfile,'out')

