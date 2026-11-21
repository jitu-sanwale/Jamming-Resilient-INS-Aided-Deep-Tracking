function L = decimateLog(L)
% array log stored in the result file at 10 ms (the 1-ms log stays with the run)
if isfield(L, 't') && numel(L.t) > 1 && median(diff(L.t)) < 5e-3
    k = 10:10:numel(L.t);
    for f = {'t', 'tbMs', 'gainSat', 'psi', 'sinrGainDb', 'blankFrac', 'gainTrueDb', 'psiTrue', 'sinrGainTrueDb'}
        if isfield(L, f{1}), L.(f{1}) = L.(f{1})(k, :); end
    end
end
end
