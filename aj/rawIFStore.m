function ok = rawIFStore(action, id, el)
%RAWIFSTORE Keep / fetch the 8-element raw IF of a campaign item.
%   rawIFStore('keep',  id, el)  moves the element folder el (8 element files,
%                                their _time.mat, the truth .mat) to
%                                J:\AJ-AS\code\campaign\rawIF\<id> and removes el
%   ok = rawIFStore('fetch', id, el)  copies a stored raw IF back to el (fast disk)
%                                for re-processing; false when none is stored
%   ok = rawIFStore('has', id)   true when the raw IF of id is stored
%   The raw IF is deterministic (generateIFdataArray with seed 100 + seed), so a
%   stored set equals a regenerated one; storing it saves the ~20 min generation
%   when array processing or the receiver is changed later.
store = fullfile(ajPath('campaign', 'rawIF'), id);
tok = regexp(id, '^(S\d)_(.+)_s(\d+)$', 'tokens', 'once');
truthName = sprintf('navic8_truth_Fighter_%s_%s.mat', tok{2}, id);
ok = false;
switch action
    case 'has'
        ok = exist(fullfile(store, truthName), 'file') == 2 && numel(dir(fullfile(store, 'navic8_e*_*.bin'))) == 8;
    case 'keep'
        if ~exist(el, 'dir'), return; end
        if ~exist(store, 'dir'), mkdir(store); end
        F = dir(fullfile(el, '*'));  F = F(~[F.isdir]);
        for f = F'
            src = fullfile(f.folder, f.name);  dst = fullfile(store, f.name);
            copyfile(src, dst);
            d = dir(dst);
            if isempty(d) || d.bytes ~= f.bytes, error('rawIFStore: copy of %s incomplete', f.name); end
        end
        rmdir(el, 's');
        ok = true;
    case 'fetch'
        if ~rawIFStore('has', id), return; end
        if ~exist(el, 'dir'), mkdir(el); end
        copyfile(fullfile(store, '*'), el);
        ok = true;
    otherwise
        error('rawIFStore: unknown action %s', action);
end
end
