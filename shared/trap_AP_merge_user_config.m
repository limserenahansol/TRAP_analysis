function C = trap_AP_merge_user_config(userC)
    C = trap_config();
    if nargin < 1 || isempty(userC) || ~isstruct(userC)
        return;
    end
    fn = fieldnames(userC);
    for k = 1:numel(fn)
        C.(fn{k}) = userC.(fn{k});
    end
end
