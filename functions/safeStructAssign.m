
function def_struct = safeStructAssign(def_struct, new_struct, exclude_fields, add_new_fields)
% 
% function targ_struct = safeStructAssign(targ_struct, ref_struct, exclude_fields, add_new_fields)
% 
% 

if ~exist('add_new_fields', 'var')
    add_new_fields = false;
end

if ~exist('exclude_fields', 'var') || isempty(exclude_fields)
    exclude_fields = {};
elseif ~iscell(exclude_fields)
    exclude_fields = {exclude_fields};
end

if ~isstruct(def_struct)
    return;
end

if iscell(new_struct) && isempty(new_struct)
    return;
elseif iscell(new_struct) && isstruct(new_struct{1})
    new_struct = new_struct{1};
end

    %override existing values of targ_struct with those in ref_struct
if isstruct(new_struct)
    names = fieldnames(new_struct);
    for i = 1 : length(names)
        if ~ismember(names{i},exclude_fields)
            if isfield(def_struct,names{i}) || add_new_fields
                def_struct.(names{i}) = new_struct.(names{i});
            end
        end
    end
elseif iscell(new_struct)
    for n=1:length(new_struct)/2
        if ~ismember(new_struct{n*2-1}, exclude_fields)
            if isfield(def_struct, new_struct{n*2-1}) || add_new_fields
                def_struct.(new_struct{n*2-1}) = new_struct{n*2};
            end
        end
    end
end





