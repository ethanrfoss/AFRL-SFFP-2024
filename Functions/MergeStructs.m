
function [StructB] = MergeStructs(StructA,StructB)

f = fieldnames(StructA);
for i = 1:length(f)
    % if ~isfield(StructB,f{i})
    %     warning('Adding Unknown Field in Structure Merging');
    % end
    StructB.(f{i}) = StructA.(f{i});
end

end