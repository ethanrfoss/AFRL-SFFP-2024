function FunctionStruct = FunctionInputs(f)

% Get Function Input Names:
FunctionInfo = functions(f);
InputNames = strtrim(strsplit(FunctionInfo.function(strfind(FunctionInfo.function,'(')+1:strfind(FunctionInfo.function,')')-1),','));

% Get Input Dimensions:
dims = InputDimensions(f);

% Make Struct:
for i = 1:length(dims)
    FunctionStruct.(InputNames{i}) = dims(i);
end

end