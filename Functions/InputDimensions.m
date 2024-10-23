
function dims = InputDimensions(f,MaxDims)

if nargin < 2
    MaxDims = 128;
end

NumInputs = nargin(f);
Inputs = cell(NumInputs,1);
dims = zeros(NumInputs,1);
for i = 1:NumInputs

    [Inputs{i+1:end}] = deal(zeros(MaxDims,1));
    j = 1;
    while j < MaxDims
        Inputs{i} = zeros(j,1);
        try
            f(Inputs{:});
            break;
        catch
            j = j + 1;
        end
    end
    if j == MaxDims
        error('Max Dimension Count Reached');
    end
    dims(i) = j;
end

end