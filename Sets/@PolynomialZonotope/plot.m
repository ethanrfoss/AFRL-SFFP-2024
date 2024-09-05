
function p = plot(PZ,varargin)

if isempty(varargin)
    varargin = {'FaceColor',[0 0.4470 0.7410],'FaceAlpha',.5,'EdgeColor',[0 0 0]};
    Options.Splits = 10;
elseif isstruct(varargin{1})
    Options = varargin{1};
    varargin(1) = [];
else
    Options.Splits = 10;
end

if PZ.dims.n == 2
    p = Plot2D(PZ,Options,varargin);
elseif PZ.dims.n == 3
    p = Plot3D(PZ,Options,varargin);
else
    error('Invalid Plotting Dimension');
end

end

function p = Plot2D(PZ,Options,varargin)

PZSplit = SplitMax(PZ);
% PolySplit = Polygon(PZSplit);
PolySplit = arrayfun(@(obj) obj.Polygon,PZSplit);
PolyFull = SetOperation(PolySplit,'union');

for i = 1:Options.Splits
    
    PZSplitNew = [];
    PolySplitNew = [];
    for j = 1:length(PZSplit)
    
        if any(sum(ismembertol(PolySplit(j).Vertices,PolyFull.Vertices,1e-10),2) == 2)
            PZSplitNew = [PZSplitNew SplitMax(PZSplit(j))];
            PolySplitNew = [PolySplitNew arrayfun(@(obj) obj.Polygon,PZSplitNew(end-1:end))];
        else
            PZSplitNew = [PZSplitNew PZSplit(j)];
            PolySplitNew = [PolySplitNew PolySplit(j)];
        end

    end
    PolyFull = SetOperation(PolySplitNew,'union');
    PZSplit = PZSplitNew;
    PolySplit = PolySplitNew;

    % hold on;
    % for k = 1:length(PolySplit)
    %     plot(PolySplit(k));
    % end

end

p = plot(PolyFull);

% hold on;
% for i = 1:length(PolySplit)
%     plot(PolySplit(i));
% end

end