
function p = plot(Poly2D,varargin)

if isempty(varargin)
    varargin = {[0 0.4470 0.7410],'FaceAlpha',.5,'EdgeColor',[0 0 0]};
end

p = fill(Poly2D.Vertices(:,1),Poly2D.Vertices(:,2),varargin{:});

end