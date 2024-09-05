
function p = plot(Z,varargin)

if isempty(varargin)
    p = plot(Polygon(Z));
else
    p = plot(Polygon(Z),varargin{:});
end

end