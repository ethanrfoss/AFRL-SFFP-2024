
function m = bij2base(varargin)

m = [varargin{end-1:-1:1}]*(varargin{end}).^(0:length(varargin)-2)';

end