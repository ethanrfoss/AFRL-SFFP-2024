classdef Zonotope

    properties
        c
        G
        id
        dims
    end

    methods

        function obj = Zonotope(c,G,id)
            
            if nargin == 0
                help("Zonotope");
                obj = [];
            end

            if nargin >= 1
                if isempty(c)
                    obj.c = zeros(size(G,1));
                    obj.dims.n = size(G,1);
                elseif size(c,2) == 1
                    obj.c = c;
                    obj.dims.n = size(c,1);
                else
                    error('Incommensurate Dimensions Specified for c');
                end
            end

            if nargin >= 2
                obj.G = G;
                if obj.dims.n ~= size(G,1)
                    error('Incommensurate Dimensions for G');
                end
                obj.dims.l = size(G,2);
            end

            if nargin >= 3
                if length(id) ~= obj.dims.l
                    error('Incommensurate Dimensions specified for zonotope');
                end
                obj.id = id;
            else
                obj.id = UniqueID(obj.dims.l);
            end

        end

    end

end
    