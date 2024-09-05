classdef PolynomialZonotope

    properties
        G
        E
        id
        dims
    end

    methods

        function obj = PolynomialZonotope(G,E,id)
            
            if nargin == 0
                help("PolynomialZonotope");
                obj = [];
                return;
            end

            if nargin >= 1
                if isa(G,'PolynomialZonotope')
                    obj.G = G.G;
                    obj.E = G.E;
                    obj.dims = G.dims;
                    obj.id = UniqueID(obj.dims.p);
                    return;
                end
                obj.G = G;
                obj.dims.n = size(G,1);
                obj.dims.h = size(G,2);
            end

            if nargin >= 2
                if size(E,2) ~= obj.dims.h
                    error('Incommensurate Dimensions between G and E');
                end
                obj.E = E;
                obj.dims.p = size(E,1);
            else
                obj.E = [0 ones(1,obj.dims.h-1)];
                obj.dims.p = 1;
            end

            if nargin >= 3
                if length(id) ~= obj.dims.p
                    error('Incommensurate Dimensions Specified for identifiers');
                end
                obj.id = id;
            else
                obj.id = UniqueID(obj.dims.p);
            end

        end

    end

end
    