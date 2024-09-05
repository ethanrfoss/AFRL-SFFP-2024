classdef Polygon2D

    properties
        Vertices
    end

    methods

        function obj = Polygon2D(Vertices)

            if size(Vertices,2) ~= 2
                error('Incorrect Dimensions Specified for Polygon');
            end
            
            obj.Vertices = Vertices;

        end

    end

end
    