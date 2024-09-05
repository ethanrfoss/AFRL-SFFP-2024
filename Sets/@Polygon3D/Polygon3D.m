classdef Polygon3D

    properties
        Vertices
        Faces
    end

    methods

        function obj = Polygon3D(Vertices,Faces)

            if size(Vertices,2) ~= 3
                error('Incorrect Dimensions Specified for Polygon');
            end
            obj.Vertices = Vertices;

            if nargin>=2
                obj.Faces = Faces;
            else
                obj.Faces = convhull(Vertices(:,1),Vertices(:,2),Vertices(:,3));
            end

        end

    end

end
    