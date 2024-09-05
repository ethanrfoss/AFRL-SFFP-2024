classdef Polytope

    properties
        Vertices
        Faces
    end

    methods

        function obj = Polytope(Vertices)

            obj.Vertices = Vertices(ismember(1:size(Vertices,1),convhulln(Vertices)),:);
            obj.Faces = convhulln(obj.Vertices);

        end

    end

end
    