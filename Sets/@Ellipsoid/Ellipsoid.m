classdef Ellipsoid

    properties
        c
        Q
    end

    methods

        function obj = Ellipsoid(Q,c)
            
            if nargin == 0
                help('Ellipsoid');
            end
            
            if size(Q,1) ~= size(Q,2)
                error('Q must be Square');
            end
            try
                chol(Q);
            catch
                error('Q Must be Positive Definite');
            end
            obj.Q = Q;

            if nargin > 1
                if size(c,1)~=size(Q,1)
                    error('Invalid Dimension for c');
                end
                obj.c = c;
            else
                obj.c = zeros(size(Q,1),1);
            end

        end

    end

end
    