classdef (Abstract) Basis
    %BASIS Virtual class for Basis construction
    
    properties (Abstract)
        size;
        centers;
        Quad;
    end
    
    methods (Static)
        [X, W] = hermquad(N);
        [X, W, reccals] = gengausslegquadrule(nmax,tol,numpoint,reccalls,maxrec)
    end
    
    methods (Abstract)
        result = Evaluation(obj, Grid, R, ind);
        result = Overlap(obj, R1, R2);
        result = W(obj, V, patch);
    end
end

