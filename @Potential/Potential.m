classdef Potential
    %PERTURBINGPOTENTIAL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        lattice
        K_list
        FourierCoefficients
    end
    
    methods
        function obj = Potential(lattice, K_list, FourierCoefficients)
            if nargin >0
                obj.lattice = lattice;
                obj.K_list = K_list;
                obj.FourierCoefficients = FourierCoefficients;
            end
        end
        
        function result = Evaluate(obj,z,k)
            result = obj.FourierCoefficients{k}(z);
        end
    end
end

