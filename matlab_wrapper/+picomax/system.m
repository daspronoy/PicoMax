



classdef system < handle & matlab.mixin.SetGet
    
    properties
        crystal {mustBeMember(crystal,{'fcc','hex'})} = 'fcc'
        a (1,1) {mustBeNumeric}
        lattice_vector (3,:) {mustBeNumeric}
        reciprocal_vector (3,:) {mustBeNumeric}

        epm {mustBeNumeric}

    end
    
    methods
        function o = system(varargin)
            % o.set(varargin{:});
        end

        function [] = initialize(o)
            switch o.crystal
                case 'fcc'
                    o.lattice_vector = [
                        [0;0.5;0.5], ...
                        [0.5;0;0.5], ...
                        [0.5;0.5;0]];
                    o.reciprocal_vector = [
                        [-1;1;1], ...
                        [1;-1;1], ...
                        [1;1;-1]];
                case 'hex'
            end

        end
    end
end

