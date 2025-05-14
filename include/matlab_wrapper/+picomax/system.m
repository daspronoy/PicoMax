



classdef system < handle & matlab.mixin.SetGet
    
    properties
        crystal = 'zb'
        a (1,1) {mustBeNumeric}
        f (1,1) {mustBeNumeric}
        u (1,1) {mustBeNumeric}
        lattice_vector (3,:) {mustBeNumeric}
        reciprocal_vector (3,:) {mustBeNumeric}

        epm {mustBeNumeric}
        vg {mustBeNumeric}
    end
    
    methods
        function o = system(varargin)
            % o.set(varargin{:});
        end

        function [] = initialize(o)
            switch o.crystal
                case 'zb'
                    o.lattice_vector = [
                        [0;0.5;0.5], ...
                        [0.5;0;0.5], ...
                        [0.5;0.5;0]];
                    o.reciprocal_vector = [
                        [-1;1;1], ...
                        [1;-1;1], ...
                        [1;1;-1]];
                case 'wz'
            end

        end
    end
end

