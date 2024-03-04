



classdef numeric < handle & matlab.mixin.SetGet
    
    properties
        % permittivity matrix dimension
        NEPS (1,1) {mustBeNumeric} = 1

        % plane wave sphere cutoff
        % TODO: change to ECUT [eV]
        GCUT (1,1) {mustBeNumeric} = 4.0
        
        % number of electronic bands
        NBAND (1,1) {mustBeNumeric} = 8

        % number of q-vectors
        NQVEC (1,1) {mustBeNumeric}
        QVEC (3,:) {mustBeNumeric}
        QPATH
        QNUM

        % number of frequency
        NFREQ (1,1) {mustBeNumeric} = 1
        DFREQ (1,1) {mustBeNumeric} = 0.1
        EPSILON (1,1) {mustBeNumeric} = 0.1
        KK (1,1) {mustBeInteger} = 0
        
        % BZ integration (KPOINT)
        KPOINT_ORDER (1,1) {mustBeNumeric} = 3
        KPOINT_SCHEME (1,1) {mustBeNumeric} = 0
        KPOINT_FILE string
    end

    properties
        % debug order
        DEBUG (1,1) {mustBeInteger} = 0

        %
        SYSTEM {mustBeText} = 'wsl'

        PATH string
        PATH_DAT {mustBeText} = '\\wsl.localhost\ubuntu-22.04\tmp\tmpepsilon.dat'
        WORKDIR 
        TMPDIR 

        % 
        PICOMAX_CMD string = "/home/mjh92/projects/picomax_try/src/picomax"
        
        % force PicoMax (c++) if possible
        FORCE_CPP {mustBeNumeric} = 0

        %
        NumHeaderLines (1,1) {mustBeInteger} = 3
    end
    
    methods
        function o = numeric(varargin)
            % o.set(varargin{:});
        end
    end
end

