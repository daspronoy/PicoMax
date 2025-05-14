



classdef numeric < handle & matlab.mixin.SetGet
    
    properties
        % permittivity matrix dimension
        neps (1,1) {mustBeNumeric} = 1

        % plane wave sphere cutoff
        % TODO: change to ECUT [eV]
        encut (1,1) {mustBeNumeric} = 200.0;
        % gcut (1,1) {mustBeNumeric} = 4.0
        
        % number of electronic bands
        nband (1,1) {mustBeNumeric} = 8
        originshift (1,1) {mustBeInteger} = 0
        nvalence (1,1) {mustBeNumeric} = 4
        
        % number of q-vectors
        nqvec (1,1) {mustBeNumeric}
        qvec (3,:) {mustBeNumeric}
        qpath
        qnum
        gammazero (3,1) {mustBeNumeric} = [0;0;0]

        % number of frequency
        nfreq (1,1) {mustBeNumeric} = 1
        dfreq (1,1) {mustBeNumeric} = 0.1
        epsilon (1,1) {mustBeNumeric} = 0.1
        kk (1,1) {mustBeInteger} = 0
        delta (1,1) {mustBeInteger} = 0

        % BZ integration (KPOINT)
        kpointfile string
        kpoint (1,1) {mustBeNumeric} = 0
        kpointorder (1,1) {mustBeNumeric} = 3

        refpoint (3,1) {mustBeNumeric} = [0;0;0]
    end

    properties
        % debug order
        debug (1,1) {mustBeInteger} = 0

        %
        system {mustBeText} = 'wsl'
        
        base_remote {mustBeText} = '/'
        base_local {mustBeText} = '//wsl.localhost/ubuntu-22.04/'
        wdir {mustBeText} = 'tmp/'
        outputfile {mustBeText} = 'tmpepsilon'


        % PATH string
        % PATH_DAT {mustBeText} = '\\wsl.localhost\ubuntu-22.04\tmp\tmpepsilon.dat'
        % WORKDIR 
        % TMPDIR 

        % picomax exe path
        picomax_exe string = "/home/mjh92/projects/picomax/src/picomax"
        
        % force PicoMax (c++) if possible
        force_cpp {mustBeNumeric} = 0

        %
        NumHeaderLines (1,1) {mustBeInteger} = 3
    end
    
    methods
        function o = numeric(varargin)
            % o.set(varargin{:});
        end
    end
end

