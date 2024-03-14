%Matlab wrapper for picomax


classdef picomax < handle & matlab.mixin.SetGet
    
    properties
        dat
        opt
        num
        sys
        tmp
    end
    
    methods
        function o = picomax(varargin)

            o.sys = picomax.system;
            o.num = picomax.numeric;
        end
        
        function [] = set(o,varargin)
            set(o.sys,varargin{:});
        end

        function [] = setopt(o,varargin)
            set(o.num,varargin{:});
        end

        function [] = run(o,option)
            
            if strcmp(option,'band') && o.num.force_cpp
                picomax.electron.band_cpp(o);

            elseif strcmp(option,'band')
                picomax.electron.band(o);

            elseif strcmp(option,'eps')
                o.permittivity;
            end
        end

        function [] = import(o,filename,option)

            % load data
            fid = fopen(filename);
            data = fread(fid,'double');
            fclose(fid);

            % parse data
            if strcmp(option,'eps')
                nfreq = o.num.nfreq;
                o.dat.omega = o.num.dfreq*(0:1:nfreq-1);
                if (~isempty(o.num.qnum) && ~isempty(o.num.qpath))
                    nq = sum(o.num.qnum)+1;
                    o.dat.qi = 0:nq-1;
                else
                    nq = 1;
                    o.dat.qi = 0;
                end
                ne = o.num.neps*(o.num.neps+1)/2;
                N = ne*nq*nfreq;

                o.dat.realepsL = permute(reshape(data(0*N+1:1*N),[nfreq ne nq]),[1 3 2]);
                o.dat.imagepsL = permute(reshape(data(1*N+1:2*N),[nfreq ne nq]),[1 3 2]);
                o.dat.realepsT = permute(reshape(data(2*N+1:3*N),[nfreq ne nq]),[1 3 2]);
                o.dat.imagepsT = permute(reshape(data(3*N+1:4*N),[nfreq ne nq]),[1 3 2]);

            elseif strcmp(option,'eband')
                
            end
        end

        function [k,f] = highsym_point(~,point,crystal)
            if strcmp(crystal,'fcc')
                [k,f] = picomax.highsym.fcc(point);
            elseif strcmp(crystal,'hex')
                [k,f] = picomax.highsym.hex(point);
            end
        end
    end
end

