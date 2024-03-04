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
            
            if strcmp(option,'band') && o.num.FORCE_CPP
                picomax.electron.band_cpp(o);

            elseif strcmp(option,'band')
                picomax.electron.band(o);

            elseif strcmp(option,'eps')
                o.permittivity;
            end
        end

           
        function [] = compute_epsM(o)
            % Compute macroscopic dielectric constants

            epsmnL = o.dat.realepsL+1i*o.dat.imagepsL;
            epsmnT = o.dat.realepsT+1i*o.dat.imagepsT;
            

            NEPS = o.neps;
            NFREQ = o.nfreq;
            NQ = sum(o.numq)+1;

            if NQ==1 && norm(o.q)==0
                epsmnL(:,1,:,:) = epsmnT(:,1,:,:);
            end


            invepsmnL = zeros(NFREQ,NQ,NEPS,NEPS);
            eigvalue_L = zeros(NFREQ,NQ,NEPS);
            % invepsmnT = zeros(NFREQ,NQ,NEPS,NEPS);
            eigvalue_T = zeros(NFREQ,NQ,NEPS);

            for f=1:NFREQ
                for q=1:NQ
                    epsmn_L = zeros(NEPS);
                    for m=1:NEPS
                        for n=1:m
                            idx = (m-1)*m/2+n;
                            epsmn_L(m,n) = epsmnL(f,q,idx);
                            if (m~=n)
                                epsmn_L(n,m) = conj(epsmn_L(m,n));
                            end
                        end
                    end
                    invepsmnL(f,q,:,:) = inv(epsmn_L);
                    tmpeig = eig(epsmn_L);
                    for m=1:NEPS
                        eigvalue_L(f,q,m) = tmpeig(m);
                    end


                    epsmn_T = zeros(NEPS);
                    for m=1:NEPS
                        for n=1:m
                            idx = (m-1)*m/2+n;
                            epsmn_T(m,n) = epsmnT(f,q,idx);
                            if (m~=n)
                                epsmn_T(n,m) = conj(epsmn_T(m,n));
                            end
                        end
                    end
                    % invepsmnT(f,:,:) = inv(epsmn_T);
                    tmpeig = eig(epsmn_T);
                    for m=1:NEPS
                        eigvalue_T(f,q,m) = tmpeig(m);
                    end
                end
            end
            
            o.dat.epsilon00 = epsmnL(:,:,1,1);
            o.dat.epsilonM = 1./invepsmnL(:,:,1,1);
            o.dat.alpha_L = eigvalue_L;
            o.dat.alpha_T = eigvalue_T;


            
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

