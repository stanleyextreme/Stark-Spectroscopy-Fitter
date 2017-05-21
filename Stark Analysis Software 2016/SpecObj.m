classdef SpecObj < handle
    properties (SetAccess = protected)
        path = 55;          %pathlength in um
        conc = 1;           %concentration in mM
        inc = 45;           %angle of incidence, normal = 0 deg
        bsln = [];          %baseline limits
        n2 = 1.4;           %refractive index of solvent
        data = [];          %container for raw data [X,I,I0]
        x = [];             %processed x-data
        y = [];             %processed y-data
        name = {};          %file path and name
        gparms = [];        %matrix of gaussian parameters [amp,center,width]
        bands = [];         %vector of transition assignments
        thefit              %simulated spectrum
        gau                 %gaussians
        prog                %progressions
    end
    methods
    % class constructor -------------------------------------------------------------
        function obj = SpecObj(), return; end
    % property SET methods ----------------------------------------------------------
        function setPath(obj,val), obj.path = val; obj.correctBaseline; end
        function setConc(obj,val), obj.conc = val; obj.correctBaseline; end
        function setInc(obj,val), obj.inc = val; obj.correctBaseline; end
        function setRefInd(obj,val), obj.n2 = val; obj.correctBaseline; end
        function setBaseline(obj,val), obj.bsln = val; obj.correctBaseline; end
        function setGauParms(obj,val), obj.gparms = val; end
        function setBands(obj,val), obj.bands = val; end
    % -------------------------------------------------------------------------------
        function setFit(obj,f,g,p)
            obj.thefit = f; obj.gau = g; obj.prog = p;
        end
    % -------------------------------------------------------------------------------
        function loadSpec(obj,x,y)
            obj.x = x;
            obj.y = y;
            if size(x,1) < size(x,2), x = x'; end;
            if size(y,1) < size(y,2), y = y'; end;
            obj.data = [x,y];
        end
    % -------------------------------------------------------------------------------
        function toExtinct(obj)
        % ---------------------------------------------------------------------------
        % Converts delta(Absorption) to delta(Extinction) based on adjusted 
        % pathlength.
        % ---------------------------------------------------------------------------
        % Correct pathlength for angle
            n1 = 1.2; 
            theta1 = obj.inc*pi/180;
            theta2 = asin(n1/obj.n2*sin(theta1));
        % Pathlength -> cm
            bcorr = obj.path*1e-4;
            bcorr = bcorr/cos(theta2);
        % Concentration -> M
            c = obj.conc*1e-3;
        % Convert to extinction and normalize to E = 300,000 V/cm
            if strcmp(class(obj),'stkSpecObj')
                if size(obj.data,2) == 3
                    A = (2*sqrt(2).*obj.data(:,2))./((2.303*bcorr*c).*obj.data(:,3));
                elseif size(obj.data,2) > 3 %for 2014 data format 6/16/14
                    d = obj.data; 
                    d(:,1) = []; %get rid of wavenumbers
                    cols = size(d,2)/2; %how many columns for i or i0?
                    i = d(:,1:cols);
                    i0 = d(:,cols+1:size(d,2));
                    A = (2*sqrt(2).*i)./((2.303*bcorr*c).*i0); % delta I/I0 gives delta extinction, energy-weighted
                    A = mean(A,2);
                else 
                    A = obj.data(:,2); %spectrum is already processed to extinction
                end
                obj.fieldNorm(A); %A is extinction normalize to 1 MV/cm
            elseif strcmp(class(obj),'ltaSpecObj')
                %if nothing is loaded into reference
                if length(unique(obj.data(:,3))) ~= 1
                    A = log10(obj.data(:,3)./obj.data(:,2));
                    obj.y = A/(bcorr*c);
                else
                    obj.y = obj.data(:,2); %spectrum is already processed to extinction
                end
            else
                disp('what kind of spectrum is this?');
            end
            obj.x = obj.data(:,1);
        end
    % -------------------------------------------------------------------------------
        function correctBaseline(obj)
        %----------------------------------------------------------------------------
        %Baseline correction
        %Correction for Stark spectra sets the mean value of the selected
        %baseline region to zero (won't handle slanted baseline well, but
        %this is discouraged because there is no evidence that a messed up
        %baseline in Stark spectra is linear) Correction for LT absorption
        %spectra is a linear fit to the selected baseline region
        %----------------------------------------------------------------------------
            obj.toExtinct; %refresh data
            if numel(obj.bsln) == 0, return; end
            if max(obj.bsln == 0) == 1, return; end;
            wl = 1e7./obj.x;
            noise = (wl < max(obj.bsln) & wl > min(obj.bsln));
            if max(noise) ~= 0
                if strcmp(class(obj),'stkSpecObj')
                    obj.y = obj.y-mean(obj.y(noise));
                elseif strcmp(class(obj),'ltaSpecObj')
                    p = polyfit(obj.x(noise),obj.y(noise),1);
                    obj.y = obj.y-polyval(p,obj.x);
                end
            else
                errordlg('No data within specified spectral region','Baseline Error')
            end
        end
    % -------------------------------------------------------------------------------
        function [g,dg,ddg] = calcGau(obj)
        %separate matrix into respective vectors
            amp = obj.gparms(:,1); mu = obj.gparms(:,2); sig = obj.gparms(:,3);            
        %resize to appropriate matrices
            [amp,x] = meshgrid(amp,obj.x);
            [mu,x] = meshgrid(mu,obj.x);
            [sig,x] = meshgrid(sig,obj.x);
            sig = 2.*sig.^2;
        %calculate gaussians & derivatives
            g = (amp.*exp(-(x-mu).^2./(sig)))./x; %individual guassians
            dg = -g.*(1./x+(x-mu).*2./(sig)); %2st derivatives
            ddg = -dg.*(2.*(x-mu)./(sig)+1./x)-g.*(2./(sig)-1./x.^2); %2nd derivative
        end
    end
    methods (Static)
        function [X,I,I0] = readLvBinary(f)
        % ---------------------------------------------------------------------------
        % Converts LabView binary files recorded with LV5 stepscan or LV7.1
        % stepscan2013 to data usable in MatLab. Data arrays are delimited by a
        % pattern of (null null null var) where var is a variable byte between
        % files, but consistent throughout the selected file. Only the wavenumber,
        % Lock-In signal and NI-DAQ signal arrays are read and output
        % INPUT ARGUMENTS:
        %   f   -> path/filename to convert
        % OUTPUT ARGUMENTS:
        %   X   -> wavenumber
        %   I   -> Lock-In signal
        %   I0  -> NI-DAQ signal
        % ---------------------------------------------------------------------------

        % Open binary file
            fid = fopen(f,'r','b');
            A = fread(fid,inf,'uint8','b');
            fclose(fid);
        % Find data start point (null null null ?)
            y = find(A == 0); 
            a1 = A(y); 
            a2 = A(y+1);
            a3 = A(y+2);
            a4 = A(y+3);
        % set index and locations
            idx = (a1 == 0 & a2 == 0 & a3 == 0 & a4 == A(4));
            loc = y(idx == 1);
        %Make sure location tag isn't part of data
            dd = diff(loc);
            i = 2;
            while i < length(dd)
                if dd(i) < dd(1)
                    loc(i+1) = [];
                    dd = diff(loc);
                    i = 2;
                else
                    i = i+1;
                end
            end
        %Read I
            i = loc(1); j = loc(2);
            y = A(i+4:j-1);
            lnth = length(y);
            I = hex2num(reshape(dec2hex(y)',16,[])');
        %Read I0
            i = loc(2); j = loc(3);
            y2 = A(i+4:j-1);
            I0 = hex2num(reshape(dec2hex(y2(1:lnth))',16,[])');
        %Read wavenum
            i = loc(5); j = loc(6);
            y5 = A(i+4:j-1);
            X = hex2num(reshape(dec2hex(y5(1:lnth))',16,[])');
        end
        function [X,I,I0] = readLvMatlab(f)
        % ---------------------------------------------------------------------------
        % Reads .mat files created in LV7.1 stepscan2013
        % INPUT ARGUMENTS:
        %   f   -> path/filename to read
        % OUTPUT ARGUMENTS:
        %   X   -> wavenumber
        %   I   -> Lock-in signal
        %   I0  -> NI-DAQ signal
        % ---------------------------------------------------------------------------

        % load file
            try
                %load(f,'i','i0');
                D = load(f);
                if size(D.i,2) > size(D.i,1)
                    D.i = D.i'; D.i0 = D.i0'; %make sure data is in column format 6/16/14
                end
                i = D.i; i0 = D.i0;
                %X = i(1,:)';
                %I = i(2,:)';
                %I0 = i0(2,:)';
                %CHANGED TO READ NEWEST DATA FORMAT 6/16/14
                X = i(:,1);
                I = i(:,2:size(i,2)); 
                I0 = i0(:,2:size(i,2));
            catch
                try %CHANGED TO READ PREPROCESSED SPECTRA 8/23/13
                    load(f,'extinct');
                    X = extinct(:,1);
                    I = extinct(:,2);
                    I0 = [];
                catch
                    errordlg('File does not contain appropriate arrays','Spec Error');
                    X = []; I = []; I0 = [];
                end
            end
        end
    end
end