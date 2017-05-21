classdef stkSpecObj < SpecObj
    properties (SetAccess = protected)
        chi = [];       %chi angle, calculated with pol and inc
        field = 1000;   %applied field in Vrms
        pol = 90;       %polarization, horizontal -> 90 deg, vertical -> 0 deg
        g = [];         %Achi*g
        dg = [];        %Bchi*dg
        ddg = [];       %Cchi*ddg
    end
    methods
    % class constructor -------------------------------------------------------------
        function obj = stkSpecObj(), return; end
    % property SET methods ----------------------------------------------------------
        function setField(obj,val), obj.field = val; obj.toExtinct; end
        function setPol(obj,val), obj.pol = val; obj.calcChi; obj.toExtinct; end
        function setDerivatives(obj,g,dg,ddg)
            obj.g = g;
            obj.dg = dg;
            obj.ddg = ddg;
        end
    % -------------------------------------------------------------------------------
        function loadSpec(obj,f)
        % ---------------------------------------------------------------------------
        % Loads either LV binary or LV-created Matlab spectrum files from stepscan
        % ---------------------------------------------------------------------------

        % Determine if file/path was passed as optional argument
            if nargin < 2
            % Get filename/path
                [fname,pname] = uigetfile('*.*','Select spectrum file');
                if fname == 0
                    return; 
                else
                    f = [pname,fname];
                end;
            end
            [~,tok] = strtok(f,'.');
        % If no extension load as binary, else load Matlab
            if numel(tok) == 0
                try
                    [X,I,I0] = stkSpecObj.readLvBinary(f);
                catch
                    errordlg('File does not contain Stark spectrum','Spec Error');
                    X = []; I = []; I0 = [];
                end
            elseif strcmp(tok,'.mat')
                    [X,I,I0] = stkSpecObj.readLvMatlab(f);
            else
                errordlg('File does not contain Stark spectrum','Spec Error');
                return;
            end
            obj.data = [X,I,I0];
        % Parse path/file -> pname, fname
            c = filesep;
            [pname,tok] = strtok(f,c);
            while numel(tok) ~= 0
                [addpath,tok] = strtok(tok,c);
                if numel(tok) ~= 0, pname = [pname,c,addpath]; end
            end
            obj.name = {pname,addpath};
        % Select last 10% of data points for baseline
            lenx = length(X);
            obj.bsln = ceil(1e7./[X(lenx-ceil(lenx*0.1)) X(lenx)]);
        % Normalize to field and process data
            obj.calcChi;
            obj.toExtinct;
            obj.correctBaseline;
        end
    % -------------------------------------------------------------------------------
        function fieldNorm(obj,extinct)
        % ---------------------------------------------------------------------------
        % Normalizes delta(Extinction) spectrum to 10^6 V/cm (8/23/13). This is
        % roughly the field applied in the actual experiment for typical
        % voltages and pathlengths.
        % INPUT ARGUMENTS:
        %   extinct ->  spectrum in extinction (M^-1 cm^-1)
        %
        % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        % WHEN YOU CHANGE THE VALUE OF THE NORMALIZATION FIELD YOU MUST
        % UPDATE THIS VALUE IN THE FITTING FUNCTIONS OF stkFitObj!!!!!!!!!!
        % IN THE STATIC METHODS simfit AND stkfit
        % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        %
        % ---------------------------------------------------------------------------

        %Field Vrms -> Vpk
            V = (obj.field*2.828)/2;
        %Pathlength um -> cm
            b = obj.path*1e-4;
            E2 = (V/b)^2;
        %Normalize to Selected field in (V/cm)^2, no f (local field factor)
            obj.y = (extinct/E2)*(1e6^2); %(1 MV/cm)^2 8/23/13
        end
    % -------------------------------------------------------------------------------
        function calcChi(obj)
        % Calculates chi based on solvent refractive index
            n1 = 1.2; 
            theta1 = obj.inc*pi/180;
            theta2 = asin(n1/obj.n2*sin(theta1));
            if obj.pol == 0
                obj.chi = 90;
            else
                obj.chi = 90-theta2*180/pi;
            end
        end
    end
end