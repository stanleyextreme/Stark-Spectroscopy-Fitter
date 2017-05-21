classdef ltaSpecObj < SpecObj
    methods
    % class constructor -----
        function obj = ltaSpecObj(), return; end
    % -------------------------------------------------------------------------------
        function loadSpec(obj,sigtype,f)
        %----------------------------------------------------------------------------
        %Loads either LV binary or LV-created Matlab spectrum files from stepscan
        %for ltaSpecObj
        %INPUT ARGUMENTS:
        %   sigtype ->  (optional) STRING spectrum type 'sig' or 'ref'   
        %   f ->        (optional) STRING path/filename
        %----------------------------------------------------------------------------

        % Determine if file/path was passed as optional argument
            if nargin == 1 %no options selected
                sigtype = 'sig';
            elseif nargin == 2 %make sure sigtype is acceptable value
                if (strcmp(sigtype,'sig') | strcmp(sigtype,'ref')) == 0
                    sigtype = 'sig';
                    disp('Unrecognized signal type. SIGNAL defaulted');
                end
            end
            if nargin ~= 3 %path/file not given, get it
                [fname,pname] = uigetfile('*.*','Select spectrum file');
                if fname == 0 %cancel button
                    return; 
                else
                    f = [pname,fname];
                end;
            end
            [~,tok] = strtok(f,'.');
        % If no extension load as binary, else load Matlab
            if numel(tok) == 0
                try
                    [X,I,~] = stkSpecObj.readLvBinary(f);
                catch
                    errordlg('File does not contain spectrum','Spec Error');
                    X = []; I = [];
                end
            elseif strcmp(tok,'.mat')
                [X,I,~] = stkSpecObj.readLvMatlab(f);
                I = mean(I,2); %for 2014 file format 6/16/14
            else
                errordlg('File does not contain spectrum','Spec Error');
                return;
            end
        % Select last 10% of data points for baseline
            lenx = length(X);
            obj.bsln = ceil(1e7./[X(lenx-ceil(lenx*0.1)) X(lenx)]);
        % Parse path/file -> pname, fname
            c = filesep;
            [pname,tok] = strtok(f,c);
            while numel(tok) ~= 0
                [addpath,tok] = strtok(tok,c);
                if numel(tok) ~= 0, pname = [pname,c,addpath]; end
            end
        %Store data in appropriate fields
            if strcmp(sigtype,'sig')
                if numel(obj.data) == 0
                    obj.data = [X,I,zeros(size(X))];
                else
                %if data was previously loaded, check that wavenumber range and
                %stepsize are the same
                    if unique(obj.data(:,1) == X)
                        obj.data(:,2) = I;
                    else
                        errordlg('Signal and Reference sizes incommensurate','Spec Error');
                        return
                    end
                end
                obj.name{1,1} = pname;
                obj.name{1,2} = addpath;
            elseif strcmp(sigtype,'ref')
                if numel(obj.data) == 0
                    obj.data = [X,zeros(size(X)),I];
                else
                    if unique(obj.data(:,1) == X)
                        obj.data(:,3) = I;
                    else
                        errordlg('Signal and Reference sizes incommensurate','Spec Error');
                        return
                    end
                end
                obj.name{2,1} = pname;
                obj.name{2,2} = addpath;
            end
        %Convert absorbance to extinction
            obj.toExtinct;
            obj.correctBaseline;
        end
    % -------------------------------------------------------------------------------
    end
end