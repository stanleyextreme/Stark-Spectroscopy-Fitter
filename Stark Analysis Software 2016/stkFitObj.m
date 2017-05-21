classdef stkFitObj < handle
    properties %(SetAccess = protected)
    % n = #transitions, m = #gaussians
        sobj = [];              %array of stkSpecObj instances
        aobj = [];              %ltaSpecObj instance
        results = [];           %nx4 matrix of tralpha, projalpha, deltamu, zeta
        Achi = [];              %nxspec matrix of Achi values
        bands = [];             %mx1 vector of band assignments
        gauparms = [];          %mx4 matrix of amp,mu,sigma
        opts = [];              %struct:MaxFunEvals,TolFun,TolX,weight
        lims = [300 650];       %wavelength range to fit
        ax = [];                %axes for plotting (1) Stk (2) Abs (3) Res (4) Opt
        kt = 1;                 %counter
        mc = [];                %Monte Carlo results
    end
    methods
    % class constructor -------------------------------------------------------------
        function obj = stkFitObj()
            obj.aobj = ltaSpecObj;
            obj.opts.MaxEval = 2000; %set fitting options
            obj.opts.weight = 0.1; %absorption weighting
            obj.opts.v = 8; %gaussians/progression
            obj.opts.TolFunAbs = 1e-6;
            obj.opts.TolFunSim = 1e-12;
            obj.opts.TolX = 1e-12;
            obj.opts.weightval = 1; %actual weighting factor
            obj.opts.MCiter = 50;
            obj.opts.comments = '';
            obj.mc.results = '';
        end     
    % -------------------------------------------------------------------------------
        function addStark(obj)
            tmpobj = stkSpecObj;
            tmpobj.loadSpec;
            if numel(tmpobj.name) ~= 0
                obj.sobj = [obj.sobj,tmpobj];
            end
        end
    % -------------------------------------------------------------------------------
        function delStark(obj,index)
            try
                obj.sobj(index) = [];
            catch exception
                errordlg('No object at specified index','Delete Error');
            end
        end
    % -------------------------------------------------------------------------------
        function xx = setRange(obj,specobj)
        %returns logical vector of indices that fall within the selected
        %limits of specobj, either obj.aobj or obj.sobj(j)
            if max(obj.lims) < 10000 %limits in wavelength
                lims = 1e7./obj.lims;
            else
                lims = obj.lims;
            end
            xx = (specobj.x < max(lims) & specobj.x > min(lims));
        end
    % -------------------------------------------------------------------------------
        function fitAbs(obj)
            obj.kt = 1;
            lam = obj.gauparms;
            lb = zeros(size(lam));
            ub = ones(size(lam)).*Inf;
            opts = optimset('lsqcurvefit');
            opts = optimset('Display','off','MaxFunEvals',obj.opts.MaxEval,...
                'TolFun',obj.opts.TolFunAbs,'TolX',obj.opts.TolX);
            y = obj.aobj.y(obj.setRange(obj.aobj))./obj.aobj.x(obj.setRange(obj.aobj));
            x = obj.aobj.x(obj.setRange(obj.aobj));simfit
            
            [beta,~,~] = lsqcurvefit(@stkFitObj.absfit,lam,x,y,lb,ub,opts,obj);
            obj.kt = 1;
            obj.gauparms = beta;
        end
    % -------------------------------------------------------------------------------
        function fitSim(obj)
            %called when the GUI Fit button is pressed
            obj.kt = 1;
            lamconvert = [reshape(obj.gauparms,[],1);...
                reshape(obj.results,[],1);...
                reshape(obj.Achi,[],1)];
        %set upper and lower bounds for fitted parameters
            lb = [zeros(size(reshape(obj.gauparms,[],1)));...
                reshape([ones(size(obj.results(:,1:2))).*0,zeros(size(obj.results(:,3:4)))],[],1);... %first set should be -Inf
                reshape(ones(size(obj.Achi)).*-Inf,[],1)]; 
            ub = [ones(size(reshape(obj.gauparms,[],1))).*Inf;...
                reshape(ones(size(obj.results(:,1:3))).*Inf,[],1);...
                reshapesimf(ones(size(obj.results(:,4))).*90,[],1);...
                reshape(ones(size(obj.Achi)).*1e-20,[],1)]; 
            ub = abs(ub./lamconvert);
        %all parameters for lsqcurvefit need to be one, and then converted
        %to appropriate values within fitting function due to variance in
        %magnitude of the parameters
            lam = ones(size(lamconvert));
        %absorption x & y(energy-weighted) data
            ax = obj.aobj.x(obj.setRange(obj.aobj));
            ay = obj.aobj.y(obj.setRange(obj.aobj)); ay = ay./ax;
        %stark x & y(energy-weighted) data
            sx = [obj.sobj.x]; ll = obj.setRange(obj.sobj(1)); sx = sx(ll,:);
            sy = [obj.sobj.y]; sy = sy(ll,:); sy = sy./sx;
        %absorption weighting
            %obj.opts.weightval = (trapz(abs(ay))./sum(trapz(abs(sy)),2))*obj.opts.weight;
            obj.opts.weightval = (trapz(abs(ay))./(sum(trapz(abs(sy)),2))./size(sy,2))*obj.opts.weight;
        %create column vector of Stark and Abs data
            x = [ax;reshape(sx,[],1)];
            y = [ay./obj.opts.weightval;reshape(sy,[],1)];
        %set fitting options
            opts = optimset('lsqcurvefit');
            opts = optimset('Display','off','MaxFunEvals',obj.opts.MaxEval,...
                'TolFun',obj.opts.TolFunSim,'TolX',obj.opts.TolX);
        %do the fit!
            [beta,resnorm] = lsqcurvefit(@stkFitObj.simfit,lam,x,y,lb,ub,opts,obj,lamconvert);
            beta = beta.*lamconvert;
        %separate fitted parameters and store in appropriate objects
            split1 = numel(obj.gauparms);
            split2 = numel(obj.results);
            obj.gauparms = reshape(beta(1:split1),[],3);
            obj.results = reshape(beta(split1+1:split1+split2),[],4);
            obj.Achi = reshape(beta(split1+split2+1:length(beta)),max(obj.bands),[]);
        end
    % -------------------------------------------------------------------------------
        function fitStark(obj)
        %linear fit for B and C coefficients using previously fitted
        %gaussians for the absorption spectrum. used only for confidence
        %interval estimation, although I'm not sure how appropriate this
        %method is. Fitted parameters and CI are printed to command line,
        %not stored
            obj.kt = 1;
            lam = reshape(obj.results,[],1);
        %stark x & y(energy-weighted) data
            sx = [obj.sobj.x]; ll = obj.setRange(obj.sobj(1)); sx = sx(ll,:);
            sy = [obj.sobj.y]; sy = sy(ll,:); sy = sy./sx;
        %create column vector of data
            x = reshape(sx,[],1);
            y = reshape(sy,[],1);
        %fitting options
            opts = optimset('lsqcurvefit');
            opts = optimset('Display','on','MaxFunEvals',obj.opts.MaxEval,...
                'TolFun',obj.opts.TolFunSim,'TolX',obj.opts.TolX);
        %nonlinear regression
            [beta,r,J,COVB,mse] = nlinfit(x,y,@(lam,x) stkFitObj.stkfit(lam,x,obj),lam);
            beta
            ci = nlparci(beta,r,'covar',COVB)
        end
    % -------------------------------------------------------------------------------
        function tryFit(obj,pl)
            if nargin == 1, pl = 1; end; %plot fit
            if pl == 1
                obj.kt = 1;
            else
                obj.kt = 2;
            end
            if numel(obj.gauparms) ~= 0
                lam = [reshape(obj.gauparms,[],1);...
                    reshape(obj.results,[],1);...
                    reshape(obj.Achi,[],1)];
            %absorption x & y(energy-weighted) data
                ax = obj.aobj.x(obj.setRange(obj.aobj));
                ay = obj.aobj.y(obj.setRange(obj.aobj)); ay = ay./ax;
            %stark x & y(energy-weighted) data
                sx = [obj.sobj.x]; ll = obj.setRange(obj.sobj(1)); sx = sx(ll,:);
                sy = [obj.sobj.y]; sy = sy(ll,:); sy = sy./sx;
            %absorption weighting
                obj.opts.weightval = (trapz(abs(ay))./(sum(trapz(abs(sy)),2))./size(sy,2))*obj.opts.weight;
            %create column vector of data
                x = [ax;reshape(sx,[],1)];
                f = stkFitObj.simfit(lam,x,obj,ones(size(lam)));
            else
                obj.plotData;
            end
        end
    % -------------------------------------------------------------------------------
        function plotData(obj)
        %determine plotting range
            if numel(obj.sobj) == 0 && numel(obj.aobj.x) ~= 0
                xl = [min(obj.aobj.x) max(obj.aobj.x)];
            elseif numel(obj.sobj) ~= 0 && numel(obj.aobj.x) == 0
                xl = [min(min([obj.sobj.x])) max(max([obj.sobj.x]))];
            elseif numel(obj.sobj) ~= 0 && numel(obj.aobj.x) ~= 0
                xl = [min([min(obj.aobj.x),min(min([obj.sobj.x]))]) max([max(obj.aobj.x),max(max([obj.sobj.x]))])];
            else
                return
            end
        %absorption data- energy-weighted
            if numel(obj.aobj.y) ~= 0
                cla(obj.ax(2));
                x = obj.aobj.x; y = obj.aobj.y; y = (y./x)./obj.opts.weightval;
                ll = (x > min(1e7./obj.lims) & x < max(1e7./obj.lims)); %indices within limits
                plot(obj.ax(2),x,y,'color',[0.8 0.8 0.8],'linestyle','-','marker','none');
                plot(obj.ax(2),x(ll),y(ll),'k');
            %axes limits
                xlim(obj.ax(2),xl); xlim(obj.ax(3),xl);
                yl = max([abs(min(y)*.1) abs(max(y)*.1)]);
                yl = [min(y)-yl max(y)+yl];
                if isfinite(yl) == 1, ylim(obj.ax(2),yl); end
                abs_data = y;
            end
        %absorption fit
            if numel(obj.aobj.thefit) ~=0
                x = obj.aobj.x; y = obj.aobj.thefit./obj.opts.weightval;
                ll = (x > min(1e7./obj.lims) & x < max(1e7./obj.lims)); %indices within limits
                plot(obj.ax(2),x(ll),y,'Color','r','LineWidth',2);
                gg = obj.aobj.gau./obj.opts.weightval;
                plot(obj.ax(2),x(ll),gg);
                cols = get(gca,'ColorOrder');
                for j = 1:max(obj.bands)
                    plot(obj.ax(2),x(ll),gg(:,obj.bands == j),'Color',cols(j,:));
                    plot(obj.ax(2),x(ll),sum(gg(:,obj.bands == j),2),'Color',cols(j,:),'LineWidth',2);
                end
            %residuals
                cla(obj.ax(3));plot(obj.ax(3),x(ll),abs_data(ll)-y,'k');
                resstd = std(abs_data(ll)-y,0);
            end
        %stark data
            if numel(obj.sobj) ~= 0
                x = [obj.sobj.x]; y = [obj.sobj.y]; y = y./x; 
                ll = (x(:,1) > min(1e7./obj.lims) & x(:,1) < max(1e7./obj.lims)); %indices within limits
                cla(obj.ax(1));
                plot(obj.ax(1),x,y,'color',[0.8 0.8 0.8],'linestyle','-','marker','none'); %(~ll,:)
                plot(obj.ax(1),x(ll,:),y(ll,:));
                plot(obj.ax(1),x(:,1),zeros(size(x(:,1))),'k');
            %axes limits
                yl = max(abs([min(min(y(ll,:))) max(max(y(ll,:)))])*.1);
                yl = [min(min(y(ll,:)))-yl max(max(y(ll,:)))+yl];
                xlim(obj.ax(1),xl);
                ylim(obj.ax(1),yl);
                stk_data = y;
                if numel(obj.sobj.thefit) ~=0
                    for j = 1:length(obj.sobj)
                        x = obj.sobj(j).x; y = obj.sobj(j).thefit; 
                        ll = (x > min(1e7./obj.lims) & x < max(1e7./obj.lims)); %indices within limits
                        cols = get(gca,'ColorOrder');
                        plot(obj.ax(1),x(ll),y,'Color',cols(j,:),'LineWidth',2);
                        plot(obj.ax(3),x(ll),stk_data(ll,j)-y,'Color',cols(j,:));
                        resstd(j+1) = std(stk_data(ll,j)-y,0);
                    end
                    plot(obj.ax(3),x(ll),zeros(size(x(ll))),'k');
                    ylim(obj.ax(3),[-2*max(resstd) 2*max(resstd)]);
                    xt = get(obj.ax(3),'XTick'); xt = 1e7./xt;
                    set(obj.ax(3),'XTickLabel',num2str(floor(xt)'));
                    set(obj.ax(2),'XTickLabel',num2str(floor(xt)'));
                end
            end
            drawnow;
        end
    % -------------------------------------------------------------------------------
        function fitMonteCarlo(obj,lbl)
            tic; obj.mc = [];
            originalResults = obj.results;
            originalAchi = obj.Achi;
            originalGau = obj.gauparms;
            for j = 1:obj.opts.MCiter
                obj.mc.ri(:,:,j) = originalResults + originalResults.*(randn(size(obj.results)).*0.5); 
                obj.mc.gi(:,:,j) = originalGau + (ones(size(originalGau)).*1000).*(randn(size(originalGau)).*0.5); 
                obj.mc.ai(:,:,j) = obj.Achi;
                obj.results = obj.mc.ri(:,:,j);
                obj.Achi = originalAchi;
                obj.gauparms = originalGau;
                if j == 1
                    str = sprintf(['Fitting: ',num2str(j)]);
                else
                    str = sprintf(['Fitting: ',num2str(j),'\n chi^2 (',num2str(j-1),'): ',num2str(chi2)]);
                end
                set(lbl,'String',str); drawnow;
                obj.fitSim;
                
                obj.mc.af(:,:,j) = obj.Achi;
                obj.mc.rf(:,:,j) = obj.results;
                obj.mc.gf(:,:,j) = obj.gauparms;
                
                x = obj.aobj.x; y = obj.aobj.y; y = (y./x)./obj.opts.weightval;
                ll = (x > min(1e7./obj.lims) & x < max(1e7./obj.lims)); %indices within limits
                absy = y(ll); afit = obj.aobj.thefit./obj.opts.weightval;
                
                x = [obj.sobj.x]; y = [obj.sobj.y]; y = y./x; 
                ll = (x(:,1) > min(1e7./obj.lims) & x(:,1) < max(1e7./obj.lims)); %indices within limits
                stky = y(ll,:); 
                sfit = [obj.sobj.thefit]; 
                obj.mc.chi2abs = sum((absy-afit).^2) ;
                obj.mc.chi2stk = sum(sum((stky-sfit).^2,2));
                chi2 = obj.mc.chi2abs+obj.mc.chi2stk;
                obj.mc.chi2(j) = chi2;
            end
            
            toc;
            
            str = sprintf(['Fitting: COMPLETE \n chi^2 (',num2str(j),'): ',num2str(chi2)]);
            set(lbl,'String',str); drawnow;
            
            obj.results = originalResults;
            obj.Achi = originalAchi;
            obj.gauparms = originalGau;
            
            obj.calcErrEstimate(lbl);
        end
    % -------------------------------------------------------------------------------
        function calcErrEstimate(obj,lbl)
            
            obj.mc.rf(:,3,:) = abs(obj.mc.rf(:,3,:));

            origRI = obj.mc.ri;
            origRF = obj.mc.rf;
            origchi2 = obj.mc.chi2;
            
            chilim = floor(log10(min(obj.mc.chi2)));
            obj.mc.ri(:,:,floor(log10(obj.mc.chi2)) > chilim) = [];
            obj.mc.rf(:,:,floor(log10(obj.mc.chi2)) > chilim) = [];
            obj.mc.chi2(floor(log10(obj.mc.chi2)) > chilim) = [];
            aTooBig = squeeze(max(max(obj.mc.af > 1e-20)));
            
            obj.mc.ri(:,:,aTooBig) = [];
            obj.mc.rf(:,:,aTooBig) = [];
           
            
            indx = [1:max(obj.bands)];
            meanTrAlpha = mean(squeeze(obj.mc.rf(:,1,:)),2)';
            stdTrAlpha = std(squeeze(obj.mc.rf(:,1,:)),0,2)';
            
            meanProjAlpha = mean(squeeze(obj.mc.rf(:,2,:)),2)';
            stdProjAlpha = std(squeeze(obj.mc.rf(:,2,:)),0,2)';
            
            meanDeltaMu = mean(squeeze(obj.mc.rf(:,3,:)),2)';
            stdDeltaMu = std(squeeze(obj.mc.rf(:,3,:)),0,2)';
            
            meanZeta = mean(squeeze(obj.mc.rf(:,4,:)),2)';
            stdZeta = std(squeeze(obj.mc.rf(:,4,:)),0,2)';
            
            str1 = sprintf('Tr(Da)[0%i]: %4.2g +/- %4.2g\n',[indx;meanTrAlpha;stdTrAlpha]); 
            str2 = sprintf('m(Da)m[0%i]: %4.2g +/- %4.2g\n',[indx;meanProjAlpha;stdProjAlpha]); 
            str3 = sprintf('Delta mu[0%i]: %4.2g +/- %4.2g\n',[indx;meanDeltaMu;stdDeltaMu]); 
            str4 = sprintf('zeta[0%i]: %4.2g +/- %4.2g\n',[indx;meanZeta;stdZeta]); 
            
            strheader = sprintf('=== ERROR ESTIMATES ========== \n');
            strsep = sprintf('------------------------------ \n');
            strfooter = sprintf('============================== \n');
            
            strUsed = sprintf('%i iterations rejected\n',...
                [obj.opts.MCiter-size(obj.mc.rf,3)]);
            
            str = [strheader,str1,strsep,str2,strsep,str3,strsep,str4,strsep,strUsed,strfooter];
            %disp(str);
            set(lbl,'String',str);
            
            obj.mc.ri = origRI;
            obj.mc.rf = origRF;
            obj.mc.chi2 = origchi2;
            
            obj.mc.results = str;
            
        end
    % -------------------------------------------------------------------------------
        function plotMCresults(obj)
            figure;
            subplot(221); set(gca,'nextplot','add'); plot(squeeze(obj.mc.ri(:,1,:))',squeeze(obj.mc.rf(:,1,:))','o');
            title('Tr \Delta\alpha'); ylabel('Final');
            subplot(222); set(gca,'nextplot','add'); plot(squeeze(obj.mc.ri(:,2,:))',squeeze(obj.mc.rf(:,2,:))','o');
            title('m \cdot \Delta\alpha \cdot m');
            subplot(223); set(gca,'nextplot','add'); plot(squeeze(obj.mc.ri(:,3,:))',squeeze(obj.mc.rf(:,3,:))','o');
            title('\Delta\mu'); ylabel('Final'); xlabel('Initial');
            subplot(224); set(gca,'nextplot','add'); plot(squeeze(obj.mc.ri(:,4,:))',squeeze(obj.mc.rf(:,4,:))','o');
            title('\zeta'); xlabel('Initial'); 
        end
    % -------------------------------------------------------------------------------
        function export(obj,s)
        % thefit is in extinction, not energy-weighted extinction
        %don't do anything if no filename is given
            if strcmp(s,'') == 1, return; end
        %EXPORT ABSORPTION FIT
            ll = obj.setRange(obj.aobj);
            x = obj.aobj.x(ll); y = obj.aobj.y(ll);
            f = obj.aobj.thefit.*x; 
            g = obj.aobj.gau.*repmat(x,1,size(obj.aobj.gau,2));
            b = obj.aobj.bands;
            for j = 1:max(b)
                gg(:,j) = sum(g(:,b == j),2);
            end
            dat = [x,y,f,gg,g];
            ii = [];
            for j = 1:max(b)
                ii = [ii,1:sum(b == j)];
            end
            ss = ['wavenum \t abs \t fit \t',sprintf('band%i\t',unique(b))];
            ss = [ss,sprintf('b%igau%i\t',[b';ii]),'\n'];
            fclose('all');
            fid = fopen([s,'_ABSFIT.txt'],'w');
            fprintf(fid,ss);
            fprintf(fid,['%e',repmat('\t%e',1,size(dat,2)-1),'\n'],dat');
            fclose(fid);
        %EXPORT STARK FITS
            for j = 1:length(obj.sobj)
                ll = obj.setRange(obj.sobj(j));
                x = obj.sobj(j).x(ll); y = obj.sobj(j).y(ll);
                f = obj.sobj(j).thefit.*x;
                g = obj.sobj(j).g.*repmat(x,1,size(obj.sobj(j).g,2));
                dg = obj.sobj(j).dg.*repmat(x,1,size(obj.sobj(j).dg,2));
                ddg = obj.sobj(j).ddg.*repmat(x,1,size(obj.sobj(j).ddg,2));
                dat = [x,y,f,g,dg,ddg];
                ss = 'wavenum\tstk\tfit'; sg = ''; sdg = ''; sddg = '';
                for k = 1:max(obj.sobj(j).bands)
                    sg = [sg,'\tg',num2str(k)];
                    sdg = [sdg,'\tdg',num2str(k)];
                    sddg = [sddg,'\tddg',num2str(k)];
                end
                fid = fopen([s,'_STKFIT',num2str(j),'.txt'],'w');
                fprintf(fid,[ss,sg,sdg,sddg,'\n']);
                fprintf(fid,['%e',repmat('\t%e',1,size(dat,2)-1),'\n'],dat');
                fclose(fid);
            end
        end
    end
    methods (Static)
        function f = absfit(lam,x,obj)
            bands = obj.bands;
            obj.aobj.setGauParms(lam); 
            obj.aobj.setBands(bands);
            [g,~,~] = stkFitObj.calcGau(x,obj.aobj);
            for j = 1:max(bands)
                gg(:,j) = sum(g(:,bands == j),2);
            end
            f = sum(g,2);
            obj.aobj.setFit(f,g,gg); %fit, gaussians, progressions
            if rem(obj.kt,150) == 1, obj.plotData; end;
            obj.kt = obj.kt + 1;
        end
    % -------------------------------------------------------------------------------
        function f = simfit(lam,x,obj,lamI)
            % simfit is called as the model used to compute the Stark data in lsqcurvefit
            lam = lam.*lamI; %all parameters
            split1 = numel(obj.gauparms);
            split2 = numel(obj.results);
            alam = reshape(lam(1:split1),[],3); %absorption parameters
            slam = reshape(lam(split1+1:split1+split2),[],4); %stark parameters
            Achi = reshape(lam(split1+split2+1:length(lam)),max(obj.bands),[]);
            tralpha = slam(:,1);
            projalpha = slam(:,2);
            deltamu = slam(:,3);
            zeta = slam(:,4);
        %fit absorption
            bands = obj.bands;
            obj.aobj.setGauParms(alam); 
            obj.aobj.setBands(bands);
            [g,~,~] = stkFitObj.calcGau(obj.aobj.x(obj.setRange(obj.aobj)),obj.aobj);
            for j = 1:max(bands)
                gg(:,j) = sum(g(:,bands == j),2);
            end
            f = sum(g,2);
            obj.aobj.setFit(f,g,gg);
        %iteratively fit stark spectra
            for j = 1:length(obj.sobj)
                obj.sobj(j).setGauParms(alam);
                obj.sobj(j).setBands(bands);
                chi = obj.sobj(j).chi;
                E = 1e8^2; %(V/m)^2, normalization set to 1e6 V/cm
                % Why are all scans set to the same field?
            %convert electrooptical parameters to B & C coeffs
                [Bchi,Cchi] = stkFitObj.calcABC(tralpha,projalpha,deltamu,zeta,chi,E);
                [g,dg,ddg] = stkFitObj.calcGau(obj.sobj(j).x(obj.setRange(obj.sobj(j))),obj.sobj(j));
            %Calculate spectrum
                AchiJ = repmat(Achi(bands,j)'.*E,size(g,1),1);
                Bchi = repmat(Bchi(bands)',size(dg,1),1);
                Cchi = repmat(Cchi(bands)',size(ddg,1),1);
                for k = 1:max(bands) %for each band
                    d0(:,k) = sum(AchiJ(:,bands == k).*g(:,bands == k),2);
                    d1(:,k) = sum(Bchi(:,bands == k).*dg(:,bands == k),2);
                    d2(:,k) = sum(Cchi(:,bands == k).*ddg(:,bands == k),2);
                end
            %store derivatives
                obj.sobj(j).setDerivatives(d0,d1,d2);
                f = sum((AchiJ.*g+Bchi.*dg+Cchi.*ddg),2);
                obj.sobj(j).setFit(f,[],[]);
            end
        %output the fitted spectra and update some plots
            f = [obj.aobj.thefit./obj.opts.weightval;reshape([obj.sobj.thefit],[],1)];
            if rem(obj.kt,150) == 1, 
                obj.plotData; 
            %plot derivatives
                x = obj.sobj(1).x(obj.setRange(obj.sobj(1)));
                y = obj.sobj(1).y(obj.setRange(obj.sobj(1)));
                d0 = obj.sobj(1).g;
                d1 = obj.sobj(1).dg;
                d2 = obj.sobj(1).ddg;
                cla(obj.ax(4),'reset');
                set(obj.ax(4),'nextplot','add','box','on');
                plot(obj.ax(4),x,y./x,'.k');
                plot(obj.ax(4),x,obj.sobj(1).thefit,'Color','r','linewidth',2);
                plot(obj.ax(4),x,d0,'c');
                plot(obj.ax(4),x,d1,'b');
                plot(obj.ax(4),x,d2,'g');
                xlim(obj.ax(4),[min(x) max(x)]);
                xlab = floor(1e7./get(obj.ax(4),'Xtick'));
                set(obj.ax(4),'XTickLabel',num2str(xlab'));
            end;
            obj.kt = obj.kt + 1;
        end
    % -------------------------------------------------------------------------------
        function f = stkfit(lam,x,obj)
            slam = reshape(lam,[],4);%.*lamI;
            Achi = obj.Achi;
            tralpha = slam(:,1);
            projalpha = slam(:,2);
            deltamu = slam(:,3);
            zeta = slam(:,4);
        %get fitted absorption
            bands = obj.bands;
            alam = obj.aobj.gparms;
        %iteratively fit stark spectra
            for j = 1:length(obj.sobj)
                obj.sobj(j).setGauParms(alam);
                obj.sobj(j).setBands(bands);
                chi = obj.sobj(j).chi;
                E = 1e8^2; %(V/m)^2, normalization set to 1e6 V/cm
                [Bchi,Cchi] = stkFitObj.calcABC(tralpha,projalpha,deltamu,zeta,chi,E);
                [g,dg,ddg] = stkFitObj.calcGau(obj.sobj(j).x(obj.setRange(obj.sobj(j))),obj.sobj(j));
                AchiJ = repmat(Achi(bands,j)'.*E,size(g,1),1);
                Bchi = repmat(Bchi(bands)',size(dg,1),1);
                Cchi = repmat(Cchi(bands)',size(ddg,1),1);
                for k = 1:max(bands)
                    d0(:,k) = sum(AchiJ(:,bands == k).*g(:,bands == k),2);
                    d1(:,k) = sum(Bchi(:,bands == k).*dg(:,bands == k),2);
                    d2(:,k) = sum(Cchi(:,bands == k).*ddg(:,bands == k),2);
                end
                f = sum((AchiJ.*g+Bchi.*dg+Cchi.*ddg),2);
                obj.sobj(j).setFit(f,[],[]);
            end
            f = reshape([obj.sobj.thefit],[],1);
            if rem(obj.kt,150) == 1, 
                obj.plotData; 
            %plot derivatives
                x = obj.sobj(1).x(obj.setRange(obj.sobj(1)));
                y = obj.sobj(1).y(obj.setRange(obj.sobj(1)));
                d0 = obj.sobj(1).g;
                d1 = obj.sobj(1).dg;
                d2 = obj.sobj(1).ddg;
                cla(obj.ax(4),'reset');
                set(obj.ax(4),'nextplot','add','box','on');
                plot(obj.ax(4),x,y./x,'.k');
                plot(obj.ax(4),x,obj.sobj(1).thefit,'Color','r','linewidth',2);
                plot(obj.ax(4),x,d0,'c');
                plot(obj.ax(4),x,d1,'b');
                plot(obj.ax(4),x,d2,'g');
                xlim(obj.ax(4),[min(x) max(x)]);
            end;
            obj.kt = obj.kt + 1;
        end
    % -------------------------------------------------------------------------------
        function [g,dg,ddg] = calcGau(xvec,obj)
        %separate matrix into respective vectors
            amp = obj.gparms(:,1); mu = obj.gparms(:,2); sig = obj.gparms(:,3); 
        %resize to appropriate matrices
            [amp,x] = meshgrid(amp,xvec);
            [mu,~] = meshgrid(mu,xvec);
            [sig,~] = meshgrid(sig,xvec);
            sig = 2.*sig.^2;
        %calculate gaussians & derivatives
            g = (amp.*exp(-(x-mu).^2./(sig)))./x;
            dg = -g.*(1./x+(x-mu).*2./(sig)); %2st derivatives
            ddg = -dg.*(2.*(x-mu)./(sig)+1./x)-g.*(2./(sig)-1./x.^2); %2nd derivative
        end
    % -------------------------------------------------------------------------------
        function [Bchi,Cchi] = calcABC(tralpha,projalpha,deltamu,zeta,chi,E)
        %convert Stark parameters to B and C coefficients
            c = 2.997924e8; % meters/sec
            h = 6.62608e-34; % J s
            debye = 3.333564e-30; % C m
            ch = c*h;
            ch2 = (ch).^2;
            toa3 = 15*ch*100/1.11265e-10*1e30; % = 2.6779e18
            
            chi2 = cos(chi.*pi./180).^2;
            zeta = cos(zeta.*pi./180).^2;
            
            Bchi = 2.5.*tralpha+(3.*chi2-1).*(1.5.*projalpha-0.5.*tralpha);
            Cchi = deltamu.^2.*(5+(3.*chi2-1).*(3.*zeta-1));
            
            btoa3 = toa3./E;
            Bchi = Bchi./btoa3;
            
            ctodebye = (E./30./ch2./1e4).*debye.^2;
            Cchi = Cchi.*ctodebye;
        end
    end
end