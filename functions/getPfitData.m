function [pfitData]= getPfitData(pfit_input,dir,prethresh95CI)
%check out http://bootstrap-software.org/psignifit/ and http://bootstrap-software.org/psignifit/publications/wichmann_hill_2001a.pdf

global is_new_pfit % shir - for using the new pfit

pfitData.valid_thresh=sum(~(pfit_input(1:end,2)==0 | pfit_input(1:end,2)==1))>1; %the threshold calculated from pfit is only valid if there are >1 data points that aren't 1 or 0
if pfitData.valid_thresh | ~exist('prethresh95CI'), thresh_prior='';, else thresh_prior=sprintf('-cosine %f %f',prethresh95CI(1),prethresh95CI(4)); end %if the fit will be based on one data point, use PREthresh as a 'prior'
pfit_input(pfit_input(1:end,2)==0,2)=1e-6; %AZ pfit program doesn't like pefect 0s
pfit_input(pfit_input(1:end,2)==1,2)=1-(1e-6); %AZ pfit program doesn't like perfect 1s
if is_new_pfit % shir - checking if new pfit is more stable while calculating fits with small N
    pfit_input(:,2) = round(pfit_input(:,2).*pfit_input(:,3));
    options.sigmoidName = 'norm';
    options.expType = 'equalAsymptote';
    result = psignifit(pfit_input,options);
    result_std = getStandardParameters(result);
    bias=result_std(1);
    thr=result_std(2);
else
    [pfit_output pfit_Full] = pfit(pfit_input,'plot_opt','no plot','shape','cumulative gaussian','n_intervals',1,'LAMBDA_EQUALS_GAMMA',1,'sens',0,'compute_stats','true','verbose','false','CONF',[0.025 0.159 0.841 0.975]); %95 percent CI and 1SD
    if pfit_output.stats.deviance.cpe > 0.95 % shir - check if the fit is good
        warning(sprintf('Fit is not good: deviance=%.2f, cpe=%.3f',pfit_output.stats.deviance.D, pfit_output.stats.deviance.cpe));
    end
end

if is_new_pfit % shir for new pfit
    pfitData.bias=bias;
    pfitData.bias95CI=0;
    pfitData.biasSD=0;
    pfitData.thresh=thr;
    pfitData.thresh95CI=0;
    % Curve
    xi = linspace(min(dir),max(dir),1000);
    y = result.psiHandle(xi);
    pfitData.xi=xi;
    pfitData.pfitcurve=y;
    
    pfitData.Rsquare=0; %need to complete. old code for Rsquare is at the end of this script
else
    % AZ TBD - trying to replace with the new version of psignifit
    priors.gamma   = 'Uniform(0,.1)';
    priors.lambda = 'Uniform(0,.1)';
    priors.m_or_a = 'None';
    priors.w_or_b = 'None';
    % Map_output=MapEstimate(pfit_input,priors,'nafc',1,'sigmoid','gauss','core','ab','gammaislambda')
    % bias_Map = Map_output.thetahat(1); thresh_Map = Map_output.thetahat(2);
    % BayesInference(pfit_input,priors,'nafc',1,'sigmoid','gauss','core','ab','gammaislambda')
    % BootstrapInference(pfit_input,priors,'nafc',1,'sigmoid','gauss','core','ab','gammaislambda')
    
    bias = pfit_output.params.est(1);
    %bias95CI = pfit_output.params.lims(:,1); % for some reason this sometimes gives NaN? Uses BCa and not percentile? replaced by my percentile measurement below
    bias95CI = prctile(pfit_Full.params.sim(:,1),[2.5 15.9 84.1 97.5]);
    under=prctile(pfit_Full.params.sim(:,1),0.25); %remove the most extreme 0.5 percent incase the bootstrap gave far outliers (due to the limited # of datapoints) which would bias the SD
    over=prctile(pfit_Full.params.sim(:,1),99.75);
    biasSD=std(pfit_Full.params.sim(pfit_Full.params.sim(:,1)>=under & pfit_Full.params.sim(:,1)<=over,1));
    thresh = pfit_output.params.est(2);
    thresh95CI = prctile(pfit_Full.params.sim(:,2),[2.5 15.9 84.1 97.5]); % for some reason thresh95CI = pfit_output.params.lims(:,2); sometimes gives NaN? Uses BCa and not percentile? replaced by my percentile measurement below
    psy_perf = [bias,thresh];
    %[bias_GaussFit thresh_GaussFit] = cum_gaussfit_max1(pfit_input); %Yongs function (quicker) see Z:\LabTools\Matlab\TEMPO_Analysis\ProtocolSpecific\MOOG\HeadingDiscrimination\Psychometric
    
    x = dir;
    xi = x(1) : 0.1 : x(end);
    %pfitcurve = cum_gaussfit(psy_perf , xi);
    pfitcurve = cum_gaussfit(psy_perf , xi) * (1 - pfit_output.params.est(3) - pfit_output.params.est(4)) + pfit_output.params.est(3); %AZ to also use the gamma and lambda
    pfitcurve_dir_only = cum_gaussfit(psy_perf , dir); %the model values only at the actual directions used (used to calc. the R-square)
    
    pfitData.bias=bias;
    pfitData.bias95CI=bias95CI;
    pfitData.biasSD=biasSD;
    pfitData.thresh=thresh;
    pfitData.thresh95CI=thresh95CI;
    pfitData.xi = xi;
    pfitData.pfitcurve=pfitcurve;
    
    % Goodness-of-fit - Rsquare calculation
    % shir - taking into account number of trials on each dir 06/2017
    mean_pfit_input = sum(pfit_input(:,2).*pfit_input(:,3))/sum(pfit_input(:,3));
    SStot = sum(((pfit_input(:,2)-mean_pfit_input).^2).*pfit_input(:,3));
    SSerr = sum(((pfit_input(:,2)-pfitcurve_dir_only(:)).^2).*pfit_input(:,3));
    pfitData.Rsquare = 1-SSerr/SStot;
    % something i learned today: not always SStot=SSreg+SSerr 12/2017 shir
end

