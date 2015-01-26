function [fun, B] = hsifcn(type, xobs, yobs, varargin)
%HSIFCN Calculate an HSI function handle
%
% [fun, B] = hsifcn(type, xobs, yobs)
% [fun, B] = hsifcn(type, xobs, yobs, p1, v1, ...)
%
% This function calculates a habitat suitability model for a particular
% species, based on count data and contemporaneous measurement of one or
% more habitat properties.  It can use one of two algorithms: suitability
% index, or logistic regression.
%
% The cdf method uses the cumulative frequency method described in Bortone,
% 2002 (pg. 257).  Each predictor variable is fit separately, with no
% consideration for interaction terms; the final suitability value is a
% simple product of each individual predictor suitability.  The logistic
% regression method uses a generalized linear model with a logit link
% function, and does consider all pairwise interaction terms between the
% various predictors.
%
% Input variables:
%
%   type:   'cdf':          Suitability index based on cumulative
%                           distribution function of presence/absence
%           'logistic':     Logistic regression of presence/absence
%           'cdfcount':     Suitability index based on cumulative
%                           distribution function of #fish 
%           'poisson':      Poisson regression of #fish
%
%   xobs:   nobs x nh array, habitat properties for each observation
%
%   yobs:   nobs x 1, number of critters observed at each observation
%
% Optional input variables:
%
%   P:      SLM prescription structure (see slmset.m) (suitability only).
%           By default, performs a piecewise cubic with 20 equally-spaced
%           knots, constrained to be monotonically increasing with a
%           minimum value of 0 and a slope of 0 at the rightmost knot.
%
%   link:   link function (logistic only) used for the GLM. ['logit']
%
%   model:  starting model for GLM stepwise process (logistic only)
%           ['quadratic']
%
%   glmpv:  cell array of parameter/value pairs for GLM model creation
%
% Output variables:
%
%   fun:    function handle, of the form fun = f(x1, x2, ...), where the
%           input variables are equally-sized arrays (of any size)
%           corresponding to the habitat properties in xobs.  The
%           logistic-regression model function returns 2 outputs (fit and
%           confidence intervals); the cdf method returns only 1 (fit).
%
%   B:      1 x 1 structure holding additional diagnostic variables
%           For suitability:
%           s:      1 x nh structure array, best-fit SLM model structures
%                   produced by slmengine for each predictor variable
%           For logistic regression:
%           mdl:    generalized linear model object

% Copyright 2013 Kelly Kearney

%----------------------
% Check input
%----------------------

if ndims(xobs) ~= 2
    error('x must be nobs x nhab array');
end

[nobs, nx] = size(xobs);

if ~isequal(size(yobs), [nobs 1])
    error('y must be nobs x 1 vector');
end

% Parse optional input variables

Opt.P = slmset('increasing', 'on', ...
               'degree', 3, ...
               'minvalue', 0, ...
               'knots', 20, ...
               'interiorknots', 'fixed', ...
               'rightslope', 0);
if strcmp(type, 'logistic')
    Opt.link = 'logit';
else
    Opt.link = 'log';
end
Opt.model = 'quadratic';
Opt.glmpv = cell(0);
Opt = parsepv(Opt, varargin);

%----------------------
% Calculate HSI
%----------------------

switch type
    
    case 'cdf'
        
        for ix = 1:nx
            
            % Calculate histograms of total observations and presence-only
            % observations
            
            x = xobs(yobs>0,ix);
            
            xedge = linspace(min(xobs(:,ix)), max(xobs(:,ix)), 100);
            nall = histc(xobs(:,ix), xedge);
            npresent = histc(x, xedge);
            
            % Final histogram is frequency-present
            
            nfreq = npresent./nall;
            isabs = nall == 0;
            nfreq(isabs) = 0; % Not sure this is the best way, probably should eliminate
            freq = cumsum(nfreq)./sum(nfreq);
            
            % Calculate spline fit to cumulative frequency histogram
            
            s(ix) = slmengine(xedge', freq, Opt.P);
        
        end
        
        % HSI function
        
        fun = @(varargin) suitcalc(s, varargin{:});
        
        B.s = s;
        
    case 'cdfcount'
        
        for ix = 1:nx
            
            % Coyne&Christensen method: no binning.  Not sure I'm
            % normalizing to 1 in the same manner; their description was
            % software-specific
            
            [xunq, yavg] = consolidator(xobs(:,ix), yobs, @mean);
            freq = cumsum(yavg)./sum(yavg);
            
            % Calculate spline fit to cumulative frequency histogram
            % (Coyne&Christensen use piecewise linear fit, but I allow any
            % slm prescription)
            
            s(ix) = slmengine(xunq', freq, Opt.P);

        end
        
        % HSI function
        
        fun = @(varargin) suitcalc(s, varargin{:});
        
        B.s = s;
        
    case 'cdfjoint'
        
       
        
        
    case 'logistic'
        
        % Fit logistic regression, starting with linear, quadratic, and
        % interaction terms; use stepwise process to throw out
        % insignificant terms 
        
        [wmsg, wid] = lastwarn;
        lastwarn('');
        
        B.mdl = GeneralizedLinearModel.stepwise(xobs, yobs>0, Opt.model, 'distribution', 'binomial', 'link', Opt.link, Opt.glmpv{:});
        
        [msg, id] = lastwarn;
        if strcmp(id, 'stats:LinearModel:RankDefDesignMat')
            % This is a hack fix... throw out the troublesome singular
            % variable.  Not sure if I'm screwing up the math here
            singvar =  B.mdl.CoefficientNames{B.mdl.Coefficients.Estimate == 0};
            B.mdl = removeTerms(B.mdl, singvar);
        end
        lastwarn(wmsg,wid);
        
        % HSI function
        
        fun = @(varargin) glmcalc(B.mdl, varargin{:});
        
        
    case 'poisson'
        
        B.mdl = GeneralizedLinearModel.stepwise(xobs, yobs, Opt.model, 'distribution', 'poisson', 'link', Opt.link, Opt.glmpv{:});

        % HSI function
        
        fun = @(varargin) glmcalc(B.mdl, varargin{:});
        
%         xint = zeros(nobs,0);
%         for ix = 1:nx
%             for iint = (ix+1):nx
%                 xint = [xint xobs(:,ix).*xobs(:,iint)];
%             end
%         end
%         x = [xobs xobs.^2 xint];
%         
%         [b, dev, Stats] = glmfit(x, yobs, 'poisson', 'link', Opt.link);
% 
%         fun = @(varargin) lrcalc(b, Opt.link, varargin{:});
%         B.b = b;
%         B.dev = dev;
%         B.Stats = Stats;
    otherwise
        error('Not a valid HSI type');
end


%----------------------
% Suitability-product 
% function
%----------------------

% Function to calculate product of all independent suitability index
% values (flexible to # of habitat variables)

function hsi = suitcalc(s, varargin) 
hsi = cell(size(s));
for is = 1:length(s)
    hsi{is} = slmeval(varargin{is}, s(is), 1)./slmpar(s(is), 'maxslope');
end
ndim = ndims(varargin{1});
hsi = cat(ndim+1, hsi{:});
hsi = prod(hsi, ndim+1);

%----------------------
% Logistic regression
% function
%----------------------

function hsi = lrcalc(b, link, p, varargin)

% Expand input to include squared and interaction terms

sz = size(varargin{1});
xfit = cellfun(@(a) a(:), varargin, 'uni', 0);
xfit = cat(2, xfit{:});

xint = zeros(size(xfit,1),0);
nx = length(varargin);

for ix = 1:nx
    for iint = (ix+1):nx
        xint = [xint xfit(:,ix).*xfit(iint)];
    end
end
xfit = [xfit xfit.^2 xint];

% Check which terms are significant, and drop any that aren't 

issig = p < 0.05;
issig(1) = true; % Always keep the constant term

% Evaluate and reshape

hsi = glmval(b(issig), xfit(:,issig(2:end)), link);
hsi = reshape(hsi, sz);

%----------------------
% GLM model calc
%----------------------

function [hsi, hsici] = glmcalc(mdl, varargin)

ciflag = false; 
if ischar(varargin{end}) && strcmp(varargin{end}, 'noci')
    ciflag = true;
    varargin = varargin(1:end-1);
end

sz = size(varargin{1});
xfit = cellfun(@(a) a(:), varargin, 'uni', 0);
xfit = cat(2, xfit{:});

if ciflag
    hsi = predict(mdl, xfit);
    hsi = reshape(hsi, sz);
    hsici = nan([size(hsi) 2]);
else
    [hsi, citmp] = predict(mdl, xfit);
    hsi = reshape(hsi, sz);
    hsici = zeros([size(hsi) 2]);
    hsici1 = reshape(citmp(:,1), sz);
    hsici2 = reshape(citmp(:,2), sz);
    hsici = cat(length(sz)+1, hsici1, hsici2);
end


