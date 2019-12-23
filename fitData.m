% OUTPUTS
%     info       - data structure with following fields:
%                     .nobs      - number of observations
%                     .nb        - number of parameters
%                     .optimizer - function minimizer used
%                     .exitflag  - see FMINSEARCH
%                     .b         - fitted parameters; note that for all the
%                                  available models, the first element of B
%                                  is a noise term for the logistic
%                                  function, the remaining elements are
%                                  parameters for the selected discount
%                                  functions. eg., for model='exp', B(2) is
%                                  the time constant of the exponential
%                                  decay.
%                     .LL        - log-likelihood evaluated at maximum
%                     .LL0       - restricted (minimal model) log-likelihood
%                     .AIC       - Akaike's Information Criterion
%                     .BIC       - Schwartz's Bayesian Information Criterion
%                     .r2        - pseudo r-squared
%                   This is a struct array if multiple models are fit.
%     p           - Estimated choice probabilities evaluated at the values
%                   delays specified by the inputs vS, vR, dS, dL. This is
%                   a cell array if multiple models are fit.
%
%     EXAMPLES
%     see TEST_FAKE_DATA_AMBIGUITTY, TEST_FAKE_DATA, TEST_JOE_DATA, and TEST_KENWAY_DATA

function [info,p] = fitData(choseLottery,refVal,lottoVal,refProb,lottoProb,b0);
thresh = 0.05;
nobs = length(choseLottery);

% Fit model, attempting to use FMINUNC first, then falling back to FMINSEARCH
if exist('fminunc','file')
   try
      optimizer = 'fminunc';
      OPTIONS = optimset('Display','off','LargeScale','off','TolFun',1e-6,'TolX',1e-10);
          [b,negLL,exitflag,convg,g,H] = fminunc(@local_negLL,b0,OPTIONS,choseLottery,refVal,lottoVal,refProb,lottoProb);
          exitflagR=exitflag;
      if exitflag ~= 1 % trap occasional linesearch failures
         optimizer = 'fminsearch';
         fprintf('FMINUNC failed to converge, switching to FMINSEARCH\n');
      end
   catch
      optimizer = 'fminsearch';
      fprintf('Problem using FMINUNC, switching to FMINSEARCH\n');
   end
else
   optimizer = 'fminsearch';
end

if strcmp(optimizer,'fminsearch')
   optimizer = 'fminsearch';
   OPTIONS = optimset('Display','off','TolCon',1e-6,'TolFun',1e-5,'TolX',1e-5,...
      'DiffMinChange',1e-4,'Maxiter',100000,'MaxFunEvals',20000);
   [b,negLL,exitflag,convg] = fminsearch(@local_negLL,b0,OPTIONS,choseLottery,refVal,lottoVal,refProb,lottoProb);
end

if exitflag ~= 1
   fprintf('Optimization FAILED, #iterations = %g\n',convg.iterations);
else
   fprintf('Optimization CONVERGED, #iterations = %g\n',convg.iterations);
end

uRef = (refVal.^b(1)).*refProb;
uLotto = (lottoVal.^b(1)).*lottoProb;
slope = b(2);
p = 1 ./ (1 + exp(slope*(uLotto-uRef)));

% Unrestricted log-likelihood
LL = -negLL;
% Restricted log-likelihood
LL0 = sum((choseLottery==1).*log(0.5) + (1 - (choseLottery==1)).*log(0.5));

% Confidence interval, requires Hessian from FMINUNC
try
    invH = inv(-H);
    se = sqrt(diag(-invH));
catch
end

info.nobs = nobs;
info.nb = length(b);
info.optimizer = optimizer;
info.exitflag = exitflag;
info.b = b;

try
    info.se = se;
    info.ci = [b'-se*norminv(1-thresh/2) b'+se*norminv(1-thresh/2)]; % Wald confidence
    info.tstat = b'./se;
catch
end

info.LL = LL;
info.LL0 = LL0;
info.AIC = -2*LL + 2*length(b);
info.BIC = -2*LL + length(b)*log(nobs);
info.r2 = 1 - LL/LL0;

%----- LOCAL FUNCTIONS
function sumerr = local_negLL(beta,choseReference,refVal,lottoVal,refProb,lottoProb);

uRef = (refVal.^beta(1)).*refProb;
uLotto = (lottoVal.^beta(1)).*lottoProb;
slope = beta(2);
p = 1 ./ (1 + exp(slope*(uLotto-uRef)));

% Trap log(0)
ind = p == 1;
p(ind) = 0.9999;
ind = p == 0;
p(ind) = 0.0001;
% Log-likelihood
err = (choseReference==1).*log(p) + (1 - (choseReference==1)).*log(1-p);
% Sum of -log-likelihood
sumerr = -sum(err);
