%openExample('stats_featured/nonparametricCDFdemo')

%% Nonparametric Estimates of Cumulative Distribution Functions and Their Inverses
% This example shows how to estimate the cumulative distribution function
% (CDF) from data in a non-parametric or semi-parametric fashion. It also
% illustrates the inversion method for generating random numbers from the
% estimated CDF. The Statistics and Machine Learning Toolbox(TM) includes
% more than two dozen random number generator functions for parametric
% univariate probability distributions. These functions allow you to
% generate random inputs for a wide variety of simulations, however, there
% are situations where it is necessary to generate random values to
% simulate data that are not described by a simple parametric family.
%
% The toolbox also includes the functions |pearsrnd| and |johnsrnd|, for
% generating random values without having to specify a parametric distribution
% from which to draw--those functions allow you to specify a distribution in
% terms of its moments or quantiles, respectively.
%
% However, there are still situations where even more flexibility is needed,
% to generate random values that "imitate" data that you have collected even
% more closely. In this case, you might use a nonparametric estimate of the
% CDF of those data, and use the inversion method to generate random values.  
% The inversion method involves generating uniform random values on the unit 
% interval, and transforming them to a desired distribution using the inverse 
% CDF for that distribution.
%
% From the opposite perspective, it is sometimes desirable to use a
% nonparametric estimate of the CDF to transform observed data onto the unit
% interval, giving them an approximate uniform distribution.
%
% The |ecdf| function computes one type of nonparametric CDF estimate, the
% empirical CDF, which is a stairstep function.  This example illustrates some
% smoother alternatives, which may be more suitable for simulating or
% transforming data from a continuous distribution.

%   Copyright 2007-2014 The MathWorks, Inc.


%%
% For the purpose of illustration, here are some simple simulated data. There
% are only 25 observations, a small number chosen to make the plots in the
% example easier to read.  The data are also sorted to simplify plotting.

rng(19,'twister');
n = 25;
x = evrnd(3,1,n,1); x = sort(x);
hist(x,-2:.5:4.5);
xlabel('x'); ylabel('Frequency');

%%
% empirical cdf
[Fi,xi] = ecdf(x);
%% Kernel Estimators for the CDF and Inverse CDF
% Instead of estimating the CDF using a piecewise linear function, you can
% perform kernel estimation using the |ksdensity| function to make a smooth
% nonparametric estimate. Though it is often used to make a nonparametric
% _density_ estimate, |ksdensity| can also estimate other functions.  For
% example, to transform your original data to the unit interval, use it to
% estimate the CDF.
F = ksdensity(x, x, 'function','cdf', 'width',.35);
stairs(xi,Fi,'r');
hold on
plot(x,F,'b.');
hold off
xlim([-2.5 5]); xlabel('x'); ylabel('F(x)');
legend({'ECDF' 'Kernel Estimates'},'location','NW');

%%
% |ksdensity| also provides a convenient way to evaluate the kernel CDF
% estimate at points other than the original data.  For example, plot the
% estimate as a smooth curve.
y = linspace(-2.5,5,1000);
Fy = ksdensity(x, y, 'function','cdf', 'width',.35);
stairs(xi,Fi,'r');
hold on
plot(y,Fy,'b-');
hold off
legend({'ECDF' 'Kernel Estimate'},'location','NW');

%%
% |ksdensity| uses a bandwidth parameter to control the amount of smoothing in
% the estimates it computes, and it is possible to let |ksdensity| choose a
% default value.  The examples here use a fairly small bandwidth to limit the
% amount of smoothing.  Even so, the kernel estimate does not follow the ECDF
% as closely as the piecewise linear estimate does. 

%%
% One way to estimate the inverse CDF using kernel estimation is to compute the
% kernel CDF estimate on a grid of points spanning the range of the original
% data, and then use the same procedure as for the piecewise linear estimate.  For
% example, to plot the inverse CDF kernel estimate as a smooth curve, simply
% swap the axes.
stairs(Fi,[xi(2:end); xi(end)],'r');
hold on
plot(Fy,y,'b-');
hold off
ylim([-2.5 5]); ylabel('x'); xlabel('F(x)');
legend({'ECDF' 'Kernel Estimate'},'location','NW');

%% 
% To transform uniform random values back to the scale of the original data, 
% interpolate using the grid of CDF estimates.
u = rand(10000,1);
Finv = @(u) interp1(Fy,y,u,'linear','extrap');
hist(Finv(u),-2.5:.25:5);
xlabel('x'); ylabel('Frequency');

%%
% Notice that the simulated data using the kernel CDF estimate has not
% completely "smoothed over" the two individual observations to the left of
% zero present in the original data.  The kernel estimate uses a fixed
% bandwidth.  With the particular bandwidth value used in this example, those
% two observations contribute to two localized areas of density, rather than a
% wide flat region as was the case with the piecewise linear estimate. In
% contrast, the kernel estimate has smoothed the data _more_ in the right tail
% than the piecewise linear estimate.

%%%%%%%
%% Here %%
%%%%%%%
%%
% Another way to generate simulated data using kernel estimation is to use
% |ksdensity| to compute an estimate of the inverse CDF directly, again using
% the |'function'| parameter.  For example, transform those same uniform values.
r = ksdensity(x, u, 'function','icdf', 'width',.35);

%%
% However, using the latter method can be time-consuming for large amounts of
% data. A simpler, but equivalent, method is to resample with replacement from
% the original data and add some appropriate random noise.
r = datasample(x,100000,'replace',true) + normrnd(0,.35,100000,1);

%%
% If you generate enough random values, a histogram of the result follows the
% kernel density estimate of the original data very closely.
binwidth = .25;
edges = -2.5:binwidth:6;
ctrs = edges(1:end-1) + binwidth./2;
counts = histc(r,edges); counts = counts(1:end-1);
bar(ctrs,counts./(sum(counts).*binwidth),1,'FaceColor',[.9 .9 .9]);
hold on
xgrid = edges(1):.1:edges(end);
fgrid = ksdensity(x, xgrid, 'function','pdf', 'width',.3);
plot(xgrid,fgrid,'k-');
hold off
xlabel('x'); ylabel('f(x)');

%%
% This semiparametric estimate has smoothed the tails of the data more than
% the center, because of the parametric model used in the tails.  In that
% sense, the estimate is more similar to the piecewise linear estimate than to
% the kernel estimate. However, it is also possible to use |paretotails| to
% create a semiparametric fit that uses kernel estimation in the center of the
% data.


%% Conclusions
% This example illustrates three methods for computing a non- or semi-parametric
% CDF or inverse CDF estimate from data.  The three methods impose different
% amounts and types of smoothing on the data.  Which method you choose depends
% on how each method captures or fails to capture what you consider the
% important features of your data.
