function [ success ] = plot_gaussian(in,out)
%PLOT_RESULTS Plotting of the GPExp simulation results
%   Required inputs:
%   in : the input structure of GPExp
%   out: the output (results) structure of GPExp

% Plot all the points on the gaussian distribution + 5% confidence intervals

alpha = 0.025;          % significance level
mu = 0;               % mean
sigma = 1;             % std
range = max(3,max(abs(out.outliers))) ;       % range (in the units of the std) for plotting
cutoff1 = norminv2(alpha, mu, sigma);
cutoff2 = norminv2(1-alpha, mu, sigma);
x = [linspace(mu-range*sigma,cutoff1,20), ...
    linspace(cutoff1,cutoff2,50), ...
    linspace(cutoff2,mu+range*sigma,20)];
y = normpdf2(x, mu, sigma);
plot(x,y)

xlo = [-range x(x<=cutoff1) cutoff1];
ylo = [0 y(x<=cutoff1) 0];
patch(xlo, ylo, 'r','FaceAlpha',0.25)

xhi = [cutoff2 x(x>=cutoff2) range];
yhi = [0 y(x>=cutoff2) 0];
patch(xhi, yhi, 'r','FaceAlpha',0.25)

xlabel('Number of std devs')
ylabel('Gaussian probability distr.')

for i = 1:length(out.outliers)
    string = num2str(i);
    xx = out.outliers(i);
    yy=normpdf2(xx,mu,sigma);
    text(xx,yy,string);
end

title('Normal error (5% sign. intervals)')

success = true; 

end