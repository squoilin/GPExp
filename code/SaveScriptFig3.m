figure('Visible','off')
alpha = 0.025;          % significance level
mu = 0;               % mean
sigma = 1;             % std
range = max(3,max(abs(handles.results.outliers))) ;       % range (in the units of the std) for plotting
cutoff1 = norminv(alpha, mu, sigma);
cutoff2 = norminv(1-alpha, mu, sigma);
x = [linspace(mu-range*sigma,cutoff1,20), ...
    linspace(cutoff1,cutoff2,50), ...
    linspace(cutoff2,mu+range*sigma,20)];
y = normpdf(x, mu, sigma);
plot(x,y);

xlo = [-range x(x<=cutoff1) cutoff1];
ylo = [0 y(x<=cutoff1) 0];
patch(xlo, ylo, 'r','FaceAlpha',0.25)

xhi = [cutoff2 x(x>=cutoff2) range];
yhi = [0 y(x>=cutoff2) 0];
patch(xhi, yhi, 'r','FaceAlpha',0.25)

xlabel('Number of standard deviations')
ylabel('Gaussian probability distribution')

for i = 1:length(handles.results.outliers)
    string = num2str(i);
    xx = handles.results.outliers(i);
    yy=normpdf(xx,mu,sigma);
    text(xx,yy,string);
end

 saveas(gcf,[handles.in.name '_fig3'],'png')
 
%  title('Detection of outliers: Normal distribution of the error, with 5% significance intervals')
 
 delete(gcf)