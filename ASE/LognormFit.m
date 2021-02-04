% Make up some data. (You should use your real data in place of x.)
x = readmatrix('515-2_75x-L-ESD.csv')
% Fit the data
parmhat = lognfit(x);
% Plot comparison of the histogram of the data, and the fit
figure
hold on
% Empirical distribution
hist(x,0:5:500);
% Fitted distribution
xt = 0:5:500;
plot(xt,1000*lognpdf(xt,parmhat(1),parmhat(2)),'r')
parmhat