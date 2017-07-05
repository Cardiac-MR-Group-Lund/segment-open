function f = gaussianpdf(x,mu,sigma)
%This function creates a gaussian distribution. The function is written in
%order to not use the built in matlab command pdf('norm',x,mu,sigma) since
%we then have to compile the statistics tool box into Segment
%
%Input:
%x: values for which the distribution should be calculated
%mu: mean
%sigma: standard deviation
%Output:
%f: gaussian distribution

%Written by Nils Lundahl
if length(mu)==1
  f=exp(-0.5 * ((x - mu)./sigma).^2) ./ (sqrt(2*pi) .* sigma);
elseif length(mu)==2
  mu1=mu(1);
  mu2=mu(2);
  sigma1=sqrt(sigma(1,1));
  sigma2=sqrt(sigma(2,2));
  ro=sigma(1,2)/(sigma1*sigma2);
  x1=x(:,1);
  x2=x(:,2);
  f=1/(2*pi*sigma1*sigma2*sqrt(1-ro*ro))*exp(-1/(2*1-ro*ro)*((x1-mu1).*(x1-mu1)/(sigma1*sigma1)+(x2-mu2).*(x2-mu2)/(sigma2*sigma2)+2*ro/(sigma1*sigma2)*(x1-mu1).*(x2-mu2)));
else
  myerror('Normal distributions can only be calculated for one or two dimensional data');
end