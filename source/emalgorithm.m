function [mu, sigma, alpha, classification]=emalgorithm(intensity, initialclassification,tol,maxiter)
%This function estimates a gaussian mixture model for the intensities using
%a modified version of the EM-algorithm with a percentage method as
%convergence criteria.
%Input:
%intensity: vector of intenisties to be classified (size N*1)
%initialclassification: vector of original classification of intenisties to the
%classes size(N*1)
%tol: the tolerance in procent for the convergence criteria. set to 0.05% if not given
%Output:
%mu: mean intenisty for the two classes (size 1*K)
%sigma: standard deviation for the two classes (size 1*K)
%alpha: the percentage of pixels each distribution consists of
% mu_plot, sigma_plot and alpha_plot is a matrix containg the values for
% mu, sigma and alpha in each iteration.
%classification: vector of final classification of intenisties to the
%classes size(N*1)

%Felicia Seemann, june 2013. Part of master theisis.

if nargin < 2
  myfailed('Too few input arguments.');
  return;
elseif nargin < 3
  tol=0.0005;
  maxiter = inf;
elseif nargin < 4
  maxiter = inf;
end

K=max(initialclassification);%number of classes
N=length(intensity); %number of pixels

%initial allocate memory for parameters
alpha=zeros(1,K);
mu=zeros(1,K);
sigma=zeros(1,K);

for k=1:K
  mu(k)=mean(intensity(initialclassification==k));
  sigma(k)=std(intensity(initialclassification==k));
  alpha(k)=sum(initialclassification==k)/N;
end

%The EM algorithm
% The EM algorithm is used for classifying the intenisties to K classes
% In the Expectation step the probabilities of
% belonging to the K classes is calculated and the intensities are
% classified to the class for which it has the highest probability.
% In the Maximization step mean and standard deviation is calculated for
% the classes.

x=linspace(0,1,4096);
classprobability=zeros(N,K);
responsibility=zeros(N,K);
t=1;

while 1
  
  %Expectation step
  for k=1:K % all classes
    classprobability(:,k)=alpha(k)*gaussianpdf(intensity,mu(k),sigma(k));
  end
  
  [~,classification(t,:)]=max(classprobability');
  
  %Maximisation step
  for k=1:K % all classes
    responsibility(:,k)=classprobability(:,k)./sum(classprobability,2);
    alpha(k)=sum(responsibility(:,k))/N;
    mu(k)=sum(responsibility(:,k).*intensity)./(sum(responsibility(:,k)));
    sigma(k)=sqrt(sum(responsibility(:,k).*(intensity-mu(k)).^2)./(sum(responsibility(:,k))));
  end
  
   %konvergenskriterie: om mindre än 0.05% procent utav pixlarna
   %byter klass inom fem iterationer

  if t>5
    for i=1:5
      C(i)=sum(classification(end-i+1,:)~=classification(end-i,:))/N;
    end
    if C<=tol
      classification=classification(end,:);
      break;
    end
  end
  if t > maxiter
    break;
  end

  t=t+1;
  
end


