%--------------------------------------------------------------------------------
function [mu sigma,alpha]=emalgorithmconstrained(intensity,classification,fixtoclass,maxiter)
%--------------------------------------------------------------------------

%This function estimates a gaussian mixture model for the intensities using a modified version of the EM-algorithm. 
%Input:
%intensity: vector of intenisties to be classified (size N*1)
%classification:
%subregions:
%maxiter:
%Output:
%mu: mean intenisty for the two classes (size 1*K)
%sigma: standard deviation for the two classes (size 1*K)
%alpha:

if nargin<2
  myfailed('Too few input arguments.');
end

K=max(classification);%number of classes
N=length(intensity);
mu=zeros(K,1);
sigma=zeros(K,1);
alpha=zeros(K,1);
if nargin<4
  maxiter=300;
end

%initailise mu, sigma and alpha from classification
for j=1:K
  classindex=classification==j;
  mu(j)=mean(intensity(classindex));
  sigma(j)=std(intensity(classindex));
  alpha(j)=sum(classindex)/N;
end
%The EM algorithm
% The EM algorithm is used for classifying the intenisties to K classes 
% In the Expectation step the probabilities of
% belonging to the K classes is calculated and the intensities are
% classified to the class for which it has the highest probability.
% In the Maximization step mean and standard deviation is calculated for
% the classes. In this modified version intensities can be forced to a spcific class and mean and standard
% deviation is updated for all classes but for those defined by fixclass
classprobability=zeros(N,K);
weight=zeros(N,K);
notfixed=(fixtoclass==0);
for j=1:K
    fixedtothisclass=fixtoclass==j;
    weight(fixedtothisclass,j)=1;
end
t=0;
while t<maxiter
	t=t+1;
  %calculate probability 
  for j=1:K  
    classprobability(:,j)=alpha(j)*gaussianpdf(intensity,mu(j),sigma(j));
  end
  %calculate responsibility ie weights
  for j=1:K
    weight(notfixed,j)=classprobability(notfixed,j)./sum(classprobability(notfixed,:),2);
  end
  for j=1:K
      alpha(j)=sum(weight(:,j))/sum(weight(:));
      mu(j)=sum(weight(:,j).*intensity)/sum(weight(:,j));
      sigma(j)=sqrt(sum(weight(:,j).*(intensity-mu(j)).*(intensity-mu(j)))/sum(weight(:,j)));
  end
end