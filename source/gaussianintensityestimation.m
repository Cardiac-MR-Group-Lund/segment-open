%--------------------------------------------------------------------------------
function [mu sigma,alpha]=gaussianintensityestimation(intensity,mu,sigma,fixclass)
%--------------------------------------------------------------------------

%This function estimates a gaussian mixture model for the intensities using a modified version of the EM-algorithm. 
%Input:
%intensity: vector of intenisties to be classified (size N*1)
%mu1: intial estiamtion of the mean value of each class (size 1*K)
%sigma1: inital estiamtion of the standard deviation of each class (sze
%1*K)
%fixclass: for a specific class mean and std can be fixed (size 1*K)
%Output:
%mu: mean intenisty for the two classes (size 1*K)
%sigma: standard deviation for the two classes (size 1*K)

if nargin<3
  disp('Gaussianmixtureestimation requires at least 3 inputs');
end

K=length(mu);%number of classes
N=length(intensity);

if nargin<4 || isempty(fixclass)
  fixclass=zeros(1,K);
end

%The EM algorithm
% The EM algorithm is used for classifying the intenisties to K classes 
% In the Expectation step the probabilities of
% belonging to the K classes is calculated and the intensities are
% classified to the class for which it has the highest probability.
% In the Maximization step mean and standard deviation is calculated for
% the classes. In this modified version intensities can be forced to a spcific class and mean and standard
% deviation is updated for all classes but for those defined by fixclass
x=linspace(0,1,4096);
classprobability=zeros(N,K);
alpha=1/3*ones(1,K);
oldalpha=zeros(1,K);
t=0;
while sum(abs(alpha-oldalpha))>10^-6 && t<300%should later be a while loop depending on an error calcualtion
	oldalpha=alpha;
	t=t+1;
  %Expectation step
  for j=1:K  
    tempgaussianpdf=alpha(j)*gaussianpdf(x,mu(j),sigma(j));
    classprobability(:,j)=fastremap(intensity,single(tempgaussianpdf));
	end

  %Maximisation step
  for j=1:K
		responsibility=classprobability(:,j)./sum(classprobability,2);
		if not(fixclass(j))
			mu(j)=sum(responsibility.*intensity)/sum(responsibility);
			sigma(j)=sqrt(sum(responsibility.*(intensity-mu(j)).*(intensity-mu(j)))/sum(responsibility));
		end
		alpha(j)=sum(responsibility)/length(responsibility);
	end
end