%--------------------------------------------------------------------------
function [mu, sigma, alpha]=gaussianmixtureestimation(intensity, initialclassification, forceinitial, fixclass)
%--------------------------------------------------------------------------

%This function estimates a gaussian mixture model for the intensities using a modified version of the EM-algorithm. 
%This version is used in t2wmarsegmentation. Autumn 2014 Jane discovered
%that this version of EM is not as correct as the one Felicia implemented
%in emalgorithm.m (or gaussianintensityestimation.m used by Jane in
%lvpeter.m). This version will however continue to be used by
%t2wmarsegmentation until cessfpmarsegmentation is implented and further evaluated
%in the data for T2-STIR. 


%Input:
%intensity: vector of intenisties to be classified (size N*1)
%initialclassification: vector of original classification of intenisties to the
%classes size(N*1)
%forceinitial: vector of zeros and ones which forces a pixel to belong to
%the initialcalssification class (size N*1)
%fixclass: for a specific class mean and std can be fixed (size 1*K)
%Output:
%mu: mean intenisty for the two classes (size 1*K)
%sigma: standard deviation for the two classes (size 1*K)

if nargin<2
  disp('Gaussianmixtureestimation requires at least 2 inputs');
end

K=max(initialclassification);%number of classes
N=length(initialclassification);

if nargin<3 || isempty(forceinitial)
  forceinitial=zeros(N,1);
end

if nargin<4 || isempty(fixclass)
  fixclass=zeros(1,K);
end
  

%error checking
for j=1:K
  if isempty(find(initialclassification==j,1))
    disp(sprintf('All classes from 1 to %d must be represented in the classification vector',K));
    return;
  end
end

%intitialise mu, sigma and alpha
%alpha is the percent of intensities belogning to each class
mu=zeros(1,K);
sigma=zeros(1,K);
alpha=zeros(1,K);
for k=1:K
  currentclass=(initialclassification==k);
  currentintensity=intensity(currentclass);
  mu(k)=mean(currentintensity);
  sigma(k)=std(currentintensity);
  alpha(k)=sum(currentclass)/length(intensity);
end

%Set weights to force initialclassification
weights=ones(N,K);
for k=1:K
  index = forceinitial & (initialclassification==k);
  weights(index,:)=0;
  weights(index,k)=1;
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
for t=1:50%should later be a while loop depending on an error calcualtion
  %Expectation step
  for k=1:K % all classes
    tempgaussianpdf=alpha(k)*gaussianpdf(x,mu(k),sigma(k));   
    classprobability(:,k)=weights(:,k).*fastremap(single(intensity),single(tempgaussianpdf));
  end

  [maxvalue,classification]=max(classprobability');

  %Maximisation step
  for k=1:K % all classes
    if not(fixclass(k))
      currentclass=classification==k;
      currentintensity=intensity(currentclass);
      mu(k)=mean(currentintensity);
      sigma(k)=std(currentintensity);
      alpha(k)=sum(currentclass)/length(intensity);
    end
  end

  %calculate error (this is not yet done since no good error measurment
  %have been found)

end
