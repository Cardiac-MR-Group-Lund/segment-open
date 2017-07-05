function [classification, c]=kmean(intensity,K,scarintensity,remoteintensity)
%This function performs the k-means algorithm for the data intensity.
%Input:
%intensity: vector with pixel intensities
%K: the number of clusters the algorithm should find
%scarintensity: vector with intensities in the scar region
%remoteintensity: intensites in the remote region
%Output:
%classification: vector with the same length as intensity with
%classification (a integer from 1 to K) for each pixel.
%c: vector with K cluster centers

%Felicia Seemann, june 2013. Theisis work.


if isempty(intensity)
  disp('Intensity vector is empty!')
  return;
end

%initialization
if nargin==4
  if K==1 % if only one cluster should be found: return the mean of all data and set all pixels to class 1
    c=mean(intensity);
    classification=ones(size(intensity));
    return;
  elseif K==2 % if two clusters should be found, let the initial guess be the mean of remote and scar intesities.
    c=[mean(remoteintensity); mean(scarintensity)];
    if (isempty(remoteintensity) || isempty(scarintensity))
      c=zeros(K,1);
      x=sort(intensity);
      index=floor(length(intensity)/K);
      for i=1:K
        c(i)=x(i*index); % take the values equally spaces from each other as initial guess.
      end
    end
  elseif K==3
    c=[mean(remoteintensity); min(scarintensity); mean(scarintensity)]; % eventuellt ta en annan gissning för scar. Percentil 25?
  end
elseif nargin ==2 %no avalibale data for remote and scar intensities.
  c=zeros(K,1);
  x=sort(intensity);
  index=floor(length(intensity)/K);
  for i=1:K
    c(i)=x(i*index); % take the values equally spaces from each other as initial guess.
  end
else
  myfailed('Wrong amout of inputs in kmean');
end


%assignment step
for i=1:K
  S(:,i)=abs(intensity-c(i)); % find the distances from all pixels to all cluster centers
end

[~,classification]=min(S,[],2); %make an initial classification such that all pixels are classified to the closest cluster center

t=1; %counter for while-loop
while 1
  %update cluster centers
  for i=1:K
    c(i)=mean(intensity(classification==i));
  end
  
  for i=1:K
    S(:,i)=abs(intensity-c(i)); % find the distances from all pixels to all cluster centers
  end
  
  [d,classification]=min(S,[],2);
  
  D(t)=sum(d.^2); %calculate distortion
  
  if t>3 && abs(mean(diff(D(end-3:end))))<10^-6 %checks if we have convergence
    break;
  end
  t=t+1;
end



