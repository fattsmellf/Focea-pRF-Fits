%This example code will load one session from a bilaterally imaged
%clear-skull mouse. The evoked response (DeltaF/F) has been calculated for
%the 800x800 pixel image.

%This mouse was tested with 31 bar stimuli, the predicted responses to the
%31 stimuli from the set of Gaussians used in the paper is loaded from the
%file PredResp_WF31.

%The pRF model is then fit to the data using linear regression, this takes
%some time (~30 mins) due to the large number of Gaussians we take a number of steps to improve processing speed;
%i)     We only fit every other pixek
%ii)    The regression model is fit for all pixels simultaneously using lscov
%iii)   The data is split into 12 slices and processed in parallel on a multi-core machine.
%iv)    Only the minimum value of SS is stored and updated

%Tested in MATAB 2019b, but should work in earlier versions that support
%parallel processing.

%M.W.Self 2020

clear all

fwhm = 2*sqrt(log(2)*2); %Converst std_dev to FWHM.

%Load the evoked responses ('evoked') to the 31 bar stimuli for an example session
load('Focea_ExampleWFData')

%Load the predicted responses to the 31 conditions for each of the Gaussians
%Here a reduced set of ~100,000 Gaussians is used for speed. In the MS we
%used >300,000 Gaussians
%Loads;
%TC - the predicted responses arranged [cond x Gaussian]
%gdets - the parameters of the Gaussians arranged [azimuth, elevation,std_dev]
load('Focea_PredictedResponses')
ngauss = size(gdets,1);

%% Downsampling
%No need to fit all pixels all due to smoothing inherent
%in WF images. Here we fit every other pixel to create 400 x 400 pixel maps.
nypix = 400; %Pixel resolution in y direction
nxpix = 400; %Pixel resolution in x direction
xselect = round(linspace(1,size(evoked,2),nxpix)); %Pixels which will be fitted
yselect = round(linspace(1,size(evoked,1),nypix));
%MAke crossed vals
[xix,yix] = meshgrid(xselect,yselect);
Cp = [reshape(xix,nxpix*nypix,1),reshape(yix,nxpix*nypix,1)]; %Subscripts of pixels to be fit
idx = sub2ind([800,800],Cp(:,2),Cp(:,1)); %index of to-be-fitted pixels in the original images

%% Reshape pixel matrix for GLM [nconds x npixels]
nconds = 31;
pixmat = zeros(nconds,nypix*nxpix);
for n = 1:nconds
    buf = evoked(:,:,n);
    pixmat(n,:) = buf(idx);
end

%% Split the data into 12 slices for effciient parallel processing
nslices = 12;
for s = 1:nslices
    %Intiliaze
    slice(s).ind = ones(nypix*nxpix,1); %Best-fitting index
    slice(s).bestbeta = NaN(1,nypix*nxpix); %Best-fitting beta
    slice(s).currentbest = inf(1,nypix*nxpix); %Currentbest SSE
end

%% Run pRF fit
d = floor(linspace(1,ngauss,nslices+1)); %Slice points
tic
parfor s = 1:nslices
    st = d(s);
    ed = d(s+1)-1;
    for gz = st:1:ed
        dm = TC(:,gz); %get design matrix for this Gaussian (no constant term!)
        [b,~,mse] = lscov(dm,pixmat); %%Fit the GLM using the matrix form inv(dm'*dm)*dm' (lscov is most efficient)
        %b are the beta values for each fitted pixel
        %mse are the mean squared errors of the fit (equivalent to sse as we always have the same number of data-points)
        mse(b(1,:)<0) = Inf; %Set the mean squared error for negative betas to Inf
        %Compare mse to current best, note down improvements in the minimum
        slice(s).ind(mse<slice(s).currentbest) = gz;
        slice(s).bestbeta(:,mse<slice(s).currentbest) = b(:,mse<slice(s).currentbest);
        slice(s).currentbest(mse<slice(s).currentbest) = mse(mse<slice(s).currentbest); %Update currentbest
    end
end
toc

%% Combine slices back into full dataset
best = zeros(nslices,nypix*nxpix);
allbeta = zeros(nypix*nxpix,nslices);
allind = zeros(nslices,nypix*nxpix);
for s = 1:nslices
    best(s,:) = slice(s).currentbest; %Arrange SSE slices into rows
    allbeta(:,s) = slice(s).bestbeta; %Store beta values
    allind(s,:) = slice(s).ind; %Store Gaussian index values
end
[bestval,bestix] = min(best);%Get the minimum value of SSE across the slices
ind = NaN(nypix*nxpix,1); %Index into Gaussians for best-fitting Gaussian
bestbeta = NaN(nypix*nxpix,1); %Fit Beta-value
for i = 1:nypix*nxpix
    ind(i) = allind(bestix(i),i);
    bestbeta(i) = allbeta(i,bestix(i));
end

%% READ in parameter values from best-fitting Gaussian
rmx = gdets(ind,1); %Azimuth
rmy = gdets(ind,2); %Elevation
rsd = gdets(ind,3); %Standard Deviation
%Convert std_dev back into FWHM
rsd_fwhm = rsd.*fwhm;

%% Calculate quality of fit (Pearson's)
rcrf = NaN(1,nxpix*nypix);
for N = 1:nxpix*nypix
    dm = TC(:,ind(N)); %Design matrix from best fitting Gaussian for this pixel
    %Make the predicted response using the GLM equation Y = X*BETA
    y = dm*bestbeta(N);
    %CAlculate Pearsons
    c = corrcoef(y,pixmat(:,N));
    rcrf(N) = c(1,2);
end

%% Example fit (best-fit pixel)
[i,j] = max(reshape(CRF,numel(CRF),1));
%Get the predicted time-course for the best-fitting Gaussian
%REMEMBER ind is the index of the best-fitting timecourse for each pixel
predR = TC(:,ind(j)).*bestbeta(j);
figure,bar(1:nconds,pixmat(:,j)),hold on,scatter(1:nconds,predR,'r','filled')

%% Convert parameters to 2D maps
%Correlation
CRF = reshape(rcrf,nypix,nxpix);
%pRF size
PRF = reshape(rsd_fwhm,nypix,nxpix);
%Azimuth and elevation
AZI = reshape(rmx,nypix,nxpix);
ELE = reshape(rmy,nypix,nxpix);

%% Show Maps
figure
subplot(2,2,1),imagesc(AZI),axis equal,axis tight,colormap jet,colorbar,caxis([-60 60]),title('Azimuth')
subplot(2,2,2),imagesc(ELE),axis equal,axis tight,colormap jet,colorbar,caxis([-40 40]),title('Elevation')
subplot(2,2,3),imagesc(PRF),axis equal,axis tight,colormap jet,colorbar,caxis([30 100]),title('pRF size')
subplot(2,2,4),imagesc(CRF),axis equal,axis tight,colormap jet,colorbar,caxis([0.5 1]),title('Correlation')
