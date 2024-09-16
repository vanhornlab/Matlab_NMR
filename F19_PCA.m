%% WVH PCA implementation for analyzing 19F 1H NMR data 
% Dependends on on rbnmr.m to read in bruker data directly 
close all
clearvars
clc
%% Read data with rbnmr
foldername = 'Data_6';%foldername containing processed NMR bruker data
A = rbnmr([foldername,'/']);
z=cell2mat(A); %converts A matrix to a series of matrixes
DataNum=length(A); %finds the number of experiments
X=[z(1:DataNum).Data]; %Makes a matrix suitable for PCA
%% extract spectra for region of interest (ROI)
xmin = -123.5 ;
xmax = -127;
X_axis = A{1,1}.XAxis;
abs_diff_xmin = abs(X_axis-xmin);
abs_diff_xmax = abs(X_axis-xmax);
[q,d]= min(abs_diff_xmin);
[r,c]= min(abs_diff_xmax);

W=X(d:c,:); %Extracts subset of data
XminPPM=A{1}.XAxis(d:c,:); %Extracts the ppm range

% [coeff,score,latent,tsquared,explained,mu]=pca(X);
% %Note: Each column of coeff contains the coefficients for one PC
% %Note: The columns are in descending order of compoenent variance
%% Plot overlay of all real raw data
figure (1)
subplot(1,2,1)
set ( gca, 'xdir', 'reverse' )
hold on
for i=1:length(A)
plot(A{i}.XAxis, A{i}.Data,'LineWidth', 2)
end
%Note: The x-axis needs to be flipped in Matlab
    title('Full Spectra for analysis')
    ylabel('Intensity', 'Interpreter', 'tex', ...
      'FontName', 'Arial','Fontsize',18)
    xlabel('^1^9F (ppm)', 'Interpreter', 'tex', ...
      'FontName', 'Arial','Fontsize',18)
  xlim([(A{1,1}.Procs.ABSF2) (A{1,1}.Procs.ABSF1)]) 
%Note: The above sets the limits in ppm based on the first spectrum
hold off

subplot(1,2,2)
set ( gca, 'xdir', 'reverse' )
hold on
for i=1:length(A)
plot(XminPPM, W(:,i),'LineWidth', 2)

end
%Note: The x-axis needs to be flipped in Matlab
    title('ROI spectra for Analysis')
    ylabel('Intensity', 'Interpreter', 'tex', ...
      'FontName', 'Arial','Fontsize',18)
    xlabel('^1^9F (ppm)', 'Interpreter', 'tex', ...
      'FontName', 'Arial','Fontsize',18)
   xlim([min(XminPPM) max(XminPPM)]) 
%Note: The above sets the limits in ppm based on the first spectrum
hold off
%% Run PCA and save PC1 and PC2

% [coeff,score,latent,tsquared,explained,mu]=pca(X);% PCA in the ROI spectra
%pc(:,1) = coeff(:,1)
%pc(:,2) = coeff(:,2)
%writematrix(pc,'pc_values.xlsx')

[coeff,score,latent,tsquared,explained,mu]=pca(W); % PCA over full spectra
pc(:,1) = coeff(:,1);
pc(:,2) = coeff(:,2);
writematrix(pc,'pc_values.xlsx')

