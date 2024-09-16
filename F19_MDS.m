ultidimensional scaling Matlab script
%Depends on rbnmr.m to read in Bruker data directly 
close all
clearvars
clc
addpath . %This should help with finding the file pathway name
%% Set up information
labels={"Apo", "chemical 1", "chemical 2"};%% Set up information

%% Read data with rbnmr
foldername = 'example_foldername';
A = rbnmr([foldername,'./'])
exptnum = length(labels)
z=cell2mat(A); 
DataNum=length(A); 
X=[z(1:DataNum).Data]; %Makes a matrix suitable for MDS and PCA

%%
xmin = -123
xmax = -128
X_axis = A{1,1}.XAxis
abs_diff_xmin = abs(X_axis-xmin)
abs_diff_xmax = abs(X_axis-xmax)
[q,d]=min(abs_diff_xmin)
[r,c]= min(abs_diff_xmax)

%%
Xx=X(d:c,:); %Extracts subset of data
XminPPM=A{1}.XAxis(d:c,:); %Extracts the ppm range

Xxx=(Xx(:,1))-(Xx); %Subtracts the control from other spectra
W=Xxx(:,2:exptnum); %delete zero

clear Xxx
clear DataNum
clear z
%%
% Creating a new folder for figures
figFolder = [foldername, '_figures']; 
if ~exist(figFolder, 'dir')
    mkdir(figFolder); 
end
%% Plot Overlay of all real raw data
figure (1)
subplot(1,3,1)
set (gca, 'xdir', 'reverse')
hold on
for i=1:length(A)
plot(A{i}.XAxis, A{i}.Data,'LineWidth', 2)
end

%Note: The x-axis needs to be flipped in Matlab
    title('Full Spectra')
    ylabel('Intensity', 'Interpreter', 'tex', ...
      'FontName', 'Arial','Fontsize',18)
    xlabel(['^1^9F (ppm)'], 'Interpreter', 'tex', ...
      'FontName', 'Arial','Fontsize',18)
  xlim([(A{1,1}.Procs.ABSF2) (A{1,1}.Procs.ABSF1)]) 
  legend(labels(1,:))
%Note: The above sets the limits in ppm based on the first spectrum
hold off

subplot(1,3,2)
set ( gca, 'xdir', 'reverse' )
hold on

for i=1:length(A)
plot(XminPPM, Xx(:,i),'LineWidth', 2)
end

%Note: The x-axis needs to be flipped in Matlab
    title('Spectral Region for Analysis')
    ylabel('Intensity', 'Interpreter', 'tex', ...
      'FontName', 'Arial','Fontsize',18)
    xlabel(['^1^9F (ppm)'], 'Interpreter', 'tex', ...
      'FontName', 'Arial','Fontsize',18)
   xlim([min(XminPPM) max(XminPPM)]) 
%Note: The above sets the limits in ppm based on the first spectrum
legend(labels(1,:))
hold off

subplot(1,3,3)
set ( gca, 'xdir', 'reverse' )
hold on

for i=1:length(A)-1
plot(XminPPM, W(:,i),'LineWidth', 2)
end

%Note: The x-axis needs to be flipped in Matlab
    title('Difference Spectra for Analysis')
    ylabel('Intensity', 'Interpreter', 'tex', ...
      'FontName', 'Arial','Fontsize',18)    
    xlabel(['^1^9F (ppm)'], 'Interpreter', 'tex', ...
      'FontName', 'Arial','Fontsize',18)
   xlim([min(XminPPM) max(XminPPM)]) 
   legend(labels(1,:))
%Note: The above sets the limits in ppm based on the first spectrum
hold off
clear i
%Note: The imaginary data is located in: for i=1:length(A)
figure(1);
% ... (your plotting code here) ...
saveas(gcf, fullfile(figFolder, 'Figure1.pdf')); % Saving the figure

%Dissimilarity matrix generation for MDS preparation
Wmds1 = pdist(W.'); %This computes the Euclidean distance 
%Note: a tranpose of W required

%classic MDS
[Y1,eigvals1]=cmdscale(Wmds1);
%% MDS
% mdscale(dissimilarities,2,'criterion','metricsstress');
[Y1a,stress1]=mdscale(Wmds1,2);
stress1;
%% Plot the cMDS values 
figure (3)
subplot(1,2,1)
title('Full Spectra')
plot(Y1(:,1), Y1(:,2),'ok', 'MarkerFaceColor', 'red')
text(Y1(:,1), Y1(:,2), labels(2:end))
    ylabel('cMDS2', 'Interpreter', 'tex', ...
      'FontName', 'Arial','Fontsize',18)
    xlabel(['cMDS1'], 'Interpreter', 'tex', ...
      'FontName', 'Arial','Fontsize',18)
subplot(1,2,2)
  bar(eigvals1, 'red')
    ylabel('Eigenvalue number', 'Interpreter', 'tex', ...
      'FontName', 'Arial','Fontsize',18)
    xlabel('Eigenvalue', 'Interpreter', 'tex', ...
      'FontName', 'Arial','Fontsize',18)
figure(3);
% ... (your plotting code here) ...
saveas(gcf, fullfile(figFolder, 'Figure3.pdf')); % Saving the figure
%% Plot the cMDS values based on different pdiff functions
figure (4)
plot(Y1(:,1), Y1(:,2),'ok', 'MarkerFaceColor', 'red','MarkerSize',9)
text(Y1(:,1), Y1(:,2),labels(2:end))
    ylabel('cMDS2', 'Interpreter', 'tex', ...
      'FontName', 'Arial','Fontsize',18)
    xlabel(['cMDS1'], 'Interpreter', 'tex', ...
      'FontName', 'Arial','Fontsize',18)
    title('Euclidean')
figure(4);
% ... (your plotting code here) ...
saveas(gcf, fullfile(figFolder, 'Figure4.pdf')); % Saving the figure
%% Goodness of cMDS
figure (5)
    bar(eigvals1, 'red')
    ylabel('Eigenvalue number', 'Interpreter', 'tex', ...
      'FontName', 'Arial','Fontsize',18)
    xlabel('Eigenvalue', 'Interpreter', 'tex', ...
      'FontName', 'Arial','Fontsize',18)
    title('Euclidean')
figure(5);
% ... (your plotting code here) ...
saveas(gcf, fullfile(figFolder, 'Figure5.pdf')); % Saving the figure 
    
%% Plot the cMDS values based on different pdiff functions
figure (6)
plot(Y1a(:,1), Y1a(:,2),'ok', 'MarkerFaceColor', 'red','MarkerSize',9)
text(Y1a(:,1), Y1a(:,2), labels(2:end))
    ylabel('MDS2', 'Interpreter', 'tex', ...
      'FontName', 'Arial','Fontsize',18)
    xlabel(['MDS1'], 'Interpreter', 'tex', ...
      'FontName', 'Arial','Fontsize',18)
    title('Euclidean')   
  
saveas(gcf, fullfile(figFolder, 'Figure6.pdf')); % Saving the figure