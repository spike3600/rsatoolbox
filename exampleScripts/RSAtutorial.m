% RSAtutorial_fmrib.m
clear;clc;close all
rootpath = '/Volumes/TOSHIBA/rsatoolbox-develop/';
addpath(genpath(rootpath))

%% start with making a random distance matrix
nStim = 100;
nVox = 200;
a = randn(nStim,nVox);% nStim stimuli, nVox voxels (random response matrix)
nDissims = nchoosek(nStim,2);
rdv = pdist(a,'Euclidean');% this is a dissimilarity vector. a vector containing the pairwise distance between patterns for each stimulus pair.

rdm = squareform(rdv);% this is the 96x96 distance matrix

%% display RDM
rsa.fig.showRDMs(rdm,1,1,[0 2],1,[],[],'gray')
% showRDMs(RDMs,figI,rankTransform01,clims,showColorbar, aspect, imagelabels, colourScheme)
%% give it a name
d.RDM = rdm;
d.name = 'exampleRDM';
rsa.fig.showRDMs(d,2,1)

% note the effect of rank transformation here. the effect is exagerated
% when using Euclidean distance. 
% if you get stripes, it basically corresponds to one condition having
% overal larger/smaller activation compared to all other conditions.
% correlation-distance RDMs don't suffer from this, but wait for a
% discussion on this later.

%% display the matrix with image labels 
% clear il
% % try add image labels
% for i=1:nStim
%     il.images(i).image = randn(20,20,3);
% end
% rsa.fig.showRDMs(rdm,1,1,[],0,[],il)


%% define categorical model RDM (1)
% often we'd like to define and tst hypothesis about categorical representational structures
% this is conceptually similar to a decoding analysis where you test whether two classes are decodable
% in RSA we do that by correlating the empirical RDM with a "categorical"
% model RDM
% the simplest case: stimuli are either from class 1 or class 2.
catVectors = randi(2,[1,nStim]);
[binRDM, nCatCrossingsRDM] = rsa.rdm.categoricalRDM(catVectors,3);
% you can use this code to create model RDMs for your analysis
%% define categorical model RDM (2)
% in case the stimuli can have diferent class membership definitions
clear catVectors
catVectors(:,1) = randi(2,[1,nStim]);
catVectors(:,2) = randi(3,[1,nStim]);
[binRDM, nCatCrossingsRDM] = rsa.rdm.categoricalRDM(catVectors,4);

%% simulate patterns with known underlying similarity structures
%                For example,
%                        clusterSpec = {20, {6, {3, 5},{2, 3}}, {4, 7}}
%                represents the following hierarchy:
%                                            |
%                                ------------------------- 20
%                                |                       |
%                          ------------- 6               |
%                          |           |              ------- 4
%                        ----- 3      --- 2           |||||||
%                        |||||        |||             |||||||
%                         [5]         [5]               [10]
%
nVoxels = 100;
clusterSpec = {20, {6, {3, 5},{2, 5}}, {4, 10}};
pats = rsa.sim.generateBetaPatterns(clusterSpec, nVoxels);

%% use the toolbox
userOptions = defineUserOptions_tutorial;% the mentality of having all params in one file

%% display and save the figure
clear d
d.RDM = squareform(pdist(pats));
d.name = 'simRDM_cat';
rsa.fig.showRDMs(d,1);
rsa.fig.handleCurrentFigure('simulatedRDM_example',userOptions);% this allows you to save the
% figure according to what you have specfied in defineUserOptions

%% RDM explorations: MDS plots
rsa.MDSConditions(d, userOptions,struct('titleString','simulated-pats MDS',...
    'figureNumber',6));

%% RDM explorations: dendrograms
rsa.dendrogramConditions(d, userOptions,...
struct('titleString', 'Dendrogram of the ground truth RDM', 'figureNumber', 7));

%% RDM explorations: dendrograms with labels
[temp{1:5}] = deal('cat1A');[temp{6:10}] = deal('cat1B');
[temp{11:20}] = deal('cat2');
userOptions.conditionLabels = temp;

rsa.dendrogramConditions(d, userOptions,...
struct('titleString', 'Dendrogram of the ground truth RDM', 'figureNumber', 8));

%% averaging RDMs improves dissimilarity SNR
% load some simulated patterns that share a true underlying structure
load(fullfile(userOptions.rootPath,'subjectRDMs.mat'));
% load the single-subject RDMs
avgSubjectRDM=mean(subjectRDMs,3);
rsa.fig.showRDMs(rsa.rdm.concatRDMs_unwrapped(subjectRDMs,avgSubjectRDM),2);
% show the RDMs for single subjects and the group-average
rsa.fig.handleCurrentFigure(fullfile(userOptions.rootPath,'simulatedSubjAndAverage'),userOptions);


%% 
% load the models
load(fullfile(userOptions.rootPath,'modelRDMs.mat'));
% display them
rsa.fig.showRDMs(modelRDMs,9);



% put the model RDMs in a cell array
for modelRDMI=1:numel(modelRDMs)
    modelRDMs_cell{modelRDMI}=modelRDMs(modelRDMI);
end
%% second order RDM correlation matrix
userOptions.RDMcorrelationType='Kendall_taua';

avgRDM.RDM = avgSubjectRDM;
avgRDM.name = 'subject-averaged RDM';
avgRDM.color = [0 0 0];
rsa.pairwiseCorrelateRDMs({avgRDM, modelRDMs}, userOptions, struct('figureNumber', 10,'fileName','RDMcorrelationMatrix'));

%% second-order MDS plot

rsa.MDSRDMs({avgRDM, modelRDMs}, userOptions, struct('titleString', 'MDS of different RDMs', 'figureNumber', 11,'fileName','2ndOrderMDSplot'));
% this looks vert cluttered, you might want to change the color of 
% different RDMs in the MDS plot:
avgRDM.color = [1 1 0];
modelRDMs(11).color = [1 .5 .5];
close;
rsa.MDSRDMs({avgRDM, modelRDMs}, userOptions, struct('titleString', 'MDS of different RDMs', 'figureNumber', 12,'fileName','2ndOrderMDSplot'));



%% run the inference
userOptions.RDMcorrelationType='Kendall_taua';
userOptions.RDMrelatednessTest = 'subjectRFXsignedRank';
userOptions.RDMrelatednessThreshold = 0.05;
userOptions.RDMrelatednessMultipleTesting = 'FDR';
userOptions.saveFiguresPDF = 1;
userOptions.candRDMdifferencesTest = 'subjectRFXsignedRank';
userOptions.candRDMdifferencesThreshold = 0.05;
userOptions.candRDMdifferencesMultipleTesting = 'FDR';
userOptions.plotpValues = '=';
userOptions.barsOrderedByRDMCorr=true;
userOptions.resultsPath = userOptions.rootPath;
userOptions.figureIndex = [13 14];
userOptions.figure1filename = 'compareRefRDM2candRDMs_barGraph_simulatedITasRef';
userOptions.figure2filename = 'compareRefRDM2candRDMs_pValues_simulatedITasRef';
stats_p_r=rsa.compareRefRDM2candRDMs(subjectRDMs, modelRDMs_cell, userOptions);

 %% MV noise normalisation
 % for this you need your pre-processed raw time series + design matrices
 Y = randn(1200,100);% 1200 volumes, 100 voxels
 X = rand(1200,12); % 12 regressors
 indices.row{1} = 1:1200;
 indices.col{1} = 1:12;
 [u,beta,~,~,~]=rsa.whiteBeta(Y,X,indices);
 %% compute cross-validated distances
 rdv = distanceLDC(u,partition,conditionVec,X);
 % partition gives the partition (e.g. run) number for each regressor (use 0 for irrelevant regressors)
 % conditionVec gives the condition number for each regressor (1:K, use the same number for repetitions of the same condition and 0 for irrelevant regressors)