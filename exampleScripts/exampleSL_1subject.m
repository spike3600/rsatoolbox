function [] = exampleSL_1subject(subID);
% this is an example script for doing RSA searchlight from the outputs of
% SPM's first-level analysis.
% you will need to change the paths, etc.
% nili.hamed@gmail.com
clear;clc
%% obtain the SL definitions
% volumetric SL
fname = fullfile(groupDir,'combinedMasks','SPM','groupMask_SPM.nii');
Vmask = spm_vol(fname);
Vmask.data = spm_read_vols(Vmask);
L = rsa.defineSearchlight({Vmask},Vmask,'sphere',[15 100]);
save(fullfile(groupDir,'combinedMasks','SPM','searchlight_100.mat'),'-struct','L');

%% run the SL
% load SL definition
fname = fullfile(groupDir,'combinedMasks','SPM','searchlight_100.mat');
L = load(fname);
subjDir= fullfile(glmDir,'firstLevelDir',subj_name{subID});
cd(subjDir);
load SPM;
% Define output images
outFiles={};
for r=1:3 % number of runs
    for k=1:66 % number of dissims
        outFiles{end+1}=sprintf('dists_run%d.nii,%d',r,k); % these are the name of the files that will be written in the subject's folder (where the SPM.mat of your subject is)
    end
end

% this part of the code extract the indices of the conditions of interest
% from the design matrix.
for r=1:3
    vec = zeros(1,size(SPM.xX.xKXs.X,2));
    u = vec;
    idx = SPM.Sess(r).col(1:12);
    vec(idx) = 1:12;
    condVecs{r} = vec;
end
rsa.runSearchlight(L,SPM.xY.VY,outFiles,@rsaFunc,'optionalParams',{SPM,condVecs});

function out=rsaFunc(Y,SPM,condVecs);
B=rsa.spm.noiseNormalizeBeta(Y,SPM);      % Get prewhitened beta weights
for r=1:3
    thisB = B(condVecs{r},:);
    RDM = pdist(real(thisB),'Euclidean');
    out = [out;RDM(:)];     % Arrange the outputs in a vector, as they are written to files
end
