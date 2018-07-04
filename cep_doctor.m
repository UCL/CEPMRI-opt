function cep_doctor
% CEP_DOCTOR Contrast Enhanced Prostate Paper doctor to check user settings
%   Checks and reports on paths, installed functions etc
%
% David Atkinson  D.Atkinson@ucl.ac.uk
%

% Copyright 2018 University College London.

needsEPGX = false ;

% Fig 1 (Ernst)
if ~exist('ssSPGR.m', 'file') 
    needsEPGX = true ;
end

% Fig 2 (SIT1vsConc)
% Same requirements as Fig 1

% Fig 3 (SNR)
% DICOM files and readable .mat files are provided for distribution 
% to enable paper figure to be regenerated. (Not all functionality of Fig3
% function is available publically).
%
% Note folder with data is selected by user during call to Fig3_SNR, hence 
% files may legitimately not be on path. Here just inform user about
% download locations 
% All versions DOI: https://doi.org/10.5281/zenodo.1305056
% v1.0.0       DOI:                 10.5281/zenodo.1305057

disp(['Fig3_SNR will ask to select a folder containing data.'])
disp([' The required .mat files are in folder ''cep_mat_files'' available from Zenodo'])
disp([' To obtain Zenodo download, copy and paste this MATLAB command:'])
disp(['    web(''https://zenodo.org/record/1305057/files/cep_paper_data.zip?download=1'',''-new'')'])


% Fig 4 (EPG)
% Needs EPG-X (already tested), plus other code:
if ~exist('flip_sweep','file')
    disp(['flip_sweep.m is required'])
end
if ~exist('sq_epg_gre','file') || ~exist('build_seq','file')
    disp(['The CEP code is required. '])
    disp([' This can be downloaded by copying and pasting this MATLAB command: '])
    disp(['   web('''',''-new'')'])

if needsEPGX == true
    disp(['The EPG-X code is required. The <a href="https://github.com/mriphysics/EPG-X">original repository</a> has been forked and the '])
    disp(['version used in this paper can be downloaded by copying and pasting this MATLAB command: '])
    disp(['   web(''https://github.com/DANAJK/EPG-X/archive/master.zip'',''-new'')'])
end

    
    
