function cep_doctor
% CEP_DOCTOR Contrast Enhanced Prostate 'doctor' to help the user 
% check settings
%   Checks and reports on paths, installed functions and data location 
%
% David Atkinson  D.Atkinson@ucl.ac.uk
%

% Copyright 2019 University College London.

% Fig 1 (Ernst)
if ~exist('ssSPGR.m', 'file') 
    needsEPGX = true ;
else
    needsEPGX = false ;
end

if needsEPGX == true
    disp(' ')
    disp(['The EPG-X code is not on your path.'])
    disp(['If you need to download the code, the <a href="https://github.com/mriphysics/EPG-X">original repository</a> has been forked and the '])
    disp(['version used for this paper can be downloaded by copying and pasting this MATLAB command: '])
    disp(['   web(''https://github.com/DANAJK/EPG-X/archive/master.zip'',''-new'')'])
end


% Fig 2 (SIT1vsConc)
% Same requirements as Fig 1


% Fig 3 (EPG)
% Needs EPG-X (already tested), plus other code:
if ~exist('qflip_sweep','file')
    disp(' ')
    disp(['qflip_sweep.m is required'])
end

if ~exist('sq_epg_gre','file') || ~exist('build_seq','file')
    disp(' ')
    disp(['The CEP code is required. '])
    disp([' This can be downloaded by copying and pasting this MATLAB command: '])
    disp(['   web('''',''-new'')'])
end


% Fig 4 (SNR)
% DICOM files and readable .mat files are provided for distribution 
% to enable paper figure to be regenerated. (Not all functionality of Fig3
% function is available publically).
%
% Note folder with data is selected by user during call to Fig3_SNR, hence 
% files may legitimately not be on path. Here just inform user about
% download locations 
% All versions DOI: https://doi.org/10.5281/zenodo.1305056
% v1.0.0       DOI:                 10.5281/zenodo.1305057  (old folder names)
% v1.0.1       DOI:                 10.5281/zenodo.3407685

disp(' ')
disp(['Fig4_SNR will ask to select a folder containing data from 20Mar2018.'])
disp('The data can be downloaded from Zenodo and then unzipped.')
disp('The zip file is called ''cep_data.zip'' and the required .mat files')
disp('are in the folder ''cep_data/cep_data_mat_files''.')
disp('To download the data held on Zenodo DOI 10.5281/zenodo.1305056 :')
disp(' Within MATLAB go to the local folder where you would like to download the data')
disp(' Copy and paste this MATLAB command:')
disp('    unzip(''https://zenodo.org/record/3407685/files/cep_data.zip?download=1'')')


    
    
