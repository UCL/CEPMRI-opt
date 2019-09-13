function Fig4_SNR(mode)
% Fig4_SNR
% Plot signal from specified pixels in TO5 gels.
%
% FIG4_SNR  
% FIG4_SNR('paper') to reproduce paper figure
% Fig4_SNR('interactive')  to manually select pixels (one off)
%
% Requires data files to be available - see file cep_doctor.m for info.
%
% David Atkinson D.Atkinson@ucl.ac.uk
%
% See also cep_doctor
%

% Copyright 2018 University College London.

% Load data.
%
% In paper mode,
%    opens figure
%    shows represenative slices
%    shows plot of signal for two gels.
% In interactive mode;
%    manually select and note pixels to plot.

if nargin< 1
    mode = 'paper' ;
end

dataset = '20Mar2018' ;


% To save the DICOMs as .mat files (run once)
%action.load  = 'loaddicom' ;
%action.savematfolder = pref_uigetdir('Fig4_SNR','matfolder') ;

% Use the following to read in the .mat files
 action.load = 'loadmat' ; 



% Feb 26 data
switch dataset
    case '26Feb2018'
        
        % From IRT1 code using 'lsqnl_mag2'  option, T1s of central gels are about
        % 1460 and 1100ms.
        % dinfoIR = datparse(dselect) ;
        % [vir, mir] = d2mat(dinfoIR,{'InversionTime','itype'},'itype',ITYPE,'op','fp') ;
        % opt = 'lsqnl_mag2' ;
        % fit2 = IRT1(vir, mir.tiVec, opt) ;
        % T1display(fit2.T1)
        %
        % T2 map from scanner shows values of about 150ms and 130ms for tubes
        % dinfo = datparse(dselect) ;
        % [vt2,mt2] = d2mat(dinfo,{'slice','itype'},'itype',7 ,'op', 'dv') ;
        % eshow(vt2)
        %
        % B0map suggests of the order 50Hz across a tube
        
        sl = 13 ; % slice to use (out of 25)
        pix1 = {123, 85 } ; % Pixel in T1 1460ms gel
        pix2 = {123 ,131 } ; % Pixel in T1 1100ms gel
        
        
        fn{1} = 'IM_0038_WIP_DCE_proset121_SENSE_1901' ; % Water Selective Excitation (single shot) 6.5/15
        dn{1} = 'WSE' ;
        fn{2} = 'IM_0034_WIP_DCESPAIRlow-high_SENSE_1701' ; % SPAIR with flip angle sweep 5.8ms/10deg
        dn{2} = 'SPAIR LH' ;
        fn{3} = 'IM_0032_WIP_DCEnoFS_SENSE_1601' ; % No fat sat (single shot)  6ms/15deg
        dn{3} = 'noFS' ;
        fn{4} = 'IM_0036_WIP_DCESPAIRnodummylow-high_SENSE_1801' ;  % 5.8ms/10deg
        dn{4} = 'No Dummy' ;
        %fn_inv = 'IM_042_WIP_DCEinv_SENSE_2101' ;
        
        fndisp = [1:4] ;
        
        dyn_noise = 1 ; % Philips Dynamic Noise feature (last dynamic
        % has no RF - just noise) set to 1 if YES, zero if NO.
        
    case '20Mar2018'
        % set values below by running first with input 'interactive'
        sl = 13 ; % Choose mid slice away from artefacts
        dyn = 5 ; % frame number of a typical dynamic (all look the same!)
        pix1 = {116 , 92} ; % row column coordinate for tube 18
        pix2 = {116 , 138} ; % row column coordinate for tube
        
        dyn_noise = 0 ; % No dynamic noise data for this scan session.
        fn{1} = 'IM_0030_DCEnoFS_FA15_SENSE_1501' ;
        dn{1} = 'No FS 5.7ms, 15\circ' ; col{1} = [1 0 0] ;
        fn{2} = 'IM_0032_DCEnoFS_FA11_SENSE_1601' ;
        dn{2} = 'FA 11' ;
        fn{3} = 'IM_0034_WIP_DCESPAIRlow-high_SENSE_1701' ;
        dn{3} = 'SPAIR 5.7ms, 10\circ' ; col{3} = [0 1 0] ;
        fn{4} = 'IM_0036_WIP_DCESPAIRnodummylow-high_SENSE_1801' ;
        dn{4} = 'SPAIR LH no dummy' ;
        fn{5} = 'IM_0038_DCE_proset121_SENSE_1901' ;
        dn{5} = 'WSE 6.9ms, 11\circ' ; col{5} = [0 0 1] ;
        fn{6} = 'IM_0040_WIP_DCEinv_SENSE_2001' ;
        dn{6} = 'inv' ;
        
        fndisp = [1 3 5 ] ; % files to display
    otherwise
        error(['Unknown dataset name'])
end

disp(['Select folder containing dowloaded dataset: ',dataset])
%folder = pref_uigetdir('Fig4_SNR','folder') ; % full filename
folder = uigetdir([],'Select folder containing downloaded dataset') ;
if isnumeric(folder) || ~exist(folder,'dir')
    warning(['Folder does not exist.'])
    return
end


switch mode
    case 'interactive'
        vd = getd(fn{1}, folder, action) ;
        
        eshow(vd(:,:,13,:),'Name','choose dynamic')
        eshow(vd(:,:,:,2),'Name','choose slice and pixels')
        
        
        disp(['Note optimal slice, representative dynamic and pixels for tubes.'])
        
    case 'paper'
        
        hf = figure('DefaultAxesFontSize',12,...
            'DefaultAxesFontWeight', 'bold', ...
            'DefaultAxesLineWidth',2, ...
            'Units','centimeters') ;
        img = [] ;
        lw = 2 ;
        
        scalef = 2.5e6 ;
        gap = 4 ;
        for ifn = 1: length(fndisp)
            vd = getd(fn{fndisp(ifn)}, folder, action) ./ scalef ;
            img = cat(2, img, vd(:,:,sl,dyn)) ;
            
            t1 = squeeze(vd(pix1{1}, pix1{2},sl,1:end-dyn_noise)) ; % exclude last dynamic if noise scan
            t2 = squeeze(vd(pix2{1}, pix2{2},sl,1:end-dyn_noise)) ;
            xt=[1:20 21+gap:40+gap] ;
            plot(xt, cat(1,t1,t2),'LineWidth',2,'Color',col{fndisp(ifn)}, 'DisplayName',dn{fndisp(ifn)})
            hold on
            
            axis([1 40+gap 0 1])
            
            grid on
        end
        xlabel({'Dynamic Frame' , 'Gel T1 1460ms (left), 1100ms (right)'})
        ax= gca ;
        ax.XTick = [1 20 21+gap 40+gap];
        ax.XTickLabel = {'1', '20', '1', '20'} ;
        
        ylabel('Signal (AU)')
        
        legend('Location', 'northwest')
        
        curr_ax = gca ;
        %oldpos = curr_ax.Position ;
        %newpos = [oldpos(1) oldpos(2) oldpos(3) oldpos(3)/length(fndisp)] ;
        colormap('gray') ;
        nc = size(colormap,1) ;
        ind_im = gray2ind(mat2gray(img,[0.2 1.8]), nc) ;
        ih = 0.32 ;
        ax_im = axes(hf, 'Units','normalized','Position', [0.06 0.14 ih*length(fndisp) ih]) ;
        ax_im.Visible = 'off' ;
        ax_im.YDir = 'reverse' ;
        hi = image(ax_im, 'CData', ind_im) ;
        axis('equal')
        [nyim, nxim, ~, ~]  = size(vd) ;
        
        for ip = 1:length(fndisp)
            text(ax_im,nxim/2+(ip-1)*nxim,nyim*0.92,dn{fndisp(ip)},'Color',col{fndisp(ip)},...
                'BackgroundColor',[0.7 0.7 0.7], 'FontSize',12, ...
                'FontWeight','Bold','HorizontalAlignment','center')
        end
        
        % gel to signal
        ta = annotation('arrow') ;
        ta.Y=[0.29 0.5] ;ta.Color=[0.7 0.7 0.2];ta.X=[0.28 0.28]; ta.LineWidth = 2;
        ta.LineStyle = ':';
        
        % gel to signal
        tb = annotation('arrow') ;
        tb.Y=[0.29 0.6] ;tb.Color=[0.7 0.7 0.2];tb.X=[0.33 0.7]; tb.LineWidth = 2;
        tb.LineStyle = ':';
        
        % fat sat artefacts
        tc = annotation('arrow') ;
        tc.X=[0.4 0.45];
        tc.Y=[0.17 0.26] ;
        tc.Color=[0.7 0.7 0.2];
        tc.LineWidth=4; tc.HeadWidth=15; tc.Color='y';
        
        
        disp('displayed')
        
        
        
    otherwise
        error('Unknown mode')
end

function vd = getd(fn, folder, action)
% Gets the volume data from DICOM or MAT files. Can save to .mat
% DICOM reading code not currently public
%
% getd(filename, folder, action)
%    action.load  must be 'loaddicom'  or 'loadmat'
%    action.savematfolder  if field is present, mat file will be saved to
%                          this folder
%
% David Atkinson  D.Atkinson@ucl.ac.uk
%

switch action.load
    case {'loaddicom'}
        ffn = fullfile(folder, fn) ;
        dinfo = datparse(ffn) ;
        [vd, md, locd] = d2mat(dinfo, {'slice','dyn'},'op','fp') ;
        vd = double(vd) ;
                
    case 'loadmat'
        ffn = fullfile(folder, fn) ;
        S = load(ffn,'vd') ;
        vd = S.vd ;
end

if isfield(action, 'savematfolder')
    ffns = fullfile(action.savematfolder, [fn '.mat']) ; 
    save(ffns,'vd')
    disp(['Saving mat-file: ',ffns])
end