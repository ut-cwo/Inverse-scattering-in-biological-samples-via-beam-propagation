%%
% This file demonstrates the multi-slice beam propagation (MSBP) method on
% various scattering samples. You can select the type of sample to be
% reconstructed. You can choose between field-based or amplitude-based
% reconstruction, different initial guess, option to use defocs diversity.
% This allows testing under various reconstruction scenarios. All data and
% parameters related to the reconstruction are included in the sample data
% .MAT file.
%
% NOTE: this code is only beta-tested. The user is encouraged to go 
% block-by-block in this code .M file to confirm that code works as intended.
%
% On a Windows computer, 'cntrl'+'enter' runs a specific block at a time
%
% Author: Jeongsoo Kim & Shwetadwip Chowdhury; Aug 16, 2025
% Thank you to Michael Chen and David Ren, for preliminary 
% versions of this code
%
% reference:
% S. Chowdhury, M. Chen, R. Eckert, D. Ren, F. Wu, N. Repina, and L. 
% Waller, "High-resolution 3D refractive index microscopy of multiple-
% scattering samples from intensity images," Optica 6, 1211-1219 (2019) 

clear; close all; clc;

%% Set the path to the 'Data' and 'core_function' folder containing raw measurements and required sub-functions

folderPath_data     = ''; % Specify the folder path here
folderPath_function = ''; % Specify the folder path here

addpath(folderPath_data);
addpath(genpath(folderPath_function));
init_unlocbox();

%% Loading raw measurements and reconstruction parameters for each sample from .MAT files

% Please select a sample to reconstruct (options: 40um_phantom, celegans, zebrafish_embryo, organoid)
% For example : sample_name = 'zebrafish_embryo';

sample_type = '40um_phantom';     

switch sample_type

    case '40um_phantom'
        load('40um_phantom_data.mat');

    case 'celegans'
        load('celegans_data.mat');

    case 'zebrafish_embryo'
        load('zebrafish_embryo_data.mat');

    case 'organoid'
        load('organoid_data.mat');

end

%% Setting parameters relevant to reconstruction conditions and physical object volume (some are drawn from .MAT file)

save_recon     = false;                 % TRUE: save result after reconstruction; FALSE: do not save results
use_gpu        = true;                 % TRUE: use GPU device; FALSE: use computer CPU
use_field      = false;                % TRUE: uses field-component for reconstruction; FALSE: uses amplitude-only component for reconstruction
RytovInitial   = false;                % TRUE: initialize using Rytov-approximation based ODT reconstruction result; FALSE: initialize with constant zeros
refocusing     = false;                 % TRUE: use digitally refoucsed amplitude data; FALSE: use amplitude measurement only from imaging plane

ps             = ps;                   % pixel size (x,y,z) in object space (micron)
lambda         = lambda;               % central wavelength (micron)
NA             = NA;                   % numerical aperture of imaging and detection lens
n_imm          = n_imm;                % refractive index of sample immersion media
z_plane        = z_plane;              % center plane of reconstruction volume, where 0 um is object volume center
pdar           = pdar;                 % padding size to avoid edge artifacts
SkipList       = SkipList;             % skip raw measurements corrupted by dust or noise in the imaging system

%% Select region of interest (ROI) from raw measurement for reconstruction (parameters are drawn from .MAT file)

FOV_size    = FOV_size;
x_start     = x_start;
y_start     = y_start;

rows        = x_start:x_start+FOV_size-1;
cols        = y_start:y_start+FOV_size-1;

Efield_amplitude_crop = Efield_amplitude(rows,cols,:);
Efield_phase_crop     = Efield_phase(rows,cols,:);

acqs    = Efield_amplitude_crop.*exp(1i.*Efield_phase_crop);

clear Efield_amplitude Efield_phase

%% Setting spatial and frequency axes and propagation kernels

N           = size(acqs,1)+2*pdar;  % lateral pixel dimension of a padded patch within the object
x           = ps*[-N/2:N/2-1];      % 1D padded axis in x
[xx,yy]     = meshgrid(x,x);        % 2D padded grid in x/y

dfx         = 1/(N*ps);             % Fourier spacing of padded axis
fx          = dfx*[-N/2:N/2-1];     % 1D padded axis in fx
[fxx,fyy]   = meshgrid(fx,fx);      % 2D padded grid in fx/fy

fx          = ifftshift(fx);        % FFT shifting Fourier axes
fxx         = ifftshift(fxx);       % FFT shifting Fourier axes
fyy         = ifftshift(fyy);       % FFT shifting Fourier axes

% setting propagation kernels and pupil support
prop_phs            = 1i*2*pi*sqrt((n_imm/lambda)^2-(fxx.^2+fyy.^2));
NA_crop             = (fxx.^2 + fyy.^2 > (NA/lambda)^2);

% converting into GPU arrays if user targets gpu-enabling
if use_gpu
    xx              = gpuArray(xx);
    yy              = gpuArray(yy);
    fyy             = gpuArray(fyy);
    fyy             = gpuArray(fyy);
    prop_phs        = gpuArray(prop_phs);
    acqs            = gpuArray(acqs);
end

%% Setting illumination k-vectors from .MAT file and accounting for system scan-angle orientation (k-vectors are drawn from .MAT file)

fx_in    = fx_illum_ref;
fy_in    = fy_illum_ref;

%% Calculate the obliquity factor for phase delay compensation due to angled illumination

kx_in = fx_in*lambda/NA;
ky_in = fy_in*lambda/NA;
k_in  = sqrt(kx_in.^2+ky_in.^2);

theta = asin(NA*k_in/n_imm);
OF    = 1./cos(theta);

%% initializing forward model measurements and initial guess of reconstructed object (parameters are drawn from .MAT file)

O         = O;                  % axial dimension size of reconstruction space
psz       = psz;                % pixel size (z) in reconstructed object space(micron)

% initialization of guess of reconstructed object (deltaRI, not RI), to be updated iteratively

if RytovInitial
    % initialize using Rytov-approximation based ODT reconstruction result
    reconObj     = Rytov_recon_init(Efield_amplitude_crop, Efield_phase_crop, n_imm, ps, psz, O, lambda, NA, fx_in, fy_in, pdar);
else
    % initialize with constant zeros
    reconObj     = single(zeros([N, N, O]));
end

if use_gpu
    reconObj     = gpuArray(reconObj);
end

clear Efield_amplitude_crop Efield_phase_crop
%% optimization params for iterative reconstruction (parameters are drawn from .MAT file)

maxiter         = maxiter;              % number of iterations to run optimization protocol for
step_size       = step_size;            % step size for gradient-based optimization protocol
regParam        = regParam;             % regularization parameter for 3D proxTV

plot_range      = [-0.02,0.04];         % contrast to be used to show the reconstruction at each iteration
cost            = zeros(maxiter,1);     % cost function to evaluate convergence

reconObj_prox   = reconObj;             % used for Nesterov acceleration protocol for faster convergence
t_k             = 1;                    % parameter used for Nesterov acceleration


%% initializing Figure windows to observe iterative process
close all;

% triframe cross-sectional views of the reconstructed object, as it undergoes iterative updates
figure('Name','Reconstruction result');
figNum = 1;
MSBP_progview(real(reconObj),figNum,plot_range,cost, 0)

pause(0.01);

%% Perform digital refocusing of the electric field from the imaging plane to the virtual focal plane (focal plane parameters are drawn from .MAT file)

if refocusing
    FocalPlane = FocalPlane;
else
    FocalPlane = [0];
end

Prop_dataset = generate_Prop_dataset(acqs, fx_in, fy_in, dfx, prop_phs, psz, FocalPlane, xx, yy, pdar);

clear acqs


%% Running iterative optimization of object volume. Variable 'reconObj' is the final 3D refractive-index reconstruction!

gpu_fail_num = 0;
tic;

iter = 0;

while true

    iter = iter + 1;

    pause(0.01);

    % randomly scramble angles and choose without replacement
    seq = randperm(length(fx_in));
    

    for f_idx = 1:length(FocalPlane)

        focus_z = FocalPlane(f_idx);
        Prop_dataset_temp = Prop_dataset(:,:,f_idx,:);

        for illum_angle = 1:length(fx_in)

            illum_idx = seq(illum_angle);

            % Skip corrupted measurements
            if ismember(illum_idx, SkipList)
                continue;
            end

            % Compute estimated exit field on the camera plane
            [efield_prop, efield, efield_vol, U_in_phs] = MultiSlice_Forward(...
                reconObj, psz, xx, yy, dfx, prop_phs, NA_crop, lambda, ...
                fx_in(illum_idx), fy_in(illum_idx), z_plane, pdar, ...
                use_gpu, focus_z, OF(illum_idx));

            % Update RI distribution using the gradient
            [reconObj, funcVal] = BPM_update(...
                reconObj, psz, efield, efield_vol, prop_phs, NA_crop, ...
                lambda, z_plane, step_size, pdar, use_field, U_in_phs, ...
                use_gpu, Prop_dataset_temp(:,:,:,illum_idx), efield_prop, ...
                focus_z, OF(illum_idx));

            % Accumulate cost
            cost(iter) = cost(iter) + gather(funcVal);
            fprintf('illum_angle: %d  iteration: %d\n', illum_angle, iter);
        end

        if f_idx<length(FocalPlane)
            reconObj = prox_tv3d(real(reconObj), regParam);
        end

    end

    % Check early stopping criterion based on relative cost reduction
    if iter > 1
        rel_change = abs(cost(iter-1) - cost(iter)) / cost(iter);
        if rel_change < 1e-4
            fprintf('Early stopping: relative cost change < 1e-4');
            break;
        end
    end

    % Stop if maximum number of iterations reached
    if iter >= maxiter
        fprintf('Reached maximum number of iterations: %d', maxiter);
        break;
    end

    % Applying non-negativity constraint
    val = mean(mean(real(reconObj(:,:,1))));
    reconObj = reconObj -val;
    Positive_mask = gather(reconObj)<0;
    reconObj(Positive_mask) = 0;

    % Prox operator is a memory-intensiver operator. If GPU crashes due to
    % memory requirements, use CPU instead. It will be slower but the
    % program won't crash.

    try
        reconObj_prox1 = prox_tv3d(real(reconObj), regParam);
    catch
        disp('running regularizer on CPU because GPU ran out of memory for this memory-intensive procedure');
        reconObj_prox1  = prox_tv3d(gather(real(reconObj)), regParam);
        reconObj_prox1  = gpuArray(reconObj_prox1);
        gpu_fail_num    = gpu_fail_num+1;   % counter to keep track of how many times GPU failed to regularize
    end

    if iter>1
        if cost(end) > cost(end-1)
            t_k   = 1;
            reconObj = reconObj_prox;
            continue;
        end
    end


    % Nesterov's update
    t_k1       = 0.5 * (1 + sqrt(1 + 4 * t_k^2));
    beta       = (t_k - 1)/t_k1;
    reconObj   = reconObj_prox1 + beta*(reconObj_prox1 - reconObj_prox);
    t_k        = t_k1;
    reconObj_prox = reconObj_prox1;
    fprintf('iteration: %d, error: %5.5e, elapsed time: %5.2f seconds\n',iter, cost(iter));


    MSBP_progview(real(reconObj), figNum, plot_range, cost, iter)
    pause(0.01);

end

toc;
close_unlocbox();

%% In case you want to save reconstruted data and relevant parameters

if save_recon
    disp('Reconstruction is done. Saving data now');
    save('reconResult.mat','reconObj','cost', 'regParam','pdar','SkipList','-v7.3');
    disp('done saving file');
end
