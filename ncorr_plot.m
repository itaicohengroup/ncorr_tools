function [mov, scalar, options] = ncorr_plot(varargin)
% NCORR_PLOT Display strain invariant from saved ncorr data
%   [mov, strain, options] = NCORR_PLOT loads data saved from an ncorr
%   analysis session and displays the first principal strain in the current
%   configuration (i.e. the trace of the 2D Euler-Almansi strain tensor) on
%   top of the corresponding image at each time point. The resulting
%   display is exported as a stack of RGB (truecolor, 8 bit) images,
%   resulting in a 4D array or movie "mov", that can be displayed using
%   "implay(mov)". The underlying principal strain data is stored in the 3D
%   array 'strain' (dimensions are Y, X, and time). The options used during
%   analysis are saved in the structure 'options'. 
%
%   [mov, strain, options] = NCORR_PLOT(Name,Value, ...) additionally
%   specifies optional analysis parameters using Name/Value pairs. Optional
%   Name/Value pairs include:   
%
%   'DataPath'      String containing the full path to the .mat file with
%                   the ncorr analysis output data. If not provided, the
%                   program will bring up a dialog box and ask the user to
%                   find the file on the system.
%
%   'ScalarFun'     Function to compute scalar strain quantity from the XX,
%                   XY, and YY components of the 2D strain tensor or,
%                   alternatively, from the dudx, dudy, dvdx, dvdy
%                   components of the deformation gradient field (strain vs
%                   deformation gradient determined based on the 'DataType'
%                   input). If 'DataType' is 'strain', then this function
%                   should take the strain components as three separate
%                   arguments in the order: XX, XY, YY. If 'DataType' is
%                   'gradu', then this function should take the deformation
%                   gradient fields as four separate arguments in the
%                   order: dudx, dudy, dvdx, dvdy. These inputs are
%                   matrices so the function should be appropriately
%                   vectorized (e.g. xx.*yy for element-wise
%                   multiplication, instead of xx*yy). The default function
%                   computes the first strain invariant: @(xx,xy,yy) xx+yy 
%   
%   'DataType'      String, either 'strain' or 'gradu' specifying if the
%                   ScalarFun input takes the 3 strain components, or the 4
%                   deformation gradients. Default is 'strain'
%
%   'Label'         String, specifying the label for the colorbar (i.e.
%                   name for the plotted scalar quantity). Default is ''.
%
%   'View'          String, either 'cur' or 'ref', specifying to display
%                   either the Euler-Almansi strain on the current image or
%                   the Green-Lagrange strain on the refrence image,
%                   respectively. The default is 'cur'. 
%
%   'Limits'        2-element array specifying the display/color limits for
%                   the displayed strain quantity. The default is to use
%                   the full range of the data. Input NaN to use the
%                   default.
%
%   'Alpha'         Double between 0 and 1 specifying the transparency of
%                   the strain map displayed on top of the underlying
%                   image. Default is 0.5
%
%   'Resolution'    Integer specifying the output display resolution used
%                   to create mov, in DPI. Default is 150.
%
% Note: This function is designed to work with analysis output from Ncorr,
% an open source 2D digital image correlation MATLAB program by Justin
% Blaber. http://www.ncorr.com/
%
% Written by: Lena Bartell
% Last updated: November 8, 2016
% Tested using Matlab R2016a with Image Processing Toolbox
%
% Example implementations:
% >> % Default settings calculatin gte first strain invariant
% >> [mov, strain,options] = ncorr_plot;
% >> implay(mov)
%
% >> % Calculate second strain invariant, which is the determinant in 2D
% >> [mov, strain,options] = ncorr_plot('StrainFun',@(xx,xy,yy)xx.*yy-xy.^2, 'Alpha', 0.75, 'Limits', [-0.05 0.05]);
% >> implay(mov)
%
% See also: ncorr implay

% set options
options = parse_inputs(varargin);

% import data to workspace
load(options.DataPath);

% get roi mask
mask = arrayfun(@(x){x.(['roi_' options.View '_formatted']).mask},...
    data_dic_save.strains);
mask = cat(4, mask{:});
mask = repmat(mask, [1 1 3 1]);

% time points
nT = size(mask, 4);

% load images
switch options.View
    case 'ref'
        impath = [reference_save.path reference_save.name];
        im = imread(impath);
        im = repmat(im, [1 1 nT]);
        
    case 'cur'
        impath = [current_save(1).path current_save(1).name];
        im = imread(impath);
        im(end, end, nT) = 0;
        
        for tt = 2:nT
            impath = [current_save(tt).path current_save(tt).name];          
            im(:,:,tt) = imread(impath);
        end
end

% other parameters
nX1 = size(im, 1);
nX2 = size(im, 2);

% setup image as truecolor rgb
im_RGB = zeros(nX1, nX2, 3, nT);
im_map = gray(256);
for tt = 1:nT
    im_RGB(:,:,:,tt) = ind2rgb(im(:,:,tt)+1, im_map);
end

% DIC analysis region positions, in pixels
step = data_dic_save.dispinfo.spacing;
x1 = (1 : step+1 : size(im, 1)) + floor(step/2);
x2 = (1 : step+1 : size(im, 2)) + floor(step/2);

% get data and calculate scalar quantity
switch options.DataType
    case 'strain'
        Exx = cat(3, data_dic_save.strains.(['plot_exx_' options.View '_formatted']));
        Exy = cat(3, data_dic_save.strains.(['plot_exy_' options.View '_formatted']));
        Eyy = cat(3, data_dic_save.strains.(['plot_eyy_' options.View '_formatted']));
        scalar = options.ScalarFun(Exx, Exy, Eyy);
        
    case 'gradu'
        dudx = cat(3, data_dic_save.deformgrad.(['plot_dudx_' options.View '_formatted']));
        dudy = cat(3, data_dic_save.deformgrad.(['plot_dudy_' options.View '_formatted']));
        dvdx = cat(3, data_dic_save.deformgrad.(['plot_dvdx_' options.View '_formatted']));
        dvdy = cat(3, data_dic_save.deformgrad.(['plot_dvdy_' options.View '_formatted']));
        scalar = options.ScalarFun(dudx, dudy, dvdx, dvdy);
end

if ~any(isfinite(options.Limits));
    % If no limits were provided, just use the full extent of the data
    options.Limits = [min(scalar(:)), max(scalar(:))];
end
E = mat2gray(scalar, options.Limits);

% setup strain data as a 4D, truecolor rgb image
E_map = parula(256);
E_RGB = zeros(size(E,1), size(E,2), 3, nT);
for tt = 1:nT
    E_RGB(:,:,:,tt) = ind2rgb(uint8(E(:,:,tt)*256), E_map);
end

% set points outside the ROI to NaN
E_RGB(~mask) = NaN;

% setup alpha mask also
E_alpha = ones(size(E_RGB)) * options.Alpha;
E_alpha(~mask) = 0;

% setup figure and axes and colorbar
fig = figure('color', 'w');
fig.Colormap = E_map;
ax = axes('parent', fig);
ax.NextPlot = 'add';
ax.XLim = [0 size(im, 2)];
ax.YLim = [0 size(im, 1)];
ax.XTick = [];
ax.YTick = [];
ax.DataAspectRatioMode = 'manual';
ax.CLim = options.Limits;
cb = colorbar;
cb.Label.String = options.Label;

% for each time point
mov = cell(nT, 1);
for tt = 1:nT
    
    % clear axis
    cla
    
    % show image
    image(im_RGB(:,:,:,tt), 'parent', ax)
    
    % overlay strains
    image(x2([1 end]), x1([1 end]), E_RGB(:,:,:,tt),...
        'parent', ax,...
        'AlphaData', E_alpha(:,:,1,tt))
    
    % save figure as an image array
    mov{tt} = print(fig, '-RGBImage', sprintf('-r%d', options.Resolution));
    
end

% finish up
close(fig)
mov = cat(4, mov{:});

function options = parse_inputs(inputs)

% parse required and optional input arguments
p = inputParser;
addParameter(p, 'DataPath', '', @(x)validateattributes(x, {'char'},{})) % path to ncorr saved data (.mat file)
addParameter(p, 'ScalarFun', @(xx,xy,yy)xx+yy, @(x)validateattributes(x, {'function_handle'},{})) % function to calculate scalar strain/gradu quantity
addParameter(p, 'DataType', 'strain', @(x)validateattributes(x, {'char'},{})) % specify strain vs gradu data, 'strain', or 'gradu';
addParameter(p, 'Label', '', @(x)validateattributes(x, {'char'},{})) % specify strain vs gradu data, 'strain', or 'gradu';
addParameter(p, 'View', 'cur', @(x)validateattributes(x, {'char'},{})) % 'cur' for current (Eulerian) view or 'ref' for reference (Lagrangian) view
addParameter(p, 'Limits', [-Inf Inf], @(x)validateattributes(x, {'numeric'}, {'size',[1 2], 'increasing'})) % strain display limits as a 2-element array [lower upper]
addParameter(p, 'Alpha', 0.5, @(x)validateattributes(x, {'numeric'},{'numel',1,'>=',0,'<=',1})) % alpha of strain displayed on top of image
addParameter(p, 'Resolution', 150, @(x)validateattributes(x,{'numeric'},{'integer'}) ) % resolution in DPI for printing figure to create movie
p.parse(inputs{:})

% extract options
options.DataPath = p.Results.DataPath;
options.ScalarFun = p.Results.ScalarFun;
options.DataType = validatestring(p.Results.DataType, {'strain', 'gradu'});
options.Label = p.Results.Label;
options.View = validatestring(p.Results.View, {'cur','ref'});
options.Limits = p.Results.Limits;
options.Alpha = p.Results.Alpha;
options.Resolution = p.Results.Resolution;

% if path to data file was not valid (or not specified), have the user find it
if ~exist(options.DataPath, 'file')
    [FileName,PathName] = uigetfile('*.mat', 'Select ncorr data file');
    options.DataPath = [PathName FileName];
end















