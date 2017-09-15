function [mov, options, scalar_data] = ncorr_plot(varargin)
% NCORR_PLOT Display scalar from saved ncorr data as a movie
%   [mov, options, scalar] = NCORR_PLOT loads data saved from an ncorr
%   analysis session and displays the first principal strain in the current
%   configuration (i.e. the trace of the 2D Euler-Almansi strain tensor) on
%   top of the corresponding image at each time point. The resulting
%   display is exported as a stack of RGB (truecolor, 8 bit) images,
%   resulting in a 4D array or movie "mov", that can be displayed using
%   "implay(mov)". The underlying principal strain data is stored in the 3D
%   array 'scalar' (dimensions are Y, X, and time). The options used during
%   analysis are saved in the structure 'options'. 
%
%   [mov, options, scalar] = NCORR_PLOT(Name,Value, ...) additionally
%   specifies optional analysis parameters using Name/Value pairs. Optional
%   Name/Value pairs include:   
%
%   'DataPath'      String containing the full path to the .mat file with
%                   the ncorr analysis output data. If not provided, the
%                   program will bring up a dialog box and ask the user to
%                   find the file on the system.
%
%   
%   'DataType'      String, 'strain', 'gradu', or 'disp', specifying if the
%                   ScalarFun input takes the 3 strain components, the 4
%                   deformation gradients, or the 2 displacements. Default
%                   is 'strain'.
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
%   'TimePoints'    1D array of time points to plot. Default is 1:nT where
%                   nT is the total number of time points. 
%
%   'SymCMap'       Single boolean. Use a center-symmertric color map.
%                   Default is false. 
%
%   'RedMap'        Single boolean. Use a red-to-white color map instead of
%                   parula. Default is false.
%
%   'Rotate90'      Single boolean. Rotate axes 90 degrees. Default is
%                   false.
%
% Note: This function is designed to work with analysis output from Ncorr,
% an open source 2D digital image correlation MATLAB program by Justin
% Blaber (http://www.ncorr.com/). To use the 'gradu' DataType option,
% however, you must use Lena Bartell's modified version of Ncorr, which
% additionally calculates and saves the displacement gradient components.
%
% Written by Lena Bartell
% October, 2016
% Matlab R2016a or 2017a with Image Processing Toolbox
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
fprintf('Loading data...'), tic
load(options.DataPath);
fprintf('done (%02f s).\n', toc)

% get time points
timepts = options.TimePoints;
if any(isnan(timepts)) % if NaN, use default: all points
    nT = length(data_dic_save.strains);
    timepts = 1:nT;    
else % use given list of time points
    nT = length(timepts);
end

% load data for first time point
im = get_image(reference_save, current_save, options, timepts(1));
[scalar, alpha_mask, ~, scalar_map] = get_scalar(data_dic_save, options, timepts(1));

% setup figure, axes and colorbar
fig = figure('color', 'w');
fig.Colormap = scalar_map;
ax = axes('parent', fig, 'fontname', 'Garamond');
ax.NextPlot = 'add';
ax.XLim = [0 size(im, 2)];
ax.YLim = [0 size(im, 1)];
ax.XTick = [];
ax.YTick = [];
ax.DataAspectRatioMode = 'manual';
cb = colorbar;
cb.Label.String = options.Label;

% DIC analysis region positions, in pixels
step = data_dic_save.dispinfo.spacing;
x1 = (1 : step+1 : size(im, 1)) + floor(step/2);
x2 = (1 : step+1 : size(im, 2)) + floor(step/2);

% setup data image hanldes
h_image = image(im, 'parent', ax);
h_scalar = image(x2([1 end]), x1([1 end]), scalar, 'parent', ax, 'AlphaData', alpha_mask);

% initialize movie
mov = cell(nT, 1);
if nargout > 2
    savescalar = true;
    scalar_data = cell(nT, 1);
else
    savescalar = false;
end

% for each time point
for ii = 1:nT
    tt = timepts(ii);
    
    % get and plot new image data, as needed
    if strcmp(options.View, 'cur')
        h_image.CData = get_image(reference_save, current_save, options, tt);
    end
    if options.Rotate90
        h_image.CData = flip(permute(h_image.CData, [2 1 3]),1);
    end
    
    % get and plot new scalar data
    [h_scalar.CData, h_scalar.AlphaData, ax.CLim, raw_scalar] = get_scalar(data_dic_save, options, tt);
    if options.Rotate90
        h_scalar.CData = flip(permute(h_scalar.CData, [2 1 3]),1);
        h_scalar.AlphaData = flip(permute(h_scalar.AlphaData, [2 1 3]),1);
        raw_scalar = flip(permute(raw_scalar, [2 1 3]),1);
    end
    if savescalar
        scalar_data{ii} = raw_scalar;
    end
    

    % save the display
    mov{ii} = print(fig, '-RGBImage', sprintf('-r%d', options.Resolution));
end

% finish up
close(fig)
mov = cat(4, mov{:});
if savescalar
    scalar_data = cat(4, scalar_data{:});
end


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
addParameter(p, 'TimePoints', NaN, @(x)validateattributes(x,{'numeric'}, {'integer', 'vector'}) ); % array of time points to plot
addParameter(p, 'SymCMap', false, @(x)validateattributes(x,{'logical'}, {'scalar'}) ); % colormap
addParameter(p, 'RedMap', false, @(x)validateattributes(x,{'logical'}, {'scalar'}) ); % colormap
addParameter(p, 'Rotate90', false, @(x)validateattributes(x,{'logical'}, {'scalar'}) ); % colormap
p.parse(inputs{:})

% extract options
options.DataPath = p.Results.DataPath;
options.ScalarFun = p.Results.ScalarFun;
options.DataType = validatestring(p.Results.DataType, {'strain', 'gradu', 'disp'});
options.Label = p.Results.Label;
options.View = validatestring(p.Results.View, {'cur','ref'});
options.Limits = p.Results.Limits;
options.Alpha = p.Results.Alpha;
options.Resolution = p.Results.Resolution;
options.TimePoints = p.Results.TimePoints;
options.SymCMap = p.Results.SymCMap;
options.RedMap = p.Results.RedMap;
options.Rotate90 = p.Results.Rotate90;

% if path to data file was not valid (or not specified), have the user find it
if ~exist(options.DataPath, 'file')
    [FileName,PathName] = uigetfile('*.mat', 'Select ncorr data file');
    options.DataPath = [PathName FileName];
end

function im = get_image(reference_save, current_save, options, tt)
% get image & return as truecolor rgb

% load image
switch options.View
    case 'ref'
        im = imread([reference_save.path reference_save.name]);
        
    case 'cur'
        im = imread([current_save(tt).path current_save(tt).name]);
end

% convert to grayscale if not already
if size(im, 3) > 1
    im = rgb2gray(im);
end

% create truecolor image
im = gray2ind(im, 2^8);
im = ind2rgb(im, gray(2^8));

function [scalar, alpha_mask, limits, scalar_map, raw_scalar] = get_scalar(data_dic_save, options, tt)

% get tensor components & calculate scalar quantity
switch options.DataType
    case 'strain'
        Exx = data_dic_save.strains(tt).(['plot_exx_' options.View '_formatted']);
        Exy = data_dic_save.strains(tt).(['plot_exy_' options.View '_formatted']);
        Eyy = data_dic_save.strains(tt).(['plot_eyy_' options.View '_formatted']);
        scalar = options.ScalarFun(Exx, Exy, Eyy);
        
    case 'gradu'
        dudx = data_dic_save.deformgrad(tt).(['plot_dudx_' options.View '_formatted']);
        dudy = data_dic_save.deformgrad(tt).(['plot_dudy_' options.View '_formatted']);
        dvdx = data_dic_save.deformgrad(tt).(['plot_dvdx_' options.View '_formatted']);
        dvdy = data_dic_save.deformgrad(tt).(['plot_dvdy_' options.View '_formatted']);
        scalar = options.ScalarFun(dudx, dudy, dvdx, dvdy);
    
    case 'disp'
        u = data_dic_save.displacements(tt).(['plot_u_' options.View '_formatted']);
        v = data_dic_save.displacements(tt).(['plot_v_' options.View '_formatted']);
        scalar = options.ScalarFun(u, v);
end

% save raw scalar values
raw_scalar = scalar;

% get/set scalar limits
limits = options.Limits;
if ~any(isfinite(limits));
    % If no limits were provided, just use the full extent of the data
    limits = [min(scalar(:)), max(scalar(:))];
end

% create truecolor image of scalar output
bit_depth = 2^8;
scalar_map = parula(bit_depth);
if options.RedMap % Use CornellRed-to-white color map instead of parula
    scalar_map = cat(1,linspace(175,255,bit_depth),linspace(27,255,bit_depth),...
        linspace(27,255,bit_depth))'/255; % cornell red: [175 27 27]
end
if options.SymCMap % Make the colormap symmetric
    m = scalar_map(1:2:end,:);
    scalar_map = cat(1, flip(m, 1), m);
end

scalar = mat2gray(scalar, options.Limits);
scalar = gray2ind(scalar, bit_depth);
scalar = ind2rgb(scalar, scalar_map);

% get alpha mask
mask = data_dic_save.strains(tt).(['roi_' options.View '_formatted']).mask;
alpha_mask = mask * options.Alpha;

% set data outside mask to NaN
scalar(repmat(~mask, [1 1 3])) = NaN;
raw_scalar(~mask) = NaN;



