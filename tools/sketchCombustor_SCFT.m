function [ myF ] = sketchCombustor_SCFT( p , varargin )
%SKETCHFLAME_L1Draws a sketch of the combustor in the physical domain (L1 coordinate system) and the image
%  domain
%
% Inputs:
%   - p      :  Struct with flame settings as returned from setUpPredefinedFlame(). Can allso be cell array of
%               structs!
%
% ////////////////////////////////////////////////////////
% // Thomas Steinbacher (steinbacher@tfd.mw.tum.de).    //
% // Created, 09.11.2017                                //
% // Last modified: 09.11.2017 by steinbacher           //
% ////////////////////////////////////////////////////////


%% Parse varargin
% Specify figure to plot in
ind = find(strcmpi(varargin,'fig'),1);
if ~isempty(ind)
  % User specified figure
  myF = varargin{ind+1};
else
  % new figure
  myF = [];
end

% Plot Stream line?
ind = find(strcmpi(varargin,'SL'),1);
if ~isempty(ind)
  % User specified stream line
  slDat = varargin{ind+1};
else
  % No stream line
  slDat = [];
end

% Plot setup in image domain?
ind = find(strcmpi(varargin,'image'),1);
if ~isempty(ind)
  % Yes
  doIm = 1;
else
  % Default: No
  doIm = 0;
end

% Zoom into image domain in order to visualize mapping?
ind = find(strcmpi(varargin,'zoomI'),1);
if ~isempty(ind)
  % Yes
  doZoom = 1;
else
  % Default: No
  doZoom = 0;
end

% Export as png to path
ind = find(strcmpi(varargin,'export2'),1);
if ~isempty(ind)
  % Yes
  expPath = varargin{ind+1};
  fileNameExt = varargin{ind+2};
else
  % Default: Current folder
  expPath = '.';
  fileNameExt = '';
end

% Include Stream lines potential flow?
ind = find(strcmpi(varargin,'doSLpot'),1);
if ~isempty(ind)
  % Yes
  doSLpot = 1;
else
  % Default: No
  doSLpot = 0;
end

% Move plotted origin for potential flow in image domain
ind = find(strcmpi(varargin,'moveOrigin4SLpot'),1);
if ~isempty(ind)
  % Yes
  move0SLpot = varargin{ind+1};
else
  % Default: No
  move0SLpot = 0;
end

% Mirror at real axis?
ind = find(strcmpi(varargin,'mirrorX1'),1);
if ~isempty(ind)
  % Yes
  mX1 = 1;
else
  % Default: No
  mX1 = 0;
end

% Should x2 size of figure be increased by a factor? Usefull when plotting mirrored and original domain
ind = find(strcmpi(varargin,'figFacX2'),1);
if ~isempty(ind)
  % Yes
  figFacX2 = varargin{ind+1};
else
  % Default: No
  figFacX2 = 1;
end

% Factor for line thicknesses
ind = find(strcmpi(varargin,'LineThick'),1);
if ~isempty(ind)
  % Specified by user
  LineThickFac = varargin{ind+1};
else
  % Default: 1
  LineThickFac = 1;
end

% Export to png?
ind = find(strcmpi(varargin,'savePNG'),1);
if ~isempty(ind)
  % Yes
  savePNG = 1;
else
  % Default: No
  savePNG = 0;
end

% Plot flame?
ind = find(strcmpi(varargin,'noFlame'),1);
if ~isempty(ind)
  % No
  plFlame = 0;
else
  % Default: Yes
  plFlame = 1;
end

% Define style for plotting flame
ind = find(strcmpi(varargin,'styleFlame'),1);
if ~isempty(ind)
  % Specified by user
  styFlm = varargin{ind+1};
else
  % Default
  styFlm = {':',2*LineThickFac,'r'};
end

%% Settings
% Check if input is cell array object,and if not convert it to one
if ~iscell(p)
  tmp = p; p = cell(1,1); p{1} = tmp;
end

% Styles {LineStyle,LineWidth,Color}
% Stlye casing
styCas = {'-',3*LineThickFac,'k'};
% Style flame  -> See above parsing varargin!
% Style symmetry line
stySym = {'-.',1*LineThickFac,'k'};
% Style stream line
stySL = {'--',2*LineThickFac,'b'};
% Potential Style stream line
stySLpot = {'-',0.7*LineThickFac,'g'};

% factor for physical domain (in order to improve plots)
fPhy = 1.2;

% compute Stream lines of potential flow if desired
if doSLpot
  % Image domain
  myBeta = [ 2 10 40 90 140 170 178 ] /180*pi;
%   xi1_SLpot = linspace(0,400,8000);
  xi1_SLpot = logspace(-2,2,1000)*4;
  xi2_SLpot = zeros(length(myBeta)+1,length(xi1_SLpot));
  x_Slpot = zeros(length(myBeta)+1,length(xi1_SLpot));
  for ii=1:length(myBeta)
    if myBeta(ii)~=pi/2
      xi2_SLpot(ii,:) = -sign(myBeta(ii)-pi/2)*tan(myBeta(ii)) * xi1_SLpot;
    else
      xi2_SLpot(ii,:) = xi1_SLpot;
    end
    % Map to physical L2 domain
    [ s ] = return_SCmap_SCFT( p{1} );
    x_Slpot(ii,:) = s.x_xi( -sign(myBeta(ii)-pi/2)*xi1_SLpot + 1i*xi2_SLpot(ii,:) );
  end
  
end

% Simplify access to R_a
if mX1
  % Mirror domain
  R_a = p{1}.R_a;
  mFac = -1;
else
  % L1
  R_a = 0;
  mFac = 1;
end

%% Plot
if ~doIm
  % Plot physical domain
  if isempty(myF)
    myF = figure('Position',[0 50 550/p{1}.Cr*0.7 550*figFacX2],'Color','w');hold on;
  else
    figure(myF)
  end
  % How much to plot left of area jump in terms of H_f?
  fracIn = 0.2;
  % How much to plot right of flame tip in terms of H_f?
  fracOut = 0.2;
  
  % Plot casings inlet
  plot([-fracIn*p{1}.H_flame 0],2*R_a+mFac*[1 1]*(p{1}.R_a-p{1}.R_flame),'LineStyle',styCas{1},'LineWidth'...
    ,styCas{2}*fPhy,'Color',styCas{3})
  
  % Plot casings combustion chamber
  plot([0 p{1}.H_flame+fracOut*p{1}.H_flame],2*R_a+mFac*[0 0],'LineStyle',styCas{1},'LineWidth'...
    ,styCas{2}*fPhy,'Color',styCas{3})
  plot([0 0],2*R_a+mFac*[(p{1}.R_a-p{1}.R_flame) 0],'LineStyle',styCas{1},'LineWidth'...
    ,styCas{2}*fPhy,'Color',styCas{3})
  
  % Plot symmetry line
  plot([-fracIn*p{1}.H_flame (1+fracOut)*p{1}.H_flame],2*R_a+mFac*[1 1]*p{1}.R_a,'LineStyle',stySym{1},'LineWidth'...
    ,stySym{2}*fPhy,'Color',stySym{3})
  
  % Plot flames
  if plFlame
    for ii=1:length(p)
      plot([0 p{ii}.H_flame],2*R_a+mFac*[(p{ii}.R_a-p{ii}.R_flame) p{ii}.R_a],'LineStyle',styFlm{1},'LineWidth'...
        ,styFlm{2}*fPhy,'Color',styFlm{3})
    end
  end
  
  % Plot Shear layer
  if ~isempty(slDat)
    % only plot up to end combusiton chamber
    inCC_SL = real(slDat.xSL)<= p{1}.H_flame+fracOut*p{1}.H_flame;
    plot(real(slDat.xSL(inCC_SL)),2*R_a+mFac*(p{1}.R_a-imag(slDat.xSL(inCC_SL))),'LineStyle',stySL{1},...
      'LineWidth',stySL{2}*fPhy,'Color',stySL{3})
%     plot(real(slDat.xSL(inCC_SL)),2*R_a+mFac*(p{1}.R_a-imag(slDat.xSL(inCC_SL))),'o--',...
%       'LineWidth',stySL{2}*fPhy,'Color',stySL{3})
  end
  
  % Include stream lines potential flow if desired
  if doSLpot
    for ii=1:length(myBeta)
      inDom_SLpot = real(x_Slpot(ii,:))>=-fracIn*p{1}.H_flame & real(x_Slpot(ii,:))<=(1+fracOut)*p{1}.H_flame;
      x1_Slpot_tmp = real(x_Slpot(ii,:));
      x2_Slpot_tmp = 2*R_a+mFac*imag(x_Slpot(ii,:));
      plot(x1_Slpot_tmp(inDom_SLpot),x2_Slpot_tmp(inDom_SLpot),'LineStyle',stySLpot{1},'LineWidth'...
        ,stySLpot{2}*fPhy,'Color',stySLpot{3} );
    end
  end
  
  set(gca,'Position',[0 0 1 1])
  axis off
  axis equal
  
  if savePNG
    print([expPath,filesep,p{1}.FlameName,'_physicalDomain',fileNameExt],'-dpng','-r200','-opengl')
  end
  
else
  % Plot image domain
  if isempty(myF)
    myF = figure('Position',[0 50 900 550],'Color','w');hold on;
  else
    figure(myF)
  end
  % How much to plot left of flame tip in terms of this value?
  fracIn = 0.1;
  % How much to plot right of shear layer end in terms of that value?
  fracOut = -0.3;
  
  % Compute flame points in L2 system
  x1F_flame = linspace(0,p{1}.L_flame,150);
  x1L2_flame = x1F_flame*cos(p{1}.alpha) + 1i*(p{1}.R_a-p{1}.R_flame + x1F_flame*sin(p{1}.alpha));
  % Map flame to image domain
  xi_flame =  SCmapInv_SCFT( x1L2_flame , p{1} );
  maxVal = max(real(xi_flame));
  % Map shear layer
  if ~isempty(slDat)
    xi_SL =  SCmapInv_SCFT( slDat.xSL , p{1} , 'L1' );
    % only plot until complex value is the same value as flame
    inCC_SL = imag(xi_SL)<=max(imag(xi_flame));
    xi_SL(~inCC_SL) = [];
    maxVal = max(real(xi_SL));
  end
  
  % Plot walls (real axis)
  plot([(1+fracIn)*(min(real(xi_flame))) 0 1 (pi/p{1}.R_flame)^2 (1+fracOut)*maxVal],[0 0 0 0 0],...
    'LineStyle',styCas{1},'LineWidth',styCas{2},'Color',styCas{3})
  
  % Plot flames
  if plFlame
    for ii=1:length(p)
      plot(xi_flame,'LineStyle',styFlm{1},'LineWidth'...
        ,styFlm{2},'Color',styFlm{3})
    end
  end
  
  % Plot Stream line
  if ~isempty(slDat)
    % Stream line
    plot(xi_SL,'LineStyle',stySL{1},'LineWidth',stySL{2},'Color',stySL{3})
    % Mirror Stream line
    plot(conj(xi_SL),'LineStyle',stySL{1},'LineWidth',stySL{2},'Color',stySL{3})
  end
  
  % Include stream lines potential flow if desired
  if doSLpot
    for ii=1:length(myBeta)
      plot(-sign(myBeta(ii)-pi/2)*xi1_SLpot+move0SLpot , xi2_SLpot(ii,:),'LineStyle',stySLpot{1},'LineWidth'...
        ,stySLpot{2}*fPhy,'Color',stySLpot{3} );
    end
  end
  
  % Limts + settings
  minXi1 = (1+fracIn)*(min(real(xi_flame)));
  maxXi1 = maxVal*1.28;
  maxXi2 = max(imag(xi_flame));
  
  if doZoom
    xlim([-5 15]);ylim([-3 7])
  else
    xlim([minXi1 maxXi1]);ylim([-1.1*0.7 1.1]*maxXi2)
  end
  set(gca,'Position',[0 0 1 1])
  axis off 
  
  if savePNG
    print([expPath,filesep,p{1}.FlameName,'_imageDomain',fileNameExt],'-dpng','-r200','-opengl')
  end
  
end


end

