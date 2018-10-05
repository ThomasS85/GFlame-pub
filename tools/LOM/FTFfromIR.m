function [ m ] = FTFfromIR( h , t_vec , varargin )
%FTFFROMIR Calculates the FTF from a given IR
%
%   Inputs:
%       - h       : IR as values over time
%       - t_vec   : Vector with time samples for h
%
%   Outputs:
%       - m       : Struct containing frequency response data
%             m.flameModel    : idpoly model of flame (FIR from estimated IR)
%             m.IR            : Estimated IR coefficients from flameModel
%             m.t_IR          : Time delay values for m.IR
%             m.mag           : Magnitude of Frequency response data (FR)
%             m.phase         : Phase of Frequency response data 
%             m.omega_vec     : Angular frequencies for FR data
%             m.comp          : FTF as complex number
%
% ////////////////////////////////////////////////////////
% // Thomas Steinbacher (steinbacher@tfd.mw.tum.de).    //
% // Created, 13.07.2015 as part of GFLAME 0.1 (ROM)    //
% // Last modified: 13.07.2015 by steinbacher           //
% ////////////////////////////////////////////////////////


%% Parse varargin
% Should phase start at 0 be forced?
ind = find(strcmpi(varargin,'noPhase0'),1);
if ~isempty(ind)
  % No
  forcePhase2zero = 0;
else
  % Yes
  forcePhase2zero = 1;
end


%% Check input
% h has to be row vector
[z,s] = size(h);
if z > s
    h = h';
end


%% Parameter
% Resolution frequency
nf = 1e3;


%% Create model and calculate IR

% Sampling time
Ts = t_vec(2) - t_vec(1);
% Maximum frequency to plot (Nyquist)
f_max = 0.5 * 1 / Ts;

% Model
m.flameModel = idpoly([],h*Ts,'Ts',Ts);

% Get IR (for comparison)
[m.IR,m.t_IR] = impulse(m.flameModel);

% Get Bode
% omega vector for plot
m.omega_vec = linspace(0,f_max,nf)*2*pi;
% Calculate mag and phase
[mag,phase] = bode(m.flameModel,m.omega_vec);
% Get vectors
m.mag = squeeze(mag);
m.phase = squeeze(phase);

% Convert phase to rad
m.phase = m.phase / 180 * pi;
% Substract first entry of phase from phase so that it start at 0
if forcePhase2zero
  m.phase = m.phase - m.phase(1);
end

% Convert to complex number
m.comp = ( cos(m.phase) + 1i*sin(m.phase) ) .* m.mag;

end

