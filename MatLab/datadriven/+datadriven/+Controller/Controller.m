classdef (Abstract) Controller < handle
    %CONTROLLER Abstract Controller class
    % 
    %   This class sets up the interface for the Abstract Controller
    % 
    %   ------------------------------------------------------------------------
    %   Copyright 2025 Vaibhav Gupta, DDMAC, EPFL (MIT License)
    %
    
    %% Abstract Methods
    methods (Abstract) 
        %FREQRESP Evaluates frequency response over a grid of frequencies
        [X, Y] = freqresp(obj, omegas)
    
        %EVALFR Evaluates frequency response at a single frequency
        [X, Y] = evalfr(obj, omega)
    
        %SS Convert from Controller object to MATLAB SS object
        K = ss(obj)
    end    
    
    %% Properties
    properties (SetAccess = protected)
        % Sample time
        %   It is specified as,
        %   * 0 for continuous-time systems.
        %   * A positive scalar representing the sampling period of a discrete-time 
        %     system in the time unit specified by the TimeUnit property.
        %   * -1 for a discrete-time system with an unspecified sample time.
        Ts (1, 1) double
        
        % Number of control signals
        nu (1, 1) double
    
        % Number of inputs to the controller
        ny (1, 1) double
    
        % Factorisation of the controller
        %   For SISO case, 'Left' and 'Right' are equivalent.
        factorisation (1, 1) datadriven.Controller.Factorisation
    end

end