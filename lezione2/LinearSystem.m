classdef LinearSystem < matlab.System
    % LINEARSYSTEM State-space plant model
    
    % Public, non-tunable properties
    properties (Nontunable)
        A;
        B;
        C;
        D;
    end

    % Internal states (automatically handled by MATLAB)
    properties (DiscreteState)
        x
    end

    % Pre-computed constants 
    properties (SetAccess = immutable)
        n
        m
    end

    properties
        x0;
    end

    methods
        % Constructor
        function obj = LinearSystem(varargin)
            setProperties(obj,nargin,varargin{:})
            obj.n = size(obj.A, 1);
            obj.m = size(obj.B, 2);
            if isempty(obj.x0)
                obj.x0 = zeros(length(obj.A),1);
            end
        end
        
        % Getters/Setters
        function set.A(obj,value)
            nn = size(value);
            if nn(1) == nn(2)
                obj.A = value;
            else
                error('A is not squared');
            end
        end

        function set.x0(obj, value)
            if length(value) ~= length(obj.A)
                error('x0 dimension error');
            else
                obj.x0 = value;
            end
        end
    end

    methods (Access = protected)
        function resetImpl(obj)
            % Initialize / Reset discrete states         
            obj.x = obj.x0;
        end

        function y = stepImpl(obj, u)            
            y = obj.C * obj.x + obj.D * u;
            obj.x = obj.A * obj.x + obj.B * u;
        end
    end
end