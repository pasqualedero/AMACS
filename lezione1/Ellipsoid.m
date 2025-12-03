classdef Ellipsoid < handle
    %ELLIPSOID class

    properties
        Q (:,:) {mustBeNumeric};
        q (:,1) {mustBeNumeric};       
    end

    properties (Dependent)
        n (1,1) {mustBeInteger};
    end

    properties (SetAccess = private)
        counterID (1,1) {mustBeNonnegative, mustBeInteger} ;
    end

    methods
        function obj = Ellipsoid(Q, optional)
            % Construct an instance of this class
            arguments
                Q 
                optional.q = zeros(size(Q,1),1) 
            end
            obj.Q = Q;
            obj.q = optional.q;
            obj.counterID = Ellipsoid.getNextID();
        end
        
        function boolean = isInternal(obj, x)
            % Verifica appartenenza Ellipsoid
            temp = (x-obj.q)'*obj.Q*(x-obj.q);
            boolean = temp <= 1;
        end

        function p = randompoint(obj)
            % Generates a random point uniformly *inside* the ellipsoid
            % using Cholensky Decomposition

            z = randn(obj.n, 1);
            z_ball = z / norm(z) * rand()^(1/obj.n);

            % Find the transformation matrix L such that L*L' = inv(Q)
            % This matrix maps the unit ball to the ellipsoid
            try
                % 'lower' ensures L*L' = inv(Q)
                L = chol(inv(obj.Q), 'lower');
            catch
                % Fallback if inv(Q) is not positive definite
                [X, E] = eig(inv(obj.Q));
                L = X * sqrt(E);
            end

            % Transform the unit ball point and shift by the center q
            p = L * z_ball + obj.q;
        end

        function [phi ,b] = supp(obj, z)
            phi=z'* obj.q+ sqrt(z'* inv(obj.Q)*z);
            b=inv(obj.Q)*z/ sqrt(z'* inv(obj.Q)*z)+ obj.q;
        end

        function points = plotEll(obj, options)
            arguments
                obj 
                options.Color {mustBeText} = 'b'
                options.Label {mustBeText} = ''
            end
            sens = 1e-4;
            angles = 0:sens:2*pi;
            points = zeros(2,length(angles));
            z = [cos(angles); sin(angles)];
            for i=1:length(angles)
                [~,points(:,i)] = obj.supp(z(:,i));
            end
            plot(points(1,:),points(2,:),'Color',options.Color);
            axis equal;
            if ~isempty(options.Label)
                text(points(1,end),points(2,end), options.Label, 'Interpreter','latex');
            end
        end

        % GET and SET methods

        function set.q(obj, value)
            if size(value,1) == obj.n
                obj.q = value;
            else
                error('q must have same size as x');
            end
        end

        function value = get.n(obj)
            value = size(obj.Q, 1);
        end

        function set.Q(obj, value)
            try chol(value);
                obj.Q = value;
            catch ME
                error('Q must be symmetric positive definite');
            end 
        end

    end

    methods (Static, Access=private)
        
        function id = getNextID()
            persistent instanceCounter;
            if isempty(instanceCounter)
                instanceCounter = 0;
            end
            id = instanceCounter;
            instanceCounter = instanceCounter + 1;   
        end

    end

end