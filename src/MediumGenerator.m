classdef MediumGenerator < handle

    properties
        idx_frac1
        separation_mult
        width_scale
        C1_amp
        C2_amp
        cont
        autoPlot = true   
        
    end

    methods
        %% =========================
        % Constructor
        %% =========================
        function obj = MediumGenerator(idx_frac1, separation_mult, ...
                                       width_scale, C1_amp, C2_amp)

            if nargin < 1 || isempty(idx_frac1), idx_frac1 = 7/50; end
            if nargin < 2 || isempty(separation_mult), separation_mult = 1; end
            if nargin < 3 || isempty(width_scale), width_scale = 1; end
            if nargin < 4 || isempty(C1_amp), C1_amp = 0.6; end
            if nargin < 5 || isempty(C2_amp), C2_amp = 1.0; end

            obj.idx_frac1 = idx_frac1;
            obj.separation_mult = separation_mult;
            obj.width_scale = width_scale;
            obj.C1_amp = C1_amp;
            obj.C2_amp = C2_amp;
        end

        %% =========================
        % BASE MEDIUM (Case 1)
        %% =========================
        function [C0, C1, x, meta] = medium1(obj, n, h, lam, r, c_base)

            [C0, C1, x, meta] = obj.baseMedium(n, h, lam, r, c_base);
            obj.cont = (C1.^2 - C0.^2)./ (C1.*C0);
            obj.plotMedium(x, C0, C1, 'Medium a');

        end


        %% =========================
        % MEDIUM 2 (bias + slope + reflection)
        %% =========================
        function [C0, C1, x, meta] = medium2(obj, n, h, lam, r, c_base)

            [C0, C1, x, meta] = obj.baseMedium(n, h, lam, r, c_base);

            % Bias
            C0 = C0 + 0.2;
            C1 = C1 + 0.2;

            % Linear trend
            C0 = C0 + x/25;
            C1 = C1 + x/25;

            % Reflection transform
            C1 = 2*C0 - C1;
            
            obj.plotMedium(x, C0, C1, 'Medium b');

        end


        %% =========================
        % MEDIUM 3 (smooth rho-based nonlinear transform)
        %% =========================
        function [C0, C1, x, meta] = medium3(obj, n, h, lam, r, c_base, amp, sigma)

            [C0, C1, x, meta] = obj.baseMedium(n, h, lam, r, c_base);
            % Bias
            %C0 = C0 + 0.2;

            % Linear trend
            C0 = C0 + x/25;

            %sigma = 0.5e-2;

            rho = @(x) amp * exp(-(x.^4) ./ (2*sigma^2));

            idx1 = meta.idx1;
            idx2 = meta.idx2;
            grid = x;

            obj.cont = rho(grid - 1 * idx1 * h) + 0* rho(grid - 1 * idx2 * h);

            C1 = 0.5 * C0 .* (obj.cont + sqrt(obj.cont.^2 + 4));
            
            obj.plotMedium(x, C0, C1, 'Medium c');

        end

        %% =========================
        % MEDIUM 4 (nonlinear increasing background)
        %% =========================
        function [C0, C1, x, meta] = medium4(obj, n, h, lam, r, c_base)
        
            if nargin < 6 || isempty(c_base)
                c_base = 1.0;
            end
        
            % Build base medium first (for anomaly placement)
            [C0, C1, x, meta] = obj.baseMedium(n, h, lam, r, c_base);
        
            % Nonlinear slow background growth
            alpha = 0.5;   % controls strength of increase
            C0 = c_base * (1 + alpha * x.^2);
        
            % Keep anomalies but scaled relative to new background
            C1 = C0;
        
            idx1 = meta.idx1;
            idx2 = meta.idx2;
        
            end1 = min(n, idx1 + meta.width1 - 1);
            end2 = min(n, idx2 + meta.width2 - 1);
        
            C1(idx1:end1) = C0(idx1:end1) * obj.C1_amp;
            C1(idx2:end2) = C0(idx2:end2) * obj.C2_amp;
            % % Update the medium properties based on the new parameters
            % C0 = C0 .* (1 + 0.1 * sin(2 * pi * x));  % Example of a sinusoidal modulation
            % C1 = C1 .* (1 + 0.1 * cos(2 * pi * x));  % Example of a cosine modulation
            obj.plotMedium(x, C0, C1, 'Medium d');

        end



        %% =========================
        % PRIVATE BASE BUILDER
        %% =========================
        function [C0, C1, x, meta] = baseMedium(obj, n, h, lam, r, c_base)

            if nargin < 6 || isempty(c_base), c_base = 1.0; end

            N = n;
            Le = (N-1) * h;
            x = linspace(0.0, Le, N).';

            C0 = c_base * ones(N,1);
            C1 = c_base * ones(N,1);

            idx1 = floor(N * obj.idx_frac1);
            lam = c_base * lam;
            sep = round(obj.separation_mult * r * lam / max(h, 1e-15));
            idx2 = idx1 + sep;

            width1 = max(1, round(obj.width_scale * r * lam / max(h, 1e-15)));
            width2 = 1;

            idx1 = max(1, min(N, idx1 + 1));
            idx2 = max(1, min(N, idx2 + 1));

            end1 = min(N, idx1 + width1 - 1);
            end2 = min(N, idx2 + width2 - 1);

            C1(idx1:end1) = obj.C1_amp;
            C1(idx2:end2) = obj.C2_amp;

            meta = struct( ...
                'idx1', idx1, ...
                'idx2', idx2, ...
                'width1', width1, ...
                'width2', width2, ...
                'sep', sep ...
            );

        end

                %% =========================
        % Internal plotting utility
        %% =========================
        function plotMedium(obj, x, C0, C1, nameStr)
        
            if ~obj.autoPlot
                return
            end
        
            figure('Color','white','Position',[200 200 900 350]);
            hold on;
        
            plot(x, C0, 'k-', 'LineWidth', 2);
            plot(x, C1, 'r-', 'LineWidth', 1.5);
        
            xlabel('x');
            ylabel('Wave speed');
            title(nameStr);
        
            legend({'C0 (background)','C1 (perturbed)'}, ...
                    'Location','best');
        
            grid on;
            box on;
        
        end


    end
end
