classdef llePlotter

    properties
        fig
        p % panel object
        t_evol
    end

    properties (Access = protected)
        linesMeanLight
        linesGenLight
        linesWF
        linesSP
    end

    methods
        function obj = llePlotter(psi, zeta, F0, t_evol)
            obj.fig = figure();
            obj.fig.Position(4) = 1.5*obj.fig.Position(4);

            p = panel();
            p.pack('v', 3)
            p(3).pack('h', 2)
            p.fontsize = 12;

            % Pump mode power and detuning profile
            p(1).select()
            yyaxis right
            fplot(zeta, [0, t_evol], 'k')
            yyaxis left
            set(gca, 'ColorOrder', lines(2))
            obj.linesMeanLight = plot(0, abs(mean(psi).^2), '.-');
            
            % Generated light panel
            p(2).select()
            set(gca, 'ColorOrder', lines(2))
            obj.linesGenLight = plot(0, sum(abs(psi - mean(psi)).^2), '.-');
            xlim([0,t_evol])

            
            % Time domain panel
            p(3, 1).select()
            obj.linesWF = plot(abs(psi));
            axis tight
            xlim([1, length(psi)])
            
            yyaxis right % plot the phase as well
            set(gca, 'ColorOrder', lines(2))
            phaseLines = plot(-tan(angle(psi)), '--');
            axis tight
            xlim([1, length(psi)])
            if ~isnumeric(zeta)
                FpumpMax = max(F0( -fminbnd(@(t) -max(F0(t)), 0, t_evol) ));
                zetaMax = zeta(-fminbnd(@(t) -abs(zeta(t)), 0, t_evol));
                ylim([-1,1] * max([1, zetaMax, FpumpMax]))
            end
            obj.linesWF = [obj.linesWF ; phaseLines];

            % Spectrum panel
            p(3, 2).select()
            obj.linesSP = plot(mag2db(abs(spectrumF(psi))));
            axis tight
            xlim([1, length(psi)])
            ylim([-200, 1])

            obj.t_evol = t_evol;
            obj.p = p;
            
            movegui(obj.fig, 'center');
        end

        function report(obj, t, psi)
            MF = mean(psi);
            
            newX = cellfun(@(X) [X, t], get(obj.linesMeanLight, 'XData'), 'UniformOutput', false);
            YVal = num2cell(abs(MF).^2, 1)';
            newY = cellfun(@(Y, y) [Y, y], get(obj.linesMeanLight, 'YData'), YVal, 'UniformOutput', false);
            set(obj.linesMeanLight, {'XData'}, newX, {'YData'}, newY);
            
            newX = cellfun(@(X) [X, t], get(obj.linesGenLight, 'XData'), 'UniformOutput', false);
            YVal = num2cell(sum(abs(psi - MF).^2), 1)';
            newY = cellfun(@(Y, y) [Y, y], get(obj.linesGenLight, 'YData'), YVal, 'UniformOutput', false);
            set(obj.linesGenLight, {'XData'}, newX, {'YData'}, newY);
            
            set(obj.linesWF, {'YData'}, num2cell([abs(psi), -tan(angle(psi))], 1)');
            set(obj.linesSP, {'YData'}, num2cell(mag2db(abs(spectrumF(psi))), 1)');
            drawnow()
        end

    end

end