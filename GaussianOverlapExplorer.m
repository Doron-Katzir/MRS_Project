function GaussianOverlapExplorer()
% GaussianOverlapExplorer
%
% Interactive visualization of two Gaussian basis functions.
%
% Shows how changing the Gaussian centers and width changes:
%   d      = |m2 - m1| / width
%   rhoB   = normalized basis overlap / cosine similarity
%   condG  = condition number of B'B
%
% This is meant as a visual helper for the two-Gaussian overlap parameter.

    %% Settings
    K = 160;
    x = linspace(0, 10, K).';

    % Default parameters
    m1_default = 4.5;
    m2_default = 5.5;
    w_default  = 0.55;

    %% Create UI
    fig = uifigure( ...
        'Name', 'Gaussian Basis Overlap Explorer', ...
        'Position', [200 200 900 520]);

    ax = uiaxes(fig, ...
        'Position', [300 120 560 360]);

    ax.XLabel.String = 'Spectral position x';
    ax.YLabel.String = 'Basis amplitude';
    ax.XLim = [0 10];
    ax.YLim = [0 1.1];
    grid(ax, 'on');
    hold(ax, 'on');

    % Sliders
    lblM1 = uilabel(fig, 'Position', [30 430 230 25]);
    sldM1 = uislider(fig, ...
        'Position', [30 410 220 3], ...
        'Limits', [0.5 9.5], ...
        'Value', m1_default);

    lblM2 = uilabel(fig, 'Position', [30 360 230 25]);
    sldM2 = uislider(fig, ...
        'Position', [30 340 220 3], ...
        'Limits', [0.5 9.5], ...
        'Value', m2_default);

    lblW = uilabel(fig, 'Position', [30 290 230 25]);
    sldW = uislider(fig, ...
        'Position', [30 270 220 3], ...
        'Limits', [0.1 2.5], ...
        'Value', w_default);

    % Output labels
    lblD = uilabel(fig, ...
        'Position', [30 200 240 30], ...
        'FontSize', 14, ...
        'FontWeight', 'bold');

    lblRho = uilabel(fig, ...
        'Position', [30 165 240 30], ...
        'FontSize', 14, ...
        'FontWeight', 'bold');

    lblCond = uilabel(fig, ...
        'Position', [30 130 240 30], ...
        'FontSize', 14, ...
        'FontWeight', 'bold');

    lblInterp = uilabel(fig, ...
        'Position', [30 40 830 60], ...
        'FontSize', 13, ...
        'WordWrap', 'on');

    btnReset = uibutton(fig, ...
        'push', ...
        'Text', 'Reset', ...
        'Position', [30 235 100 30], ...
        'ButtonPushedFcn', @(~,~) resetValues());

    %% Initial plot objects
    b1 = gaussianBasis(x, m1_default, w_default);
    b2 = gaussianBasis(x, m2_default, w_default);

    h1 = plot(ax, x, b1, 'LineWidth', 2.5);
    h2 = plot(ax, x, b2, 'LineWidth', 2.5);

    hM1 = xline(ax, m1_default, '--', 'm_1', 'LineWidth', 1.5);
    hM2 = xline(ax, m2_default, '--', 'm_2', 'LineWidth', 1.5);

    legend(ax, {'b_1', 'b_2'}, 'Location', 'northeast');

    %% Add slider callbacks
    sldM1.ValueChangingFcn = @(src,event) updatePlot(event.Value, sldM2.Value, sldW.Value);
    sldM2.ValueChangingFcn = @(src,event) updatePlot(sldM1.Value, event.Value, sldW.Value);
    sldW.ValueChangingFcn  = @(src,event) updatePlot(sldM1.Value, sldM2.Value, event.Value);

    sldM1.ValueChangedFcn = @(~,~) updatePlot(sldM1.Value, sldM2.Value, sldW.Value);
    sldM2.ValueChangedFcn = @(~,~) updatePlot(sldM1.Value, sldM2.Value, sldW.Value);
    sldW.ValueChangedFcn  = @(~,~) updatePlot(sldM1.Value, sldM2.Value, sldW.Value);

    updatePlot(m1_default, m2_default, w_default);

    %% Nested functions

    function updatePlot(m1, m2, w)
        b1 = gaussianBasis(x, m1, w);
        b2 = gaussianBasis(x, m2, w);

        B = [b1 b2];
        G = B.' * B;

        d = abs(m2 - m1) / w;

        rhoB = (b1.' * b2) / sqrt((b1.' * b1) * (b2.' * b2));
        condG = cond(G);

        h1.YData = b1;
        h2.YData = b2;

        hM1.Value = m1;
        hM2.Value = m2;

        lblM1.Text = sprintf('Center 1, m_1 = %.2f', m1);
        lblM2.Text = sprintf('Center 2, m_2 = %.2f', m2);
        lblW.Text  = sprintf('Width, w = %.2f', w);

        lblD.Text    = sprintf('d = |m_2 - m_1| / w = %.3f', d);
        lblRho.Text  = sprintf('\\rho_B = %.4f', rhoB);
        lblCond.Text = sprintf('cond(B''B) = %.3g', condG);

        if d < 0.75
            interp = 'Very small d: the two peaks are almost on top of each other. The basis functions are nearly dependent, so OLS separation becomes very unstable.';
        elseif d < 1.5
            interp = 'Small-to-moderate d: the peaks strongly overlap. They are still linearly independent, but noise can strongly contaminate the estimated amplitudes.';
        elseif d < 3
            interp = 'Moderate d: there is visible overlap, but the two basis directions are easier to distinguish.';
        else
            interp = 'Large d: the peaks are well separated. The basis overlap is low and OLS is numerically stable.';
        end

        lblInterp.Text = interp;

        title(ax, sprintf('Gaussian basis overlap: d = %.2f, \\rho_B = %.3f', d, rhoB));
        drawnow limitrate;
    end

    function resetValues()
        sldM1.Value = m1_default;
        sldM2.Value = m2_default;
        sldW.Value  = w_default;
        updatePlot(m1_default, m2_default, w_default);
    end
end

function b = gaussianBasis(x, center, width)
% gaussianBasis
% Gaussian basis normalized to max 1.

    b = exp(-(x - center).^2 ./ (2 * width^2));
    b = b ./ max(b);
end