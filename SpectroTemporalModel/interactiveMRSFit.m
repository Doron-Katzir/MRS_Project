function fig = interactiveMRSFit(basisFids, t, mask, ppmFit, sFit, ...
                                  metaboliteNames, thetaFit, cFit, betaFit)
% INTERACTIVEMRSFIT  Interactive complex-spectrum metabolite explorer.
%
%   fig = interactiveMRSFit(basisFids, t, mask, ppmFit, sFit, ...
%             metaboliteNames, thetaFit, cFit, betaFit)
%
% Display channels: Real / Imaginary / Both.  The *calculations* always use
% the full complex spectrum, irrespective of which channel is displayed.
%
% All metabolite amplitudes, linewidths and shifts start at the supplied
% original-fit values. Moving a slider changes only that metabolite; the
% other metabolites and the global Gaussian, phase and baseline are fixed.
% There is deliberately NO refitting during interaction.
%
% IMPORTANT: To start from an actual COMPLEX optimized fit, thetaFit,
% cFit and betaFit must come from fitting the complex data, not the old
% real-only objective. See the corresponding profiledResidual update.
%
% Inputs:
%   basisFids        N x M complex, unbroadened time-domain metabolite bases
%   t                N x 1 time coordinates (seconds)
%   mask             N x 1 logical fit-region mask
%   ppmFit           K x 1 ppm at mask
%   sFit             K x 1 COMPLEX measured spectrum at mask
%   metaboliteNames  M names
%   thetaFit         [gamma(1:M); delta(1:M); g; phi]
%                    gamma: Lorentzian FWHM (Hz); delta: shift (Hz)
%                    g: global Gaussian FWHM (Hz); phi: global phase (rad)
%   cFit             M x 1 real, nonnegative metabolite coefficients
%   betaFit          2 x 1 real or complex baseline coefficients
%                    constant + slope on the masked frequency grid
%
% Call after the final Strategy B fit for one patient:
%   interactiveMRSFit(basisFids, t, mask, ppmFit, sFit, ...
%       metaboliteNames, thetaInitial, cHatB, betaHatB);

    % ------------------------- Validate data -------------------------
    [N, M] = size(basisFids);
    t = t(:);
    mask = logical(mask(:));
    ppm = ppmFit(:);
    y = sFit(:);  % IMPORTANT: retain the complex measured signal
    names = string(metaboliteNames(:));
    thetaFit = thetaFit(:);
    cFit = cFit(:);
    betaFit = betaFit(:);

    % Convert four real baseline coefficients to two complex coefficients
    if numel(betaFit) == 4
        betaFit = betaFit(1:2) + 1i*betaFit(3:4);
    end

    assert(numel(t) == N && numel(mask) == N, ...
        't and mask must match the number of basis FID samples.');
    assert(nnz(mask) == numel(y) && numel(ppm) == numel(y), ...
        'ppmFit and sFit must match the true mask entries.');
    assert(numel(names) == M && numel(cFit) == M, ...
        'metaboliteNames and cFit must contain M entries.');
    assert(numel(thetaFit) == 2*M+2, ...
        'thetaFit must contain [gamma(1:M); delta(1:M); g; phi].');
    assert(numel(betaFit) == 2, ...
        'betaFit must have two complex or real coefficients.');
    assert(norm(y) > 0, 'The measured spectrum is zero.');
    assert(all(isfinite(cFit)) && all(cFit >= 0), ...
        'cFit must be real, finite and nonnegative.');

    gammaFit = thetaFit(1:M);
    deltaFit = thetaFit(M+1:2*M);
    g = thetaFit(2*M+1);
    phi = thetaFit(2*M+2);

    % Same baseline shape as profiledResidual: constant + linear slope.
    K = numel(y);
    baselineBasis = [ones(K,1), linspace(-1,1,K).'];
    baseline = baselineBasis * betaFit;  % potentially COMPLEX

    % Shared spectral transformations held fixed throughout exploration.
    gaussian = exp(-((pi*g*t).^2)/(4*log(2)));
    phase = exp(1i*phi);

    gammaNow = gammaFit;
    deltaNow = deltaFit;
    cNow = cFit;
    Bnow = complex(zeros(K,M));
    for k = 1:M
        Bnow(:,k) = makeColumn(k);
    end
    Boriginal = Bnow;
    sOriginal = Boriginal*cFit + baseline;
    originalError = 100 * norm(y-sOriginal) / norm(y);

    % If a fitted coefficient was zero, create a nonzero exploratory scale.
    % For those metabolites only, 0x = original fit and 1x is artificial.
    amplitudeRef = cFit;
    for k = 1:M
        if amplitudeRef(k) <= 0
            amplitudeRef(k) = 0.30*norm(y) / max(norm(Boriginal(:,k)), eps);
        end
    end

    selected = find(strcmpi(names,'Glu'),1);
    if isempty(selected), selected = 1; end

    % ------------------------- Colors / styles ----------------------
    % Saturated colors keep measured and original fit visible even where
    % they almost coincide. Original is DASHED and drawn over current.
    dataColor  = [0.00 0.40 0.95];  % blue: measured
    fitColor   = [1.00 0.34 0.04];  % orange: modified fit
    refColor   = [0.69 0.05 0.72];  % vivid purple: original fit
    residColor = [0.00 0.48 0.63];
    otherColor = [0.58 0.61 0.68];  % visible gray, not near-white
    baseColor  = [0.02 0.48 0.37];  % teal: constant+slope baseline
    selColor   = [0.86 0.12 0.25];  % red: selected metabolite

    % --------------------------- Figure -----------------------------
    fig = uifigure('Name','MRS Interactive Explorer - Full Complex Fit', ...
                   'Position',[70 65 1500 900], 'Color',[1 1 1]);
    mainGrid = uigridlayout(fig,[1 2]);
    mainGrid.ColumnWidth = {'1x',345};
    mainGrid.Padding = [12 12 12 12];
    mainGrid.ColumnSpacing = 12;

    plotGrid = uigridlayout(mainGrid,[3 1]);
    plotGrid.Layout.Row = 1;
    plotGrid.Layout.Column = 1;
    plotGrid.RowHeight = {'1x','2x','1.5x'};
    plotGrid.Padding = [0 0 0 0];
    plotGrid.RowSpacing = 8;

    axResidual = uiaxes(plotGrid); axResidual.Layout.Row = 1;
    axFit      = uiaxes(plotGrid); axFit.Layout.Row = 2;
    axContrib  = uiaxes(plotGrid); axContrib.Layout.Row = 3;

    residualOrig = y-sOriginal;

    hold(axResidual,'on');
    hResidualR = plot(axResidual,ppm,real(residualOrig),'-', ...
        'Color',residColor,'LineWidth',1.6);
    hResidualI = plot(axResidual,ppm,imag(residualOrig),':', ...
        'Color',residColor,'LineWidth',1.7);
    hOriginalResidualR = plot(axResidual,ppm,real(residualOrig),'--', ...
        'Color',refColor,'LineWidth',1.4);
    hOriginalResidualI = plot(axResidual,ppm,imag(residualOrig),'-.', ...
        'Color',refColor,'LineWidth',1.4);
    yline(axResidual,0,':','Color',[0.35 0.35 0.35], ...
        'HandleVisibility','off');
    ylabel(axResidual,'Residual');

    hold(axFit,'on');
    % Draw modified spectrum first, so the measured line and dashed
    % original line are not completely hidden by overlapping fit curves.
    hCurrentFitR = plot(axFit,ppm,real(sOriginal),'-', ...
        'Color',fitColor,'LineWidth',2.0);
    hCurrentFitI = plot(axFit,ppm,imag(sOriginal),':', ...
        'Color',fitColor,'LineWidth',2.0);
    hMeasuredR = plot(axFit,ppm,real(y),'-', ...
        'Color',dataColor,'LineWidth',1.45);
    hMeasuredI = plot(axFit,ppm,imag(y),':', ...
        'Color',dataColor,'LineWidth',1.45);
    % Keep the original fit dashed and high-contrast, drawn last.
    hOriginalFitR = plot(axFit,ppm,real(sOriginal),'--', ...
        'Color',refColor,'LineWidth',1.85);
    hOriginalFitI = plot(axFit,ppm,imag(sOriginal),'-.', ...
        'Color',refColor,'LineWidth',1.85);
    ylabel(axFit,'Spectral intensity');
    title(axFit,'Measured, original fit, and interactive fit');

    hold(axContrib,'on');
    hOtherR = gobjects(M,1);
    hOtherI = gobjects(M,1);
    for k = 1:M
        contribution = Bnow(:,k)*cNow(k);
        hOtherR(k) = plot(axContrib,ppm,real(contribution),'-', ...
            'Color',otherColor,'LineWidth',0.9,'HandleVisibility','off');
        hOtherI(k) = plot(axContrib,ppm,imag(contribution),':', ...
            'Color',otherColor,'LineWidth',0.9,'HandleVisibility','off');
    end
    % Legend proxies: one label for *all* gray curves, not 18 entries.
    hOtherLegendR = plot(axContrib,nan,nan,'-', ...
        'Color',otherColor,'LineWidth',1.2);
    hOtherLegendI = plot(axContrib,nan,nan,':', ...
        'Color',otherColor,'LineWidth',1.2);
    hBaselineR = plot(axContrib,ppm,real(baseline),'--', ...
        'Color',baseColor,'LineWidth',1.7);
    hBaselineI = plot(axContrib,ppm,imag(baseline),'-.', ...
        'Color',baseColor,'LineWidth',1.7);
    hSelectedR = plot(axContrib,ppm,real(Bnow(:,selected)*cNow(selected)),'-', ...
        'Color',selColor,'LineWidth',2.3);
    hSelectedI = plot(axContrib,ppm,imag(Bnow(:,selected)*cNow(selected)),':', ...
        'Color',selColor,'LineWidth',2.3);
    ylabel(axContrib,'Individual contribution');
    xlabel(axContrib,'Chemical shift (ppm)');

    axesArray = [axResidual,axFit,axContrib];
    for k = 1:numel(axesArray)
        set(axesArray(k),'XDir','reverse');
        grid(axesArray(k),'on');
        xlim(axesArray(k),[min(ppm),max(ppm)]);
    end
    axResidual.XTickLabel = {};
    axFit.XTickLabel = {};
    linkaxes(axesArray,'x');

    % ------------------------- Control panel ------------------------
    ctrl = uigridlayout(mainGrid,[16 1]);
    ctrl.Layout.Row = 1;
    ctrl.Layout.Column = 2;
    ctrl.RowHeight = {30,31,20,31,22,46,22,46,22,46,28,34,34,34,41,'1x'};
    ctrl.Padding = [8 8 8 8];
    ctrl.RowSpacing = 6;

    uilabel(ctrl,'Text','Explore fitted metabolites', ...
        'FontSize',17,'FontWeight','bold');

    dd = uidropdown(ctrl,'Items',cellstr(names),'ItemsData',1:M, ...
        'Value',selected,'ValueChangedFcn',@selectMetabolite);

    uilabel(ctrl,'Text','Display channel (fit is always complex)');
    viewDropdown = uidropdown(ctrl, ...
        'Items',{'Real','Imaginary','Both'},'Value','Both', ...
        'ValueChangedFcn',@viewChanged);

    uilabel(ctrl,'Text','Amplitude (multiplier of reference)');
    ampGroup = uigridlayout(ctrl,[1 2]);
    ampGroup.ColumnWidth = {'1x',83};
    ampGroup.Padding = [0 0 0 0];
    initialMultiplier = cNow(selected)/amplitudeRef(selected);
    ampSlider = uislider(ampGroup,'Limits',[0 4], ...
        'Value',min(4,initialMultiplier), ...
        'ValueChangingFcn',@amplitudeChanging, ...
        'ValueChangedFcn',@amplitudeChanged);
    ampSlider.Layout.Column = 1;
    ampNumber = uieditfield(ampGroup,'numeric','Limits',[0 4], ...
        'Value',min(4,initialMultiplier),'ValueDisplayFormat','%.2f', ...
        'ValueChangedFcn',@amplitudeEntered);
    ampNumber.Layout.Column = 2;

    uilabel(ctrl,'Text','Extra Lorentzian FWHM (Hz)');
    gammaGroup = uigridlayout(ctrl,[1 2]);
    gammaGroup.ColumnWidth = {'1x',83};
    gammaGroup.Padding = [0 0 0 0];
    gammaSlider = uislider(gammaGroup,'Limits',[0 120], ...
        'Value',gammaNow(selected), ...
        'ValueChangingFcn',@linewidthChanging, ...
        'ValueChangedFcn',@linewidthChanged);
    gammaSlider.Layout.Column = 1;
    gammaNumber = uieditfield(gammaGroup,'numeric','Limits',[0 120], ...
        'Value',gammaNow(selected),'ValueDisplayFormat','%.2f', ...
        'ValueChangedFcn',@linewidthEntered);
    gammaNumber.Layout.Column = 2;

    uilabel(ctrl,'Text','Frequency shift (Hz)');
    shiftGroup = uigridlayout(ctrl,[1 2]);
    shiftGroup.ColumnWidth = {'1x',83};
    shiftGroup.Padding = [0 0 0 0];
    shiftSlider = uislider(shiftGroup,'Limits',[-50 50], ...
        'Value',deltaNow(selected), ...
        'ValueChangingFcn',@shiftChanging, ...
        'ValueChangedFcn',@shiftChanged);
    shiftSlider.Layout.Column = 1;
    shiftNumber = uieditfield(shiftGroup,'numeric','Limits',[-50 50], ...
        'Value',deltaNow(selected),'ValueDisplayFormat','%.2f', ...
        'ValueChangedFcn',@shiftEntered);
    shiftNumber.Layout.Column = 2;

    showOthers = uicheckbox(ctrl,'Text','Show other metabolite contributions', ...
        'Value',true,'ValueChangedFcn',@othersChanged);

    uibutton(ctrl,'Text','Reset selected metabolite', ...
        'ButtonPushedFcn',@resetSelected);
    uibutton(ctrl,'Text','Reset all metabolites', ...
        'ButtonPushedFcn',@resetAll);
    uibutton(ctrl,'Text','Export current figure...', ...
        'ButtonPushedFcn',@exportSnapshot);

    coeffLabel = uilabel(ctrl,'Text','', ...
        'FontSize',12,'FontWeight','bold');
    note = uitextarea(ctrl,'Editable','off', ...
        'Value',{'Both real and imaginary residuals enter the error metric.', ...
                 'The green baseline is constant + slope, kept fixed.', ...
                 'Gray traces are individual OTHER metabolite contributions.', ...
                 'No refitting occurs when changing a metabolite.', ...
                 'Global Gaussian and global phase stay fixed.', ...
                 'For an originally zero coefficient, 1x is artificial.'});
    note.FontSize = 11;

    refreshControls();
    refreshPlots();
    refreshView();

    % --------------------- Spectral recomputation -------------------
    function b = makeColumn(j)
        transform = exp(-pi*t*gammaNow(j) + 1i*2*pi*t*deltaNow(j));
        fullSpec = phase * fftshift(fft(basisFids(:,j) .* transform .* gaussian));
        b = fullSpec(mask);  % IMPORTANT: retain complex spectral column
    end

    % --------------------------- Callbacks --------------------------
    function selectMetabolite(src,~)
        selected = src.Value;
        refreshControls();
        refreshPlots();
        refreshView();
    end

    function viewChanged(~,~)
        refreshView();
    end

    function refreshControls()
        factor = cNow(selected)/amplitudeRef(selected);
        ampSlider.Value = min(4,max(0,factor));
        ampNumber.Value = min(4,max(0,factor));
        gammaSlider.Value = gammaNow(selected);
        gammaNumber.Value = gammaNow(selected);
        shiftSlider.Value = deltaNow(selected);
        shiftNumber.Value = deltaNow(selected);
    end

    function amplitudeChanging(~,event), setAmplitude(event.Value); end
    function amplitudeChanged(src,~), setAmplitude(src.Value); end
    function amplitudeEntered(src,~), setAmplitude(src.Value); end
    function setAmplitude(multiplier)
        cNow(selected) = max(0,multiplier)*amplitudeRef(selected);
        ampSlider.Value = multiplier;
        ampNumber.Value = multiplier;
        refreshPlots();
    end

    function linewidthChanging(~,event), setLinewidth(event.Value); end
    function linewidthChanged(src,~), setLinewidth(src.Value); end
    function linewidthEntered(src,~), setLinewidth(src.Value); end
    function setLinewidth(value)
        gammaNow(selected) = value;
        gammaSlider.Value = value;
        gammaNumber.Value = value;
        Bnow(:,selected) = makeColumn(selected);
        refreshPlots();
    end

    function shiftChanging(~,event), setShift(event.Value); end
    function shiftChanged(src,~), setShift(src.Value); end
    function shiftEntered(src,~), setShift(src.Value); end
    function setShift(value)
        deltaNow(selected) = value;
        shiftSlider.Value = value;
        shiftNumber.Value = value;
        Bnow(:,selected) = makeColumn(selected);
        refreshPlots();
    end

    function othersChanged(~,~)
        refreshView();
    end

    % Update visibility and legends only on view/selection changes.
    function refreshView()
        rVisible = any(strcmp(viewDropdown.Value,{'Real','Both'}));
        iVisible = any(strcmp(viewDropdown.Value,{'Imaginary','Both'}));
        hReal = [hResidualR hOriginalResidualR hCurrentFitR hMeasuredR ...
                 hOriginalFitR hSelectedR hBaselineR];
        hImag = [hResidualI hOriginalResidualI hCurrentFitI hMeasuredI ...
                 hOriginalFitI hSelectedI hBaselineI];
        set(hReal,'Visible',onOff(rVisible));
        set(hImag,'Visible',onOff(iVisible));
        for j = 1:M
            visibleOther = showOthers.Value && (j ~= selected);
            hOtherR(j).Visible = onOff(visibleOther && rVisible);
            hOtherI(j).Visible = onOff(visibleOther && iVisible);
        end
        % NaN-data proxy lines are visible to the legend, not on the plot.
        hOtherLegendR.Visible = 'on';
        hOtherLegendI.Visible = 'on';

        resHandles = gobjects(0); resNames = {};
        fitHandles = gobjects(0); fitNames = {};
        contribHandles = gobjects(0); contribNames = {};
        if rVisible
            resHandles = [resHandles hResidualR hOriginalResidualR];
            resNames = [resNames {'Modified real','Original real'}];
            fitHandles = [fitHandles hMeasuredR hCurrentFitR hOriginalFitR];
            fitNames = [fitNames {'Measured real','Modified real','Original real'}];
            contribHandles = [contribHandles hSelectedR];
            contribNames = [contribNames {'Selected real'}];
            if showOthers.Value
                contribHandles = [contribHandles hOtherLegendR];
                contribNames = [contribNames {'Other metabolites real'}];
            end
            contribHandles = [contribHandles hBaselineR];
            contribNames = [contribNames {'Fixed baseline real'}];
        end
        if iVisible
            resHandles = [resHandles hResidualI hOriginalResidualI];
            resNames = [resNames {'Modified imag','Original imag'}];
            fitHandles = [fitHandles hMeasuredI hCurrentFitI hOriginalFitI];
            fitNames = [fitNames {'Measured imag','Modified imag','Original imag'}];
            contribHandles = [contribHandles hSelectedI];
            contribNames = [contribNames {'Selected imag'}];
            if showOthers.Value
                contribHandles = [contribHandles hOtherLegendI];
                contribNames = [contribNames {'Other metabolites imag'}];
            end
            contribHandles = [contribHandles hBaselineI];
            contribNames = [contribNames {'Fixed baseline imag'}];
        end
        legend(axResidual,resHandles,resNames,'Location','northeast');
        legend(axFit,fitHandles,fitNames,'Location','northeast');
        legend(axContrib,contribHandles,contribNames,'Location','northeast');
    end

    function str = onOff(value)
        if value, str = 'on'; else, str = 'off'; end
    end

    function refreshPlots()
        sCurrent = Bnow*cNow + baseline;
        residual = y-sCurrent;
        hCurrentFitR.YData = real(sCurrent);
        hCurrentFitI.YData = imag(sCurrent);
        hResidualR.YData = real(residual);
        hResidualI.YData = imag(residual);
        for j = 1:M
            contribution = Bnow(:,j)*cNow(j);
            hOtherR(j).YData = real(contribution);
            hOtherI(j).YData = imag(contribution);
        end
        selectedContribution = Bnow(:,selected)*cNow(selected);
        hSelectedR.YData = real(selectedContribution);
        hSelectedI.YData = imag(selectedContribution);
        % Do not rebuild all three legends on every slider tick.

        currentError = 100*norm(residual)/norm(y);
        title(axResidual,sprintf( ...
            'Complex residual | Current %.2f%% | Original %.2f%%', ...
            currentError,originalError));
        title(axContrib,sprintf( ...
            '%s | c=%.4g | linewidth=%.2f Hz | shift=%+.2f Hz', ...
            char(names(selected)),cNow(selected), ...
            gammaNow(selected),deltaNow(selected)));
        coeffLabel.Text = sprintf('Coefficient: %.5g   (fitted: %.5g)', ...
                                   cNow(selected),cFit(selected));
    end

    function resetSelected(~,~)
        cNow(selected) = cFit(selected);
        gammaNow(selected) = gammaFit(selected);
        deltaNow(selected) = deltaFit(selected);
        Bnow(:,selected) = Boriginal(:,selected);
        refreshControls();
        refreshPlots();
    end

    function resetAll(~,~)
        cNow = cFit;
        gammaNow = gammaFit;
        deltaNow = deltaFit;
        Bnow = Boriginal;
        refreshControls();
        refreshPlots();
    end

    function exportSnapshot(~,~)
        [filename,folder] = uiputfile( ...
            {'*.png','PNG image (*.png)';'*.pdf','PDF document (*.pdf)'}, ...
            'Export interactive MRS figure','MRS_Complex_Fit_Explorer.png');
        if isequal(filename,0), return; end
        try
            exportapp(fig,fullfile(folder,filename));
        catch err
            uialert(fig,err.message,'Could not export figure');
        end
    end
end
