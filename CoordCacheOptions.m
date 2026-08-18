function cacheOptions = CoordCacheOptions(cfg, varargin)
%CoordCacheOptions Resolve the shared per-patient .coord cache settings.

    cacheOptions = struct();
    cacheOptions.enabled = true;
    cacheOptions.forceRefresh = false;
    cacheOptions.directory = "";

    rootDir = "";
    coordDir = "";

    if HasFieldOrProperty(cfg, "paths")
        cfgPaths = GetFieldOrProperty(cfg, "paths");

        if HasFieldOrProperty(cfgPaths, "rootDir")
            rootDir = string(GetFieldOrProperty(cfgPaths, "rootDir"));
        end
        if HasFieldOrProperty(cfgPaths, "coordDir")
            coordDir = string(GetFieldOrProperty(cfgPaths, "coordDir"));
        end
    end

    if HasFieldOrProperty(cfg, "load")
        cfgLoad = GetFieldOrProperty(cfg, "load");
        if HasFieldOrProperty(cfgLoad, "coordDir")
            coordDir = string(GetFieldOrProperty(cfgLoad, "coordDir"));
        end
    end

    if HasFieldOrProperty(cfg, "covariance")
        cfgCovariance = GetFieldOrProperty(cfg, "covariance");
        if HasFieldOrProperty(cfgCovariance, "coordRoot")
            coordDir = string(GetFieldOrProperty(cfgCovariance, "coordRoot"));
        elseif HasFieldOrProperty(cfgCovariance, "coordDir")
            coordDir = string(GetFieldOrProperty(cfgCovariance, "coordDir"));
        end
    end

    if HasFieldOrProperty(cfg, "cache")
        cfgCache = GetFieldOrProperty(cfg, "cache");

        if HasFieldOrProperty(cfgCache, "coord")
            cfgCoordCache = GetFieldOrProperty(cfgCache, "coord");

            if HasFieldOrProperty(cfgCoordCache, "enabled")
                cacheOptions.enabled = logical(GetFieldOrProperty(cfgCoordCache, "enabled"));
            end
            if HasFieldOrProperty(cfgCoordCache, "forceRefresh")
                cacheOptions.forceRefresh = logical(GetFieldOrProperty(cfgCoordCache, "forceRefresh"));
            end
            if HasFieldOrProperty(cfgCoordCache, "directory")
                cacheOptions.directory = string(GetFieldOrProperty(cfgCoordCache, "directory"));
            end
        end
    end

    p = inputParser;
    p.KeepUnmatched = true;
    p.addParameter('coordCacheEnabled', cacheOptions.enabled, ...
        @(x) islogical(x) && isscalar(x));
    p.addParameter('coordCacheForceRefresh', cacheOptions.forceRefresh, ...
        @(x) islogical(x) && isscalar(x));
    p.addParameter('coordCacheDirectory', cacheOptions.directory, ...
        @(x) ischar(x) || isstring(x));
    p.addParameter('coordRoot', coordDir, @(x) ischar(x) || isstring(x));
    p.addParameter('coordDir', coordDir, @(x) ischar(x) || isstring(x));
    parse(p, varargin{:});

    cacheOptions.enabled = logical(p.Results.coordCacheEnabled);
    cacheOptions.forceRefresh = logical(p.Results.coordCacheForceRefresh);
    cacheOptions.directory = string(p.Results.coordCacheDirectory);

    if strlength(string(p.Results.coordRoot)) > 0
        coordDir = string(p.Results.coordRoot);
    elseif strlength(string(p.Results.coordDir)) > 0
        coordDir = string(p.Results.coordDir);
    end

    if strlength(cacheOptions.directory) == 0
        if strlength(rootDir) == 0 && strlength(coordDir) > 0
            rootDir = string(fileparts(coordDir));
        end
        if strlength(rootDir) > 0
            cacheOptions.directory = string(fullfile(rootDir, "GeneratedCache", "CoordParsing"));
        end
    end

    if cacheOptions.enabled && strlength(cacheOptions.directory) == 0
        error(['Coord caching is enabled, but no cache directory could be resolved. ', ...
            'Set cfg.cache.coord.directory or cfg.paths.rootDir.']);
    end
end


function tf = HasFieldOrProperty(value, name)
    tf = (isstruct(value) && isfield(value, name)) || ...
        (isobject(value) && isprop(value, name));
end


function valueOut = GetFieldOrProperty(valueIn, name)
    valueOut = valueIn.(name);
end
