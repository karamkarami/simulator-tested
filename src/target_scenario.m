function targets = target_scenario(cfg)
%TARGET_SCENARIO  The targets of the scenario selected by cfg.scenario.type.
%   'six'       six_targets_scenario.m: six targets around the radar
%   'challenge' challenge_targets_scenario.m: thirteen targets that are hard for the radar (C11 - C13 maneuver)
%   'random'    random_targets_scenario.m: cfg.scenario.random.n targets (default 10)
%               drawn in the limits of cfg.scenario.random with the seed
%               cfg.scenario.random.seed ([] = cfg.sim.seed); the ranges are kept
%               inside 95% of the instrumented range (cfg.rx.maxRange_m)
switch lower(cfg.scenario.type)
    case 'six'
        targets = six_targets_scenario();
    case 'challenge'
        targets = challenge_targets_scenario();
    case 'random'
        r = cfg.scenario.random;
        seed = r.seed;
        if isempty(seed), seed = cfg.sim.seed; end
        r.range_m = min(r.range_m, 0.95 * cfg.rx.maxRange_m);
        if ~(r.n >= 1 && r.n == round(r.n)) || r.range_m(1) > r.range_m(2) || r.rcs_m2(1) <= 0 || ...
                r.rcs_m2(1) > r.rcs_m2(2) || r.alt_m(1) > r.alt_m(2) || r.speed_mps(1) > r.speed_mps(2)
            error('cfg.scenario.random: n must be a positive integer and every [min max] ordered (RCS > 0).');
        end
        targets = random_targets_scenario(r, seed);
    otherwise
        error('cfg.scenario.type must be ''six'', ''challenge'' or ''random''.');
end
end
