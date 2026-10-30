function budgets = print_snr_budget(sys, targets)
%PRINT_SNR_BUDGET  Print the expected SNR of every target after every processing block.
%   Each target is evaluated at the time the antenna boresight crosses it, with
%   the CPI centred on it (see snr_budget.m).
budgets = [];
if isempty(targets), return; end
fprintf('\nExpected SNR per processing block [dB] (when the beam crosses each target)\n');
fprintf('%-4s %6s %6s %7s %6s %5s %4s | %7s %7s %6s %6s %6s %8s | %6s %6s\n', 'Tgt', 't[s]', 'RCS', ...
    'R[km]', 'El', 'Pulse', 'Beam', 'ADC', 'DDC', 'DBF', 'MF', 'MTI', 'Doppler', 'Req.', 'Margin');
for i = 1:numel(targets)
    t = illumination_time(sys, targets(i));
    b = snr_budget(sys, targets(i), t);
    fprintf('%-4s %6.2f %6.1f %7.2f %6.2f %5s %4d | %7.1f %7.1f %6.1f %6.1f %6.1f %8.1f | %6.1f %6.1f\n', ...
        b.name, t, b.rcs_m2, b.range_m/1e3, b.el_deg, b.pulse, b.beam, b.adc, b.ddc, b.dbf, ...
        b.mf, b.mti, b.doppler, b.required, b.margin);
    budgets = [budgets, b]; %#ok<AGROW>
end
fprintf(['ADC: one element, one real sample | DDC: one element, one complex sample | ', ...
         'Doppler: per row, %d row(s) (%.1f independent looks) integrated non-coherently\n'], ...
         sys.cpi.nRows, sys.dop.nLooks);
if isfield(sys.cfg.sim, 'swerling') && sys.cfg.sim.swerling == 1
    model = 'Swerling 1 target';
else
    model = 'Albersheim, steady target';
end
fprintf('Req.: SNR per row for Pd = 0.9 at Pfa = %.0e (%s)\n\n', sys.cfar.pfa, model);
end
