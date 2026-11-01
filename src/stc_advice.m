function s = stc_advice(cfg)
%STC_ADVICE  What the STC can do against an ADC overload (text for the warnings).
%   s = STC_ADVICE(cfg) returns 'switch the STC on (...)' when cfg.rx.stc.enable
%   is false, else 'raise the STC attenuation (...)'.
if cfg.rx.stc.enable
    s = 'raise the STC attenuation (cfg.rx.stc.curve, cfg.rx.stc.curveLong)';
else
    s = 'switch the STC on (cfg.rx.stc.enable = true)';
end
end
