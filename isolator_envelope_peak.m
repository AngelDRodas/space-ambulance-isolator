function [a_pk, T, DR_rigid] = isolator_envelope_peak(P, pulse)
%ISOLATOR_ENVELOPE_PEAK  Peak of the synthetic input pulse for the chosen basis.
% Revision 2026-11-20 (Angel): the design pulse no longer comes from Apollo.
%
% P.pulse_basis = 'impulse'  (floor mount, current baseline, Angel 2026-11-20 PM)
%   The isolators sit on the capsule floor with no seat attenuation, so the floor
%   pulse is unknown. For a given velocity change, the isolated response grows as
%   the pulse gets shorter and levels off at the impulsive limit. The bounding
%   pulse is therefore a pulse of fixed short duration P.T_impulse that removes Vi
%   (it shows a solution exists for any pulse; it is not a design pulse):
%   half-sine a_pk = 1 + pi*Vi/(2*g*T), rectangle a_pk = 1 + Vi/(g*T).
%   At Vi = 8.9 m/s and T = 5 ms: half-sine 286 g.
%
% P.pulse_basis = 'envelope'  (seat mount)
%
% A crewed capsule must keep its seated crew within the Brinkley limits for
% healthy crew: NASA-STD-3001 Vol 2 Rev D [V2 6070] limits the injury risk
% criterion to beta <= 1 with the low-risk Dynamic Response limits
% (+X 35, -X -28; OCHMO-TB-024 Rev D, 2025; Somers et al. 2020). The worst
% landing the vehicle may deliver at a seat is therefore the pulse that, on a
% RIGID mount, drives the dynamic response exactly to that limit. This
% function finds the peak of that pulse by bisection, for the pulse shape
% given and the touchdown velocity P.Vi. The pulse always removes Vi:
%   half-sine T = pi*Vi/(2*A),  rectangle T = Vi/A,  A = (a_pk - 1)*g.
% Under the isolator, the same pulse must give beta <= 1 against the
% ill/injured limits (about 75 % of healthy).
%
%   [a_pk, T, DR_rigid] = isolator_envelope_peak(P)          uses P.pulse
%   [a_pk, T, DR_rigid] = isolator_envelope_peak(P, 'rect')  another shape
%
% With P.pulse_basis = 'fixed' (or measured data) it returns P.a_fixed.
% DR_rigid is returned for the 'envelope' basis only.
% Envelope results at 8.9 m/s: half-sine 29.7 g (49.7 ms), rectangle 25.4 g (37.2 ms).
if nargin < 2, pulse = P.pulse; end
T = NaN;  DR_rigid = NaN;
if strcmp(P.pulse_basis, 'fixed') || strcmp(pulse, 'data')
    a_pk = P.a_fixed;
    return
end
if strcmp(P.pulse_basis, 'impulse')
    T = P.T_impulse;
    switch pulse
        case 'halfsine', a_pk = 1 + pi*P.Vi/(2*P.g*T);
        case 'rect',     a_pk = 1 + P.Vi/(P.g*T);
        otherwise, error('Unknown pulse shape %s', pulse);
    end
    return
end

g   = P.g;
bt  = P.BTAB.(P.brinkley_axis);
BR.wn = bt{1};  BR.zeta = bt{2};  BR.lim = bt{3};      % HEALTHY limits [+ -]
BR.sgn = P.brinkley_sign;  BR.inc1g = P.brinkley_1g;

Q = P;  Q.pulse = pulse;
lo = 1.5;  hi = 100;                                     % [g] bracket
for it = 1:45
    Q.a_in_pk = (lo + hi)/2;
    [~, a, a0, ~, dt] = isolator_make_input(Q);
    if rigid_beta(a, a0, dt, g, BR) > 1, hi = Q.a_in_pk; else, lo = Q.a_in_pk; end
end
a_pk = lo;                                               % largest peak with beta <= 1
Q.a_in_pk = a_pk;
[~, a, a0, ~, dt, T] = isolator_make_input(Q);
[~, DR] = rigid_beta(a, a0, dt, g, BR);
DR_rigid = max(DR);
end

function [beta, DR] = rigid_beta(a, a0, dt, g, BR)
% Brinkley response of a crew member on a rigid mount (input = capsule reading),
% same exact discrete filter as the vibrations model.
    if BR.inc1g, A_in = BR.sgn*a;  A0 = BR.sgn*a0;
    else,        A_in = BR.sgn*(a - a0);  A0 = 0;
    end
    w  = BR.wn;  zt = BR.zeta;
    Ac = [0 1; -w^2 -2*zt*w];
    B  = [0; 1];
    Ad = expm(Ac*dt);
    Bd = Ac \ ((Ad - eye(2))*B);
    tr = trace(Ad);
    M2 = Ad - tr*eye(2);
    Cd = [w^2 0];
    du = A_in - A0;
    du = [(du(1:end-1) + du(2:end))/2; du(end)];        % interval average
    DR = (A0 + filter([0, Cd*Bd, Cd*M2*Bd], [1, -tr, det(Ad)], du))/g;
    beta = max([max(DR)/BR.lim(1), min(DR)/BR.lim(2), 0]);
end
