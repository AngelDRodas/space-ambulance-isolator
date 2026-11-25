# Space Ambulance Stretcher Isolator Models

MATLAB models used to size a passive spring–damper isolator for a stretcher carrying an ill or injured patient in a crewed return capsule at splashdown. They support the thesis *[thesis title]* ([author], [university], [year]).

The patient lies supine with the back toward the heat shield. The landing load is chest-to-back (+X, eyeballs in), and injury risk is judged with the Brinkley dynamic response criterion against the ill or injured limits (β ≤ 1, DR_X ≤ +26.3 / ≥ −21.0), together with a 0.5 m stroke limit.

## Models

| Script | What it does |
|---|---|
| `isolator_inputs.m` | All inputs in one place: masses, landing pulse, criteria, sweep ranges, pitch–heave geometry and numerical settings. **Edit inputs here only.** |
| `isolator_make_input.m` | Builds the floor accelerometer record (synthetic pulse or measured data) and the time grid shared by every model. |
| `isolator_envelope_peak.m` | Sets the pulse peak for the chosen pulse basis (see below). |
| `ShockAttentuationSpaceAmbulanceBrinkley.m` | Single-mass model, mechanical-vibrations formulation, solved as an exact discrete filter. Sweeps spring rate and damping and selects a design. |
| `ShockAttentuationSpaceAmbulanceBondGraphs.m` | The same single-mass model derived as a bond graph and integrated with fixed-step RK4. Independent of the script above. |
| `compare_models.m` | Cross-validates the two single-mass formulations over the full sweep. |
| `PitchHeaveIsolator.m` | Pitch–heave model: rigid stretcher on head and foot isolators, with Brinkley evaluated at the patient's centre of mass. Searches equal and separately tuned isolators and runs its own verification checks. |

The two single-mass scripts are kept separate on purpose. Two independent derivations that agree are the verification. Merging them would let a shared mistake pass unnoticed.

### Equations

Single-mass model, with `q` the isolator compression and `a(t)` the floor accelerometer reading (1 g at rest, so gravity needs no separate term):

```
m q'' + c q' + k q = m a(t)          felt acceleration  a_p = (k q + c q')/m
```

The Brinkley X body model (ωn = 62.8 rad/s, ζ = 0.20) is driven by `a_p` and gives DR_X and β. The pitch–heave model adds the pitch of the rigid stretcher, with the felt acceleration resolved at the patient's centre of mass.

## How to run

Requires MATLAB R2018b or later (base functions only, no toolboxes). It also runs in GNU Octave.

1. Put all `.m` files in one folder.
2. Set the inputs in `isolator_inputs.m`.
3. Run `ShockAttentuationSpaceAmbulanceBrinkley.m` (writes `results_vibrations.mat`).
4. Run `ShockAttentuationSpaceAmbulanceBondGraphs.m` (writes `results_bondgraph.mat`).
5. Run `compare_models.m` (writes `model_comparison.csv` and `model_comparison_table.tex`).
6. Run `PitchHeaveIsolator.m` (writes `results_pitchheave.mat`).

## Pulse basis

`P.pulse_basis` in `isolator_inputs.m`:

- `'impulse'` (default): a 5 ms half-sine that removes the touchdown velocity (8.9 m/s, 286 g peak). For a linear isolator, any pulse that stays above 1 g and removes the same velocity gives a response no larger than this impulsive limit. Passing it therefore shows that a solution exists for any floor pulse. It is a bounding case, not a design pulse.
- `'envelope'`: the half-sine that drives a crew member on a rigid mount to the healthy Brinkley limit DR_X = +35 (29.7 g, 49.7 ms at 8.9 m/s). Use this for a stretcher mounted at a seat interface.
- `'fixed'`: a fixed peak `P.a_fixed`.

To use a measured floor accelerometer record, set `P.pulse = 'data'` and point `P.data_file` to a two-column CSV (time in s, acceleration).

## Results (impulse basis, GNU Octave 8.4)

| Design | Spring rate | Damping | Worst-case margin |
|---|---|---|---|
| Pitch–heave, equal isolators (selected) | 19,162 N/m each | 884 N·s/m each | 9.3 % |
| Single-mass optimum | 38,195 N/m | 1,550 N·s/m | 11.8 % |

For the selected design, the lightest patient (102.9 kg isolated mass) has β = 0.90 and the heaviest (170.5 kg) has a stroke of 0.45 m. Without an isolator, β = 1.67.

The two single-mass formulations agree to better than 4 × 10⁻⁵ % over 80,000 cases, and the difference falls by a factor of four when the time step is halved.

## Assumptions and limitations

Rigid capsule floor, single axis (chest-to-back), patient rigid on the stretcher (vacuum mattress not modelled), linear springs and dampers, and isolators that stay attached on rebound. The Brinkley limits were derived for seated, harnessed crew. See the thesis for the full discussion.

## Citation

If you use this code, please cite the archived version: [[Zenodo DOI](https://doi.org/10.5281/zenodo.24518527)].

## License

MIT
