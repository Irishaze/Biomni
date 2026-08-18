def analyze_scaffold_degradation_kinetics(
    time_points,
    mass_retention_data,
    degradation_model="first_order",
    output_dir="./",
):
    """Fits scaffold/coating mass-retention decay over time, for one or more layers.

    Fits a simplified degradation model to mass-retention-vs-time data for a
    biodegradable vascular scaffold (e.g. PCL/PLA/PLGA), deriving a decay
    constant, half-life, and R^2 per layer. These are simplified engineering
    approximations intended for exploratory design comparisons, not validated
    predictors of in vivo degradation.

    Parameters
    ----------
    time_points : list or numpy.ndarray
        Time points at which mass retention was measured, in weeks
    mass_retention_data : list, numpy.ndarray, or dict
        Measured mass retention (%) at each time point. Pass a 1D list/array
        for a single layer, or a dict mapping layer name to its own
        list/array of mass retention values (measured at the same
        time_points) for multiple layers
    degradation_model : str, optional
        Degradation model to fit: "first_order" (exponential decay,
        mass(t) = 100 * exp(-k*t)) or "autocatalytic" (logistic decay used as
        a simplified approximation of autocatalytic hydrolysis,
        mass(t) = 100 / (1 + exp(k*(t - t_mid)))) (default: "first_order")
    output_dir : str, optional
        Directory to save output files (default: "./")

    Returns
    -------
    str
        Research log summarizing each layer's fitted decay parameters,
        half-life, R^2, and the saved output file paths

    """
    import os
    from datetime import datetime

    import matplotlib

    matplotlib.use("Agg")  # Use non-interactive backend
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    from scipy.optimize import curve_fit

    os.makedirs(output_dir, exist_ok=True)

    if degradation_model not in ("first_order", "autocatalytic"):
        return f"Error: Unknown degradation_model '{degradation_model}'. Use 'first_order' or 'autocatalytic'."

    time_points = np.array(time_points, dtype=float)

    if isinstance(mass_retention_data, dict):
        layers = {name: np.array(vals, dtype=float) for name, vals in mass_retention_data.items()}
    else:
        layers = {"Layer 1": np.array(mass_retention_data, dtype=float)}

    def first_order(t, k):
        return 100 * np.exp(-k * t)

    def autocatalytic(t, k, t_mid):
        return 100 / (1 + np.exp(k * (t - t_mid)))

    log = []
    log.append(f"Scaffold Degradation Kinetics Analysis - {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    log.append(f"Model: {degradation_model}")
    log.append(
        "NOTE: This is a simplified engineering approximation for exploratory design "
        "comparisons, not a validated predictor of in vivo degradation.\n"
    )

    results = []
    fitted_curves = {}

    for layer_name, mass_data in layers.items():
        log.append(f"Layer: {layer_name}")
        try:
            if degradation_model == "first_order":
                params, _ = curve_fit(first_order, time_points, mass_data, p0=[0.05], bounds=(0, np.inf))
                k = params[0]
                y_pred = first_order(time_points, k)
                half_life = np.log(2) / k if k > 0 else float("inf")
                t_mid = None
            else:
                t_mid_guess = time_points[np.argmin(np.abs(mass_data - 50))]
                params, _ = curve_fit(
                    autocatalytic,
                    time_points,
                    mass_data,
                    p0=[0.5, t_mid_guess],
                    bounds=([0, 0], [np.inf, np.inf]),
                )
                k, t_mid = params
                y_pred = autocatalytic(time_points, k, t_mid)
                half_life = t_mid

            ss_total = np.sum((mass_data - np.mean(mass_data)) ** 2)
            ss_residual = np.sum((mass_data - y_pred) ** 2)
            r2 = 1 - (ss_residual / ss_total) if ss_total > 0 else float("nan")

            fitted_curves[layer_name] = (mass_data, y_pred)
            results.append(
                {
                    "Layer": layer_name,
                    "Model": degradation_model,
                    "k": k,
                    "t_mid_weeks": t_mid,
                    "Half_life_weeks": half_life,
                    "R2": r2,
                }
            )
            log.append(f"  - Decay constant k = {k:.4f}")
            if t_mid is not None:
                log.append(f"  - Fitted midpoint t_mid = {t_mid:.2f} weeks")
            log.append(f"  - Half-life = {half_life:.2f} weeks")
            log.append(f"  - R^2 = {r2:.4f}")
        except Exception as e:
            log.append(f"  - Fitting failed: {str(e)}")
            results.append(
                {
                    "Layer": layer_name,
                    "Model": degradation_model,
                    "k": None,
                    "t_mid_weeks": None,
                    "Half_life_weeks": None,
                    "R2": None,
                }
            )
        log.append("")

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    plt.figure(figsize=(10, 6))
    for layer_name, (mass_data, y_pred) in fitted_curves.items():
        (line,) = plt.plot(time_points, mass_data, "o", label=f"{layer_name} (data)")
        order = np.argsort(time_points)
        plt.plot(
            time_points[order],
            y_pred[order],
            "--",
            color=line.get_color(),
            label=f"{layer_name} (fit)",
        )
    plt.xlabel("Time (weeks)")
    plt.ylabel("Mass Retention (%)")
    plt.title("Scaffold Mass Retention vs. Time")
    plt.legend()
    plt.grid(True, linestyle="--", alpha=0.7)
    plot_path = os.path.join(output_dir, f"degradation_kinetics_{timestamp}.png")
    plt.savefig(plot_path, dpi=300, bbox_inches="tight")
    plt.close()

    results_df = pd.DataFrame(results)
    csv_path = os.path.join(output_dir, f"degradation_kinetics_results_{timestamp}.csv")
    results_df.to_csv(csv_path, index=False)

    log.append("Files Generated:")
    log.append(f"  - Plot: {plot_path}")
    log.append(f"  - Results table: {csv_path}")

    return "\n".join(log)


def analyze_vascular_graft_mechanical_compliance(
    pressure_data,
    diameter_data,
    systolic_pressure=120,
    diastolic_pressure=80,
    burst_pressure=None,
    target_compliance_range=(4.4, 5.9),
    target_safety_factor=3.0,
    output_dir="./",
):
    """Computes vascular graft compliance, distensibility, and safety factor.

    Computes compliance (%/100mmHg, the conventional unit in the vascular
    graft literature since typical pulse pressure is ~40 mmHg) and a linear
    pressure-diameter distensibility slope from pressure-diameter test data,
    and (if burst pressure is supplied) a mechanical safety factor, comparing
    both against target values for a compliance-matched, mechanically safe
    graft.

    Parameters
    ----------
    pressure_data : list or numpy.ndarray
        Applied intraluminal pressure at each measurement, in mmHg
    diameter_data : list or numpy.ndarray
        Measured graft diameter at each corresponding pressure, in mm
    systolic_pressure : float, optional
        Physiological systolic pressure used for the compliance calculation,
        in mmHg (default: 120)
    diastolic_pressure : float, optional
        Physiological diastolic pressure used for the compliance calculation,
        in mmHg (default: 80)
    burst_pressure : float, optional
        Measured graft burst pressure, in mmHg. If provided, used to compute
        a mechanical safety factor (burst_pressure / systolic_pressure)
        (default: None)
    target_compliance_range : tuple of float, optional
        Target compliance range for a compliance-matched graft, in
        %/100mmHg (default: (4.4, 5.9), spanning reported values for
        saphenous vein (~4.4, the clinical gold-standard graft) to native
        femoral artery (~5.9))
    target_safety_factor : float, optional
        Minimum acceptable ratio of burst pressure to systolic pressure
        (default: 3.0)
    output_dir : str, optional
        Directory to save output files (default: "./")

    Returns
    -------
    str
        Research log with computed compliance, distensibility, safety
        factor, pass/fail flags against the target values, and saved output
        file paths

    """
    import os
    from datetime import datetime

    import matplotlib

    matplotlib.use("Agg")  # Use non-interactive backend
    import matplotlib.pyplot as plt
    import numpy as np

    os.makedirs(output_dir, exist_ok=True)

    pressure_data = np.array(pressure_data, dtype=float)
    diameter_data = np.array(diameter_data, dtype=float)

    if len(pressure_data) < 2:
        return "Error: At least 2 pressure-diameter data points are required."

    order = np.argsort(pressure_data)
    pressure_sorted = pressure_data[order]
    diameter_sorted = diameter_data[order]

    log = []
    log.append(f"Vascular Graft Mechanical Compliance Analysis - {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    log.append(
        "NOTE: This is a simplified engineering approximation for exploratory design "
        "comparisons, not a validated predictor of in vivo mechanical behavior.\n"
    )

    try:
        slope, intercept = np.polyfit(pressure_sorted, diameter_sorted, 1)
    except Exception as e:
        return f"Error fitting pressure-diameter relationship: {str(e)}"

    diameter_at_systolic = slope * systolic_pressure + intercept
    diameter_at_diastolic = slope * diastolic_pressure + intercept

    if diameter_at_diastolic <= 0 or systolic_pressure <= diastolic_pressure:
        return "Error: Invalid pressure range or non-physical fitted diameter at diastolic pressure."

    # Fractional diameter change per mmHg, expressed as %/100mmHg (the
    # conventional unit in the vascular graft literature): one factor of 100
    # converts the fractional strain to a percentage, the second rescales
    # from "per 1 mmHg" to "per 100 mmHg" (typical pulse pressure ~40 mmHg,
    # so "per mmHg" values are inconveniently small and not how compliance
    # is reported in practice).
    compliance = (
        (diameter_at_systolic - diameter_at_diastolic)
        / diameter_at_diastolic
        / (systolic_pressure - diastolic_pressure)
        * 1e4
    )

    compliance_pass = target_compliance_range[0] <= compliance <= target_compliance_range[1]

    log.append(f"Linear fit: diameter (mm) = {slope:.6f} * pressure (mmHg) + {intercept:.4f}")
    log.append(f"Distensibility slope (dD/dP): {slope:.6f} mm/mmHg")
    log.append(f"Fitted diameter at diastolic pressure ({diastolic_pressure} mmHg): {diameter_at_diastolic:.3f} mm")
    log.append(f"Fitted diameter at systolic pressure ({systolic_pressure} mmHg): {diameter_at_systolic:.3f} mm")
    log.append(f"Compliance: {compliance:.3f} %/100mmHg")
    log.append(
        f"Target compliance range: {target_compliance_range[0]}-{target_compliance_range[1]} %/100mmHg -> "
        f"{'PASS' if compliance_pass else 'FAIL'}"
    )

    safety_factor = None
    if burst_pressure is not None:
        safety_factor = burst_pressure / systolic_pressure
        safety_pass = safety_factor >= target_safety_factor
        log.append(f"Burst pressure: {burst_pressure} mmHg")
        log.append(f"Safety factor (burst pressure / systolic pressure): {safety_factor:.2f}")
        log.append(
            f"Target safety factor: >= {target_safety_factor} -> {'PASS' if safety_pass else 'FAIL'}"
        )
    else:
        log.append("Burst pressure not provided; safety factor not computed.")

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    plt.figure(figsize=(8, 6))
    plt.plot(pressure_sorted, diameter_sorted, "o", label="Measured data")
    fit_line_p = np.linspace(min(pressure_sorted.min(), diastolic_pressure), max(pressure_sorted.max(), systolic_pressure), 100)
    plt.plot(fit_line_p, slope * fit_line_p + intercept, "--", label="Linear fit")
    plt.axvline(diastolic_pressure, color="gray", linestyle=":", alpha=0.7, label="Diastolic")
    plt.axvline(systolic_pressure, color="black", linestyle=":", alpha=0.7, label="Systolic")
    plt.xlabel("Pressure (mmHg)")
    plt.ylabel("Diameter (mm)")
    plt.title("Pressure-Diameter Relationship")
    plt.legend()
    plt.grid(True, linestyle="--", alpha=0.7)
    plot_path = os.path.join(output_dir, f"pressure_diameter_{timestamp}.png")
    plt.savefig(plot_path, dpi=300, bbox_inches="tight")
    plt.close()

    log.append(f"\nFiles Generated:\n  - Plot: {plot_path}")

    return "\n".join(log)


def assess_graft_diameter_thrombosis_risk(
    graft_diameter_mm,
    blood_flow_rate_ml_s,
    blood_viscosity_cp=3.5,
    surface_hydrophobicity_score=None,
    diameter_sweep_range=(2, 10),
    output_dir="./",
):
    """Estimates diameter-dependent wall shear stress and a heuristic thrombosis-risk score.

    Computes wall shear stress (Poiseuille flow) and surface-area-to-volume
    ratio for a cylindrical vascular graft at the given diameter, and
    combines them into a heuristic patency-risk score reflecting the known
    scaling effect where smaller-caliber grafts experience both higher wall
    shear stress and higher surface contact per unit blood volume. This is a
    simplified heuristic for exploratory design comparisons, not a validated
    clinical predictor of thrombosis.

    Parameters
    ----------
    graft_diameter_mm : float
        Inner diameter of the graft, in mm
    blood_flow_rate_ml_s : float
        Volumetric blood flow rate through the graft, in mL/s
    blood_viscosity_cp : float, optional
        Blood dynamic viscosity, in centipoise (default: 3.5, a typical
        whole-blood value)
    surface_hydrophobicity_score : float, optional
        Optional surface hydrophobicity score from 0 (hydrophilic) to 1
        (hydrophobic), applied as a linear multiplier on the risk score
        (risk *= 1 + score). If None, this modifier is not applied
        (default: None)
    diameter_sweep_range : tuple of float, optional
        Diameter range (mm) to sweep for the comparison plot (default: (2, 10))
    output_dir : str, optional
        Directory to save output files (default: "./")

    Returns
    -------
    str
        Research log with computed wall shear stress, surface-area-to-volume
        ratio, heuristic risk score and category, and saved output file
        paths

    """
    import os
    from datetime import datetime

    import matplotlib

    matplotlib.use("Agg")  # Use non-interactive backend
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd

    os.makedirs(output_dir, exist_ok=True)

    if graft_diameter_mm <= 0 or blood_flow_rate_ml_s <= 0:
        return "Error: graft_diameter_mm and blood_flow_rate_ml_s must be positive."

    def wall_shear_stress_dyne_cm2(diameter_mm, flow_ml_s, viscosity_cp):
        mu_poise = viscosity_cp * 0.01
        radius_cm = diameter_mm / 20.0
        return 4 * mu_poise * flow_ml_s / (np.pi * radius_cm**3)

    def sa_to_volume_ratio_per_cm(diameter_mm):
        diameter_cm = diameter_mm / 10.0
        return 4 / diameter_cm

    tau = wall_shear_stress_dyne_cm2(graft_diameter_mm, blood_flow_rate_ml_s, blood_viscosity_cp)
    sa_v = sa_to_volume_ratio_per_cm(graft_diameter_mm)

    ref_diameter = max(diameter_sweep_range)
    tau_ref = wall_shear_stress_dyne_cm2(ref_diameter, blood_flow_rate_ml_s, blood_viscosity_cp)
    sa_v_ref = sa_to_volume_ratio_per_cm(ref_diameter)

    risk_score = (tau / tau_ref) + (sa_v / sa_v_ref)
    if surface_hydrophobicity_score is not None:
        risk_score *= 1 + surface_hydrophobicity_score

    if risk_score < 3:
        risk_category = "Low"
    elif risk_score < 6:
        risk_category = "Moderate"
    else:
        risk_category = "High"

    log = []
    log.append(f"Graft Diameter-Dependent Thrombosis Risk Assessment - {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    log.append(
        "NOTE: This is a simplified heuristic for exploratory design comparisons, not a "
        "validated clinical predictor of thrombosis risk.\n"
    )
    log.append(f"Graft diameter: {graft_diameter_mm} mm")
    log.append(f"Blood flow rate: {blood_flow_rate_ml_s} mL/s")
    log.append(f"Blood viscosity: {blood_viscosity_cp} cP")
    log.append(f"Wall shear stress: {tau:.2f} dyne/cm^2")
    log.append(f"Surface-area-to-volume ratio: {sa_v:.3f} cm^-1")
    log.append(
        f"Reference values at {ref_diameter} mm (same flow rate): "
        f"shear stress = {tau_ref:.2f} dyne/cm^2, SA/V = {sa_v_ref:.3f} cm^-1"
    )
    if surface_hydrophobicity_score is not None:
        log.append(f"Surface hydrophobicity modifier applied: score = {surface_hydrophobicity_score}")
    log.append(f"Heuristic patency risk score: {risk_score:.2f} (relative to reference diameter = 2.0)")
    log.append(f"Risk category: {risk_category}")

    diam_sweep = np.linspace(diameter_sweep_range[0], diameter_sweep_range[1], 50)
    tau_sweep = wall_shear_stress_dyne_cm2(diam_sweep, blood_flow_rate_ml_s, blood_viscosity_cp)
    sa_v_sweep = sa_to_volume_ratio_per_cm(diam_sweep)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    fig, ax1 = plt.subplots(figsize=(9, 6))
    color1 = "tab:red"
    ax1.set_xlabel("Graft Diameter (mm)")
    ax1.set_ylabel("Wall Shear Stress (dyne/cm^2)", color=color1)
    ax1.plot(diam_sweep, tau_sweep, color=color1, label="Wall shear stress")
    ax1.scatter([graft_diameter_mm], [tau], color=color1, zorder=5)
    ax1.tick_params(axis="y", labelcolor=color1)

    ax2 = ax1.twinx()
    color2 = "tab:blue"
    ax2.set_ylabel("Surface-Area-to-Volume Ratio (cm^-1)", color=color2)
    ax2.plot(diam_sweep, sa_v_sweep, color=color2, label="SA/V ratio")
    ax2.scatter([graft_diameter_mm], [sa_v], color=color2, zorder=5)
    ax2.tick_params(axis="y", labelcolor=color2)

    plt.title("Diameter-Dependent Shear Stress and Surface-Area-to-Volume Ratio")
    fig.tight_layout()
    plot_path = os.path.join(output_dir, f"thrombosis_risk_scaling_{timestamp}.png")
    plt.savefig(plot_path, dpi=300, bbox_inches="tight")
    plt.close()

    sweep_df = pd.DataFrame(
        {
            "Diameter_mm": diam_sweep,
            "Wall_Shear_Stress_dyne_cm2": tau_sweep,
            "SA_to_V_ratio_per_cm": sa_v_sweep,
        }
    )
    csv_path = os.path.join(output_dir, f"thrombosis_risk_scaling_{timestamp}.csv")
    sweep_df.to_csv(csv_path, index=False)

    log.append("\nFiles Generated:")
    log.append(f"  - Plot: {plot_path}")
    log.append(f"  - Sweep data: {csv_path}")

    return "\n".join(log)


def generate_multilayer_graft_design_report(
    layer_specs,
    target_compliance_range=(4.4, 5.9),
    target_safety_factor=3.0,
    simulation_weeks=52,
    output_dir="./",
):
    """Generates a structured multilayer biodegradable vascular graft design report.

    Projects per-layer mass retention over time (first-order decay from each
    layer's target half-life), combines thickness-weighted layer
    contributions into an overall mechanical-support timeline, and flags the
    time at which the combined support falls below the critical level
    implied by the target safety factor (i.e. safety factor < 1.0, assuming
    the graft is designed to have exactly target_safety_factor at
    implantation). This is a simplified engineering projection for
    exploratory design comparisons, not validated for clinical or
    regulatory use.

    Parameters
    ----------
    layer_specs : list of dict
        One dict per layer, each with keys: "name" (str), "target_half_life_weeks"
        (float), "thickness_um" (float), and optionally "materials" (str,
        free-text description of the layer's material composition)
    target_compliance_range : tuple of float, optional
        Target compliance range for the overall graft, in %/100mmHg,
        reported in the design summary only (not computed from data in this
        function; see analyze_vascular_graft_mechanical_compliance)
        (default: (4.4, 5.9), spanning saphenous vein (~4.4) to native
        femoral artery (~5.9))
    target_safety_factor : float, optional
        Minimum acceptable combined mechanical-support safety factor,
        assumed to equal the safety factor at implantation (t=0) (default: 3.0)
    simulation_weeks : int, optional
        Number of weeks to project mass retention / mechanical support over
        (default: 52)
    output_dir : str, optional
        Directory to save output files (default: "./")

    Returns
    -------
    str
        Research log summarizing the design, the projected time at which
        combined mechanical support drops below the safety factor threshold
        (if any within simulation_weeks), and saved output file paths

    """
    import os
    from datetime import datetime

    import matplotlib

    matplotlib.use("Agg")  # Use non-interactive backend
    import matplotlib.pyplot as plt
    import numpy as np

    os.makedirs(output_dir, exist_ok=True)

    if not isinstance(layer_specs, list) or len(layer_specs) == 0:
        return "Error: layer_specs must be a non-empty list of layer dicts."

    for i, layer in enumerate(layer_specs):
        for required_key in ("name", "target_half_life_weeks", "thickness_um"):
            if required_key not in layer:
                return f"Error: layer_specs[{i}] is missing required key '{required_key}'."

    total_thickness = sum(layer["thickness_um"] for layer in layer_specs)
    if total_thickness <= 0:
        return "Error: total layer thickness must be positive."

    t = np.linspace(0, simulation_weeks, 200)
    combined_support = np.zeros_like(t)
    per_layer_curves = {}

    for layer in layer_specs:
        k = np.log(2) / layer["target_half_life_weeks"]
        mass_t = 100 * np.exp(-k * t)
        weight = layer["thickness_um"] / total_thickness
        combined_support += weight * mass_t
        per_layer_curves[layer["name"]] = mass_t

    safety_factor_t = target_safety_factor * combined_support / 100

    below_critical_idx = np.where(safety_factor_t < 1.0)[0]
    time_below_critical = t[below_critical_idx[0]] if len(below_critical_idx) > 0 else None

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    plt.figure(figsize=(10, 6))
    for name, mass_t in per_layer_curves.items():
        plt.plot(t, mass_t, "--", alpha=0.6, label=f"{name} (mass retention)")
    plt.plot(t, combined_support, "k-", linewidth=2, label="Combined mechanical support")
    plt.axhline(100 / target_safety_factor, color="red", linestyle=":", label="Critical support level (SF=1.0)")
    plt.xlabel("Time (weeks)")
    plt.ylabel("Mass Retention / Mechanical Support (%)")
    plt.title("Multilayer Graft: Projected Mass Retention and Combined Mechanical Support")
    plt.legend()
    plt.grid(True, linestyle="--", alpha=0.7)
    plot_path = os.path.join(output_dir, f"multilayer_design_{timestamp}.png")
    plt.savefig(plot_path, dpi=300, bbox_inches="tight")
    plt.close()

    report_lines = []
    report_lines.append("=" * 80)
    report_lines.append("MULTILAYER BIODEGRADABLE VASCULAR GRAFT DESIGN REPORT")
    report_lines.append("=" * 80)
    report_lines.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    report_lines.append(
        "\nNOTE: Simplified engineering projection for exploratory design comparisons. "
        "Not validated for clinical or regulatory use.\n"
    )
    report_lines.append("LAYER SPECIFICATIONS")
    report_lines.append("-" * 80)
    for layer in layer_specs:
        weight_pct = 100 * layer["thickness_um"] / total_thickness
        report_lines.append(f"- {layer['name']}")
        report_lines.append(f"    Target half-life: {layer['target_half_life_weeks']} weeks")
        report_lines.append(f"    Thickness: {layer['thickness_um']} um ({weight_pct:.1f}% of total thickness)")
        if "materials" in layer:
            report_lines.append(f"    Materials: {layer['materials']}")
    report_lines.append("")
    report_lines.append("DESIGN TARGETS")
    report_lines.append("-" * 80)
    report_lines.append(f"Target compliance range: {target_compliance_range[0]}-{target_compliance_range[1]} %/100mmHg")
    report_lines.append(f"Target safety factor at implantation: {target_safety_factor}")
    report_lines.append("")
    report_lines.append("PROJECTED MECHANICAL SUPPORT TIMELINE")
    report_lines.append("-" * 80)
    report_lines.append(f"Simulation horizon: {simulation_weeks} weeks")
    report_lines.append(f"Combined support at t=0: {combined_support[0]:.1f}% (safety factor = {target_safety_factor})")
    report_lines.append(f"Combined support at t={simulation_weeks}w: {combined_support[-1]:.1f}%")
    if time_below_critical is not None:
        report_lines.append(
            f"FLAG: Combined support drops below the critical level (safety factor < 1.0) at "
            f"t = {time_below_critical:.1f} weeks"
        )
    else:
        report_lines.append("Combined support stays above the critical level (safety factor >= 1.0) throughout")

    txt_path = os.path.join(output_dir, f"multilayer_design_report_{timestamp}.txt")
    with open(txt_path, "w") as f:
        f.write("\n".join(report_lines))

    log = []
    log.append(f"Multilayer Graft Design Report Generated - {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    log.append(f"Number of layers: {len(layer_specs)}")
    for layer in layer_specs:
        log.append(f"  - {layer['name']}: target half-life {layer['target_half_life_weeks']} weeks, "
                    f"thickness {layer['thickness_um']} um")
    log.append(f"Combined support at t=0: {combined_support[0]:.1f}%")
    log.append(f"Combined support at t={simulation_weeks}w: {combined_support[-1]:.1f}%")
    if time_below_critical is not None:
        log.append(f"FLAG: Drops below critical safety factor (1.0) at t = {time_below_critical:.1f} weeks")
    log.append("\nFiles Generated:")
    log.append(f"  - Design report: {txt_path}")
    log.append(f"  - Plot: {plot_path}")

    return "\n".join(log)
