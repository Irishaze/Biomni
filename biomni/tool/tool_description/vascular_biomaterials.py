description = [
    {
        "description": "Fits scaffold/coating mass-retention decay over time, for one or more layers, and derives decay constant, half-life, and R^2 per layer.",
        "name": "analyze_scaffold_degradation_kinetics",
        "optional_parameters": [
            {
                "default": "first_order",
                "description": "Degradation model to fit: 'first_order' (exponential decay) or 'autocatalytic' (logistic decay approximation of autocatalytic hydrolysis)",
                "name": "degradation_model",
                "type": "str",
            },
            {
                "default": "./",
                "description": "Directory to save output files",
                "name": "output_dir",
                "type": "str",
            },
        ],
        "required_parameters": [
            {
                "default": None,
                "description": "Time points at which mass retention was measured, in weeks",
                "name": "time_points",
                "type": "List[float] or numpy.ndarray",
            },
            {
                "default": None,
                "description": "Measured mass retention (%) at each time point. Pass a 1D list/array for a single layer, or a dict mapping layer name to its own list/array for multiple layers",
                "name": "mass_retention_data",
                "type": "List[float] or numpy.ndarray or dict",
            },
        ],
    },
    {
        "description": "Computes vascular graft compliance, distensibility, and mechanical safety factor from pressure-diameter test data, compared against target values.",
        "name": "analyze_vascular_graft_mechanical_compliance",
        "optional_parameters": [
            {
                "default": 120,
                "description": "Physiological systolic pressure used for the compliance calculation, in mmHg",
                "name": "systolic_pressure",
                "type": "float",
            },
            {
                "default": 80,
                "description": "Physiological diastolic pressure used for the compliance calculation, in mmHg",
                "name": "diastolic_pressure",
                "type": "float",
            },
            {
                "default": None,
                "description": "Measured graft burst pressure, in mmHg. If provided, used to compute a mechanical safety factor",
                "name": "burst_pressure",
                "type": "float",
            },
            {
                "default": "(4.4, 5.9)",
                "description": "Target compliance range for a compliance-matched graft, in %/100mmHg (spans saphenous vein ~4.4 to native femoral artery ~5.9)",
                "name": "target_compliance_range",
                "type": "tuple of float",
            },
            {
                "default": 3.0,
                "description": "Minimum acceptable ratio of burst pressure to systolic pressure",
                "name": "target_safety_factor",
                "type": "float",
            },
            {
                "default": "./",
                "description": "Directory to save output files",
                "name": "output_dir",
                "type": "str",
            },
        ],
        "required_parameters": [
            {
                "default": None,
                "description": "Applied intraluminal pressure at each measurement, in mmHg",
                "name": "pressure_data",
                "type": "List[float] or numpy.ndarray",
            },
            {
                "default": None,
                "description": "Measured graft diameter at each corresponding pressure, in mm",
                "name": "diameter_data",
                "type": "List[float] or numpy.ndarray",
            },
        ],
    },
    {
        "description": "Estimates diameter-dependent wall shear stress and a heuristic thrombosis-risk score for a cylindrical vascular graft.",
        "name": "assess_graft_diameter_thrombosis_risk",
        "optional_parameters": [
            {
                "default": 3.5,
                "description": "Blood dynamic viscosity, in centipoise",
                "name": "blood_viscosity_cp",
                "type": "float",
            },
            {
                "default": None,
                "description": "Optional surface hydrophobicity score from 0 (hydrophilic) to 1 (hydrophobic), applied as a linear multiplier on the risk score",
                "name": "surface_hydrophobicity_score",
                "type": "float",
            },
            {
                "default": "(2, 10)",
                "description": "Diameter range (mm) to sweep for the comparison plot",
                "name": "diameter_sweep_range",
                "type": "tuple of float",
            },
            {
                "default": "./",
                "description": "Directory to save output files",
                "name": "output_dir",
                "type": "str",
            },
        ],
        "required_parameters": [
            {
                "default": None,
                "description": "Inner diameter of the graft, in mm",
                "name": "graft_diameter_mm",
                "type": "float",
            },
            {
                "default": None,
                "description": "Volumetric blood flow rate through the graft, in mL/s",
                "name": "blood_flow_rate_ml_s",
                "type": "float",
            },
        ],
    },
    {
        "description": "Generates a structured multilayer biodegradable vascular graft design report, projecting per-layer mass retention and combined mechanical-support timeline.",
        "name": "generate_multilayer_graft_design_report",
        "optional_parameters": [
            {
                "default": "(4.4, 5.9)",
                "description": "Target compliance range for the overall graft, in %/100mmHg, reported in the design summary only (spans saphenous vein ~4.4 to native femoral artery ~5.9)",
                "name": "target_compliance_range",
                "type": "tuple of float",
            },
            {
                "default": 3.0,
                "description": "Minimum acceptable combined mechanical-support safety factor, assumed to equal the safety factor at implantation (t=0)",
                "name": "target_safety_factor",
                "type": "float",
            },
            {
                "default": 52,
                "description": "Number of weeks to project mass retention / mechanical support over",
                "name": "simulation_weeks",
                "type": "int",
            },
            {
                "default": "./",
                "description": "Directory to save output files",
                "name": "output_dir",
                "type": "str",
            },
        ],
        "required_parameters": [
            {
                "default": None,
                "description": "One dict per layer, each with keys 'name' (str), 'target_half_life_weeks' (float), 'thickness_um' (float), and optionally 'materials' (str)",
                "name": "layer_specs",
                "type": "List[dict]",
            },
        ],
    },
]
