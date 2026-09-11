"""
tree_stages.py  v3.0
====================
Picea abies development stages — calibrated with Norwegian field data.

MAJOR UPDATE v3.0: Integrated field measurements from Shore (2025):
  "Identifying patterns of growth in Picea abies from soil moisture
  conditions in the boreal forest of Norway", Wageningen University MSc Thesis.

KEY IMPROVEMENTS FROM PAPER:
  1. Norwegian boreal soil VG parameters from Shore (2025) Table 3.1
     (replaces generic loam/sand/clay — mean values from 20 Norwegian plots)
  2. Soil depth constraint: Norwegian soils 7-80 cm, median 9-30 cm
     → taproot lmax reduced 55→30 cm, sinker lmax 70→35 cm, zmin = -80 cm
  3. REW < 0.40 = water stress threshold (Shore 2025, Granier et al. 1999)
  4. O-layer (organic horizon) explicitly parameterised
     Shore (2025): O-layer holds 26-47% of total root length in forest soils
  5. MPdiff concept: rate of change of matric potential is strongest
     soil-based growth predictor (Shore 2025, correlation -0.39 with dGROpc)
  6. Norwegian monitored plots = 35-49 years old → Pole_OMT / Pole_VT stages

Sources:
  Shore (2025)              — Norwegian Picea abies soil moisture and growth
  Kalliokoski et al. (2008) — Finnish Picea abies field measurements
  Finér et al. (2007)       — Fine root turnover (rlt ≈ 200 days)
  Stokes et al. (2009)      — Plate-sinker root architecture
"""

# ======================================================================
# OBSERVED TREE STAGES
# ======================================================================

TREE_STAGES = [
    {
        "stage"        : "Sapling",
        "age_mean"     : 15,
        "age_std"      : 0.8,
        "d0_mean"      : 7.9,
        "d0_std"       : 0.41,
        "dbh_mean"     : 4.7,
        "dbh_std"      : 0.33,
        "height_mean"  : 4.1,
        "height_std"   : 0.53,
        "site_type"    : "Myrtillus MT",
        "sim_time"     : 15 * 365,
    },
    {
        "stage"        : "Pole_MT",
        "age_mean"     : 28,
        "age_std"      : 3.5,
        "d0_mean"      : 11.7,
        "d0_std"       : 2.07,
        "dbh_mean"     : 7.8,
        "dbh_std"      : 1.65,
        "height_mean"  : 10.2,
        "height_std"   : 2.37,
        "site_type"    : "Myrtillus MT",
        "sim_time"     : 28 * 365,
    },
    {
        "stage"        : "Pole_OMT",
        # Shore (2025) reference plots: 35-49 years, richest sites
        "age_mean"     : 34,
        "age_std"      : 3.3,
        "d0_mean"      : 18.2,
        "d0_std"       : 3.06,
        "dbh_mean"     : 13.0,
        "dbh_std"      : 2.45,
        "height_mean"  : 13.5,
        "height_std"   : 1.74,
        "site_type"    : "Oxalis-Myrtillus OMT",
        "sim_time"     : 34 * 365,
    },
    {
        "stage"        : "Pole_VT",
        # Shore (2025) dry plots: 35-49 years, continental dry sites
        "age_mean"     : 41,
        "age_std"      : 3.7,
        "d0_mean"      : 14.7,
        "d0_std"       : 2.64,
        "dbh_mean"     : 10.2,
        "dbh_std"      : 2.11,
        "height_mean"  : 9.5,
        "height_std"   : 1.25,
        "site_type"    : "Vaccinium VT",
        "sim_time"     : 41 * 365,
    },
    {
        "stage"        : "Mature",
        "age_mean"     : 55,
        "age_std"      : 0.6,
        "d0_mean"      : 29.2,
        "d0_std"       : 1.40,
        "dbh_mean"     : 21.8,
        "dbh_std"      : 1.12,
        "height_mean"  : 21.3,
        "height_std"   : 2.14,
        "site_type"    : "Myrtillus MT",
        "sim_time"     : 55 * 365,
    },
]


# ======================================================================
# NORWEGIAN BOREAL FOREST SOIL LAYERS (v3.0 NEW)
# Source: Shore (2025), Table 3.1 — MvG parameters from 20 Norwegian plots.
#
# Three-layer Norwegian boreal profile:
#   O-layer  :  0 to  -8 cm — organic horizon (humus/mor)
#   Mineral  : -8 to -35 cm — mineral topsoil + subsoil
#   Till     :-35 to -80 cm — till/weathered bedrock
#
# MvG parameter means from Table 3.1:
#   Topsoil: θsat=0.72, θres=0.03, n=1.27, α=15.0  (organic-dominated)
#   Subsoil: θsat=0.57, θres=0.04, n=1.30, α=8.0   (mineral sand)
#
# KEY: Norwegian soils are SHALLOW. Shore (2025) Fig 3.6:
#   Median total soil depth per plot: 9-30 cm
#   Maximum recorded depth: 80 cm
# ======================================================================

NORWEGIAN_SOIL_LAYERS = [
    {
        "label"    : "OLayer",
        "z_top"    :   0.0,
        "z_bottom" :  -8.0,    # O-layer: 5-15 cm in Norwegian boreal
        # From Shore (2025) Table 3.1 topsoil averages
        "theta_r"  : 0.03,
        "theta_s"  : 0.74,     # High porosity of organic material
        "alpha"    : 15.0,     # Fast drainage (ranges 2.5-28.4 across plots)
        "n"        : 1.27,
        "Ks"       : 120.0,
        "z_nappe"  : -55.0,
        "theta_fc" : 0.60,     # Field capacity (pF 1.8 ≈ 6 kPa)
        "theta_pwp": 0.10,     # Permanent wilting point (pF 4.2 ≈ 1500 kPa)
    },
    {
        "label"    : "Mineral",
        "z_top"    :  -8.0,
        "z_bottom" : -35.0,    # Norwegian median total depth ~20-30 cm
        # From Shore (2025) Table 3.1 subsoil averages
        "theta_r"  : 0.04,
        "theta_s"  : 0.57,
        "alpha"    :  8.0,     # Slower drainage (ranges 1.4-25.5)
        "n"        : 1.30,
        "Ks"       :  55.0,
        "z_nappe"  : -90.0,
        "theta_fc" : 0.44,
        "theta_pwp": 0.08,
    },
    {
        "label"    : "Till",
        "z_top"    : -35.0,
        "z_bottom" : -80.0,    # Maximum Norwegian soil depth
        # Till/rock fragments — very few roots penetrate here
        "theta_r"  : 0.02,
        "theta_s"  : 0.18,
        "alpha"    :  1.0,
        "n"        : 1.15,
        "Ks"       :  1.5,
        "z_nappe"  : -120.0,
        "theta_fc" : 0.14,
        "theta_pwp": 0.04,
    },
]

# REW threshold for water stress (Shore 2025, Granier et al. 1999)
REW_STRESS_THRESHOLD = 0.40


# ======================================================================
# SITE MODIFIERS
# ======================================================================

SITE_MODIFIERS = {
    "Myrtillus MT": {
        "z_nappe_offset"  :   0.0,
        "absorption_rate" : 10e-6,
        "n_pockets_base"  : 60,
        "rew_stress_freq" : "low",
        "description"     : "Average moisture, moderate fertility (Norwegian reference)",
    },
    "Oxalis-Myrtillus OMT": {
        "z_nappe_offset"  : +25.0,
        "absorption_rate" : 12e-6,
        "n_pockets_base"  : 80,
        "rew_stress_freq" : "very_low",
        "description"     : "Rich moist Norwegian forest (OMT type)",
    },
    "Vaccinium VT": {
        "z_nappe_offset"  : -35.0,
        "absorption_rate" :  8e-6,
        "n_pockets_base"  : 40,
        "rew_stress_freq" : "high",
        "description"     : "Poor dry Norwegian forest (VT, continental)",
    },
}


# ======================================================================
# D0 → ROOT SYSTEM SCALING (with Norwegian soil depth constraint)
# ======================================================================

def d0_to_root_extent(d0_cm: float, soil_depth_cm: float = 25.0) -> dict:
    """
    Estimate root dimensions from D0. Norwegian soil depth constrains depth.

    Parameters
    ----------
    d0_cm        : root collar diameter (cm)
    soil_depth_cm: actual Norwegian soil depth (Shore 2025: median 9-30 cm)
    """
    lateral    = 18.0 * d0_cm
    tap_allom  =  3.5 * d0_cm
    sink_allom =  4.0 * d0_cm
    # Soil constraint: roots stopped by bedrock/till
    tap_depth  = min(tap_allom,  soil_depth_cm * 1.1)
    sink_depth = min(sink_allom, soil_depth_cm * 0.9)
    return {
        "lateral_spread_cm" : round(lateral,    1),
        "taproot_depth_cm"  : round(tap_depth,  1),
        "sinker_depth_cm"   : round(sink_depth, 1),
        "max_depth_cm"      : round(tap_depth + sink_depth, 1),
        "lmax_lateral"      : round(lateral,    1),
        "lmax_sinker"       : round(sink_depth, 1),
        "soil_constrained"  : tap_allom > soil_depth_cm * 1.1,
    }


# ======================================================================
# MAIN CONFIG BUILDER
# ======================================================================

def get_sim_config_for_stage(stage_name: str,
                              soil_depth_cm: float = 25.0) -> dict:
    """
    Return SimConfig overrides for a development stage.
    v3.0: Norwegian soil depth replaces -200 cm domain.
    """
    stage = next((s for s in TREE_STAGES if s["stage"] == stage_name), None)
    if stage is None:
        raise KeyError(f"Stage '{stage_name}' not found. "
                       f"Available: {[s['stage'] for s in TREE_STAGES]}")

    site    = SITE_MODIFIERS[stage["site_type"]]
    extents = d0_to_root_extent(stage["d0_mean"], soil_depth_cm)

    margin   = extents["lateral_spread_cm"] * 1.25
    sim_days = stage["sim_time"]
    dt       = max(50, sim_days // 200)

    domain_area_units = (2 * margin / 100) ** 2
    n_pockets = int(site["n_pockets_base"] * (1.0 + domain_area_units * 0.3))
    n_pockets = max(40, min(n_pockets, 250))

    return {
        "sim_time"        : sim_days,
        "dt"              : dt,
        "absorption_rate" : site["absorption_rate"],
        "n_pockets"       : n_pockets,
        "xmin"            : -margin,
        "xmax"            :  margin,
        "ymin"            : -margin,
        "ymax"            :  margin,
        # v3.0 FIX: -80 cm = Norwegian max soil depth (Shore 2025, Fig 3.6)
        # Previous v1/v2 used -200 cm which is WRONG for Norwegian forest soils
        "zmin"            : -80.0,
        "zmax"            :  5.0,
        "depth_limit"     : -75.0,
        "output_dir"      : f"results/LayeredSoil_{stage_name}",
    }


# ======================================================================
# REW CALCULATION (v3.0 NEW)
# Source: Shore (2025), Chapter 3.4
# REW = (theta - theta_PWP) / (theta_FC - theta_PWP)
# REW < 0.40 → water-limiting for Picea abies
# ======================================================================

def calculate_rew(theta: float, layer: dict) -> float:
    """Relative Extractable Water. REW < 0.4 = water stress (Shore 2025)."""
    theta_fc  = layer.get("theta_fc",  layer["theta_s"] * 0.85)
    theta_pwp = layer.get("theta_pwp", layer["theta_r"] * 3.0)
    theta_pwp = max(theta_pwp, layer["theta_r"] + 0.01)
    available = max(0.0, theta_fc - theta_pwp)
    if available < 1e-6:
        return 0.0
    return max(0.0, min(1.0, (theta - theta_pwp) / available))


# ======================================================================
# FEATURE EXTRACTOR
# ======================================================================

def extract_features(stage_name: str, results: dict) -> dict:
    import numpy as np
    stage   = next(s for s in TREE_STAGES if s["stage"] == stage_name)
    site    = SITE_MODIFIERS[stage["site_type"]]
    extents = d0_to_root_extent(stage["d0_mean"])
    tips    = results["final_tips"]

    depths = [z for _, _, _, z in tips] if tips else [0.0]
    radii  = [(x**2+y**2)**0.5 for _, x, y, _ in tips] if tips else [0.0]

    # Norwegian layers
    counts = {"OLayer": 0, "Mineral": 0, "Till": 0}
    for _, _, _, z in tips:
        if   z >= -8.0:  counts["OLayer"]  += 1
        elif z >= -35.0: counts["Mineral"] += 1
        else:            counts["Till"]    += 1

    n = max(len(tips), 1)
    return {
        "age_years"              : stage["age_mean"],
        "d0_cm"                  : stage["d0_mean"],
        "dbh_cm"                 : stage["dbh_mean"],
        "height_m"               : stage["height_mean"],
        "site_type"              : stage["site_type"],
        "stage_name"             : stage_name,
        "expected_lateral_cm"    : extents["lateral_spread_cm"],
        "soil_constrained"       : extents["soil_constrained"],
        "z_nappe_offset"         : site["z_nappe_offset"],
        "rew_stress_freq"        : site["rew_stress_freq"],
        "n_tips_total"           : len(tips),
        "max_depth_cm"           : round(abs(min(depths)), 2),
        "mean_depth_cm"          : round(abs(float(np.mean(depths))), 2),
        "max_lateral_cm"         : round(max(radii), 2) if radii else 0.0,
        "pct_tips_olayer"        : round(100 * counts["OLayer"]  / n, 1),
        "pct_tips_mineral"       : round(100 * counts["Mineral"] / n, 1),
        "pct_tips_till"          : round(100 * counts["Till"]    / n, 1),
        "lateral_allometry_ratio": round(max(radii)/extents["lateral_spread_cm"], 3)
                                   if radii else 0.0,
    }


# ======================================================================
# ROOT DEPTH DISTRIBUTION TABLE (Kalliokoski 2008)
# ======================================================================

ROOT_DEPTH_DISTRIBUTION = {
    "pasture": {
        "sapling"     : {"0_10cm": 0.60, "10_50cm": 0.24, "50_100cm": 0.15,
                         "100cm+": 0.01, "total_root_length_m_m2": 847,
                         "live_root_mass_g_m2": 108, "dead_root_mass_g_m2": 65,
                         "specific_root_length_m_g": 7.8},
        "pole_stage"  : {"0_10cm": 0.36, "10_50cm": 0.40, "50_100cm": 0.17,
                         "100cm+": 0.06, "total_root_length_m_m2": 1797,
                         "live_root_mass_g_m2": 197, "dead_root_mass_g_m2": 105,
                         "specific_root_length_m_g": 9.1},
        "young_timber": {"0_10cm": 0.45, "10_50cm": 0.35, "50_100cm": 0.04,
                         "100cm+": 0.17, "total_root_length_m_m2": 1762,
                         "live_root_mass_g_m2": 238, "dead_root_mass_g_m2": 168,
                         "specific_root_length_m_g": 7.4},
    },
    "forest": {
        "sapling"     : {"O_layer": 0.26, "0_10cm": 0.34, "10_50cm": 0.09,
                         "50_100cm": 0.16, "100cm+": 0.16,
                         "total_root_length_m_m2": 709, "live_root_mass_g_m2": 154,
                         "dead_root_mass_g_m2": 488, "specific_root_length_m_g": 4.6},
        "pole_stage"  : {
            # Shore (2025): Norwegian pole-stage plots
            # O-layer fraction = 47% (highest recorded)
            "O_layer": 0.47, "0_10cm": 0.24, "10_50cm": 0.17,
            "50_100cm": 0.12, "100cm+": 0.01,
            "total_root_length_m_m2": 1348, "live_root_mass_g_m2": 227,
            "dead_root_mass_g_m2": 234, "specific_root_length_m_g": 5.9},
        "young_timber": {"O_layer": 0.40, "0_10cm": 0.15, "10_50cm": 0.38,
                         "50_100cm": 0.07, "100cm+": 0.00,
                         "total_root_length_m_m2": 757, "live_root_mass_g_m2": 181,
                         "dead_root_mass_g_m2": 394, "specific_root_length_m_g": 4.2},
    },
}

STAGE_TO_TABLE_STAGE = {
    "Sapling" : "sapling",
    "Pole_MT" : "pole_stage",
    "Pole_OMT": "pole_stage",
    "Pole_VT" : "pole_stage",
    "Mature"  : "young_timber",
}
SITE_TO_TABLE_SOIL = {
    "Myrtillus MT"         : "forest",
    "Oxalis-Myrtillus OMT" : "forest",
    "Vaccinium VT"         : "forest",
}


# ======================================================================
# VALIDATION
# ======================================================================

def validate_depth_distribution(final_tips: list,
                                 stage_name: str,
                                 site_type: str) -> dict:
    """Compare simulated depth distribution against field data."""
    if not final_tips:
        print("[WARNING] No tips to validate.")
        return {"well_calibrated": False}

    table_stage  = STAGE_TO_TABLE_STAGE.get(stage_name, "sapling")
    table_soil   = SITE_TO_TABLE_SOIL.get(site_type,   "forest")
    observed_raw = ROOT_DEPTH_DISTRIBUTION[table_soil][table_stage]

    mineral_keys  = ["0_10cm", "10_50cm", "50_100cm", "100cm+"]
    mineral_total = sum(observed_raw.get(k, 0.0) for k in mineral_keys)
    observed = {k: round(observed_raw.get(k, 0.0) / max(mineral_total, 1e-9), 3)
                for k in mineral_keys}

    total = len(final_tips)
    sim_counts = {
        "0_10cm"   : sum(1 for _, _, _, z in final_tips if -10.0 <= z <=   0.0),
        "10_50cm"  : sum(1 for _, _, _, z in final_tips if -50.0 <= z < -10.0),
        "50_100cm" : sum(1 for _, _, _, z in final_tips if -100.0<= z < -50.0),
        "100cm+"   : sum(1 for _, _, _, z in final_tips if          z < -100.0),
    }
    nor_counts = {
        "O-layer(0-8cm)"  : sum(1 for _, _, _, z in final_tips if  -8.0 <= z <= 0.0),
        "Mineral(8-35cm)" : sum(1 for _, _, _, z in final_tips if -35.0 <= z < -8.0),
        "Till(35-80cm)"   : sum(1 for _, _, _, z in final_tips if -80.0 <= z < -35.0),
    }
    sim_fracs  = {k: round(v / total, 3) for k, v in sim_counts.items()}
    errors     = {k: round(abs(sim_fracs[k] - observed[k]), 3) for k in mineral_keys}
    mean_error = round(sum(errors.values()) / len(errors), 3)

    def flag(e): return "PASS " if e<=0.08 else ("WARN " if e<=0.15 else "FAIL ⚠")

    print(f"\n{'='*68}")
    print(f"  Depth Distribution Validation — {stage_name} | {site_type}")
    print(f"  {'Zone':<14} {'Simulated':>10} {'Observed':>10} {'Error':>8}  Status")
    print(f"  {'-'*52}")
    for k in mineral_keys:
        print(f"  {k:<14} {sim_fracs[k]:>10.3f} {observed[k]:>10.3f} "
              f"{errors[k]:>8.3f}  {flag(errors[k])}")
    print(f"  {'-'*52}")
    print(f"  {'Mean error':<14} {'':>20} {mean_error:>8.3f}  "
          f"{'CALIBRATED' if mean_error<0.10 else 'NEEDS TUNING'}")

    o_obs = observed_raw.get("O_layer", 0.0) * 100
    o_sim = 100 * nor_counts["O-layer(0-8cm)"] / total
    print(f"\n  Norwegian soil layers (Shore 2025):")
    print(f"    O-layer: {o_sim:.0f}% simulated vs {o_obs:.0f}% observed "
          f"({'OK' if abs(o_sim-o_obs)<10 else 'ADJUST subType 5 rate'})")
    for lbl, cnt in nor_counts.items():
        print(f"    {lbl:<18}: {cnt:4d} tips ({100*cnt/total:.1f}%)")

    print(f"\n  Field data: RL={observed_raw['total_root_length_m_m2']} m/m²  "
          f"LRM={observed_raw['live_root_mass_g_m2']} g/m²  "
          f"SRL={observed_raw['specific_root_length_m_g']} m/g")
    print(f"{'='*68}\n")

    return {
        "simulated"       : sim_fracs,
        "observed"        : observed,
        "errors"          : errors,
        "mean_error"      : mean_error,
        "well_calibrated" : mean_error < 0.10,
        "norwegian_layers": {k: v/total for k, v in nor_counts.items()},
    }


# ======================================================================
# QUICK SUMMARY
# ======================================================================
if __name__ == "__main__":
    print(f"\n{'='*74}")
    print(f"  Picea abies Stages — Shore (2025) Norwegian context")
    print(f"  Norwegian monitored plots: 35-49 years (Pole_OMT / Pole_VT)")
    print(f"  {'Stage':<12} {'Age':>5} {'D0':>6} {'SimDays':>8} "
          f"{'Lateral':>9} {'Depth(soil)':>12}")
    print(f"  {'-'*60}")
    for s in TREE_STAGES:
        ext = d0_to_root_extent(s["d0_mean"], soil_depth_cm=25.0)
        print(f"  {s['stage']:<12} {s['age_mean']:>5} "
              f"{s['d0_mean']:>6.1f} {s['sim_time']:>8d} "
              f"{ext['lateral_spread_cm']:>8.0f}cm "
              f"{ext['max_depth_cm']:>10.0f}cm")

    print(f"\n  Norwegian soil layers (Shore 2025, Table 3.1 averages):")
    for L in NORWEGIAN_SOIL_LAYERS:
        print(f"    {L['label']:<10} ({L['z_top']:+5.0f}→{L['z_bottom']:+5.0f} cm)"
              f"  θs={L['theta_s']:.2f}  α={L['alpha']:5.1f}  n={L['n']:.2f}"
              f"  REW_FC={L['theta_fc']:.2f}  REW_PWP={L['theta_pwp']:.2f}")

    print(f"\n  REW stress threshold: {REW_STRESS_THRESHOLD} (Shore 2025)")
    print(f"\n  SimConfig for 'Pole_MT' (soil_depth=25cm):")
    ov = get_sim_config_for_stage("Pole_MT", soil_depth_cm=25.0)
    for k, v in ov.items():
        print(f"    {k:<24} = {v}")