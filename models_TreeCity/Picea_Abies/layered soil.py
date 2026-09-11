"""
root_growth_predictor.py  v4.0
================================
Picea abies root growth — calibrated for Norwegian boreal forest conditions.

MAJOR UPDATE v4.0: Integrated field data from Shore (2025):
  "Identifying patterns of growth in Picea abies from soil moisture
  conditions in the boreal forest of Norway", Wageningen University MSc Thesis.

────────────────────────────────────────────────────────────────────────
KEY CHANGES FROM v3 → v4:

SOIL PROFILE (most impactful change):
  v3: Generic loam/sand/clay to -200 cm depth
  v4: Norwegian boreal O-layer/Mineral/Till to -80 cm
      Parameters from Shore (2025) Table 3.1 (mean of 20 Norwegian plots)
      Shore (2025) Fig 3.6: median soil depth 9-30 cm, maximum 80 cm
      Using -200 cm was completely unrealistic for Norwegian forest soils.

DEPTH CONSTRAINTS:
  zmin:              -200 cm → -80 cm (Norwegian max soil depth)
  depth_limit:       -180 cm → -75 cm
  Taproot lmax (XML):  55 cm → 30 cm (soil constraint)
  Sinker lmax (XML):   70 cm → 35 cm (soil constraint)

REW-BASED HYDROTROPISM (new feature from Shore 2025):
  Shore (2025): REW < 0.40 = water-limiting conditions for Picea abies
  New REW calculation in LayeredSoil.getREW()
  LateralPlagiotropism and GravHydroTropism now check REW stress.
  When REW < 0.4: hydrotropism weight multiplied by 2.0 (roots chase water).
  MPdiff concept: when moisture is increasing (soil rewetting),
  roots grow toward wetter zones (Shore 2025: strongest soil-based predictor).

O-LAYER FINE ROOTS (new subType 5):
  Shore (2025): organic O-layer holds 26-47% of root length in forest soils.
  New XML subType 5 (surface_fine_roots) grows into O-layer.
  SurfaceFineTropism: upward bias + horizontal spread for O-layer colonisation.
  Python setTropism() assigns SurfaceFineTropism to subType 5.

SEASONAL CONTEXT:
  Shore (2025): Norwegian growing season 100-200 days (latitude-dependent).
  Daylength is strongest predictor of growth probability (effect size ~-3).
  VPD is second strongest predictor.
  These affect ABOVE-GROUND growth — root growth is more continuous.

────────────────────────────────────────────────────────────────────────

Norwegian soil profile (Shore 2025, Table 3.1):
    O-layer  :   0 →  -8 cm : Organic horizon (θs=0.74, α=15.0)
    Mineral  :  -8 → -35 cm : Mineral topsoil/subsoil (θs=0.57, α=8.0)
    Till     : -35 → -80 cm : Till/bedrock (θs=0.18, α=1.0)

References:
  Shore (2025) — Norwegian Picea abies, Wageningen University MSc Thesis
  Kalliokoski et al. (2008) — Finnish Picea abies field measurements
  Stokes et al. (2009) — Plate-sinker architecture
  Finér et al. (2007) — Fine root turnover (rlt ≈ 200 days)
  Granier et al. (1999) — REW threshold (0.40)
"""

import os
import sys
import random
import pathlib
import time
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
import pandas as pd

try:
    from tree_stages import (
        get_sim_config_for_stage, validate_depth_distribution,
        TREE_STAGES, SITE_MODIFIERS, NORWEGIAN_SOIL_LAYERS,
        calculate_rew, REW_STRESS_THRESHOLD,
    )
    _STAGES_AVAILABLE = True
except ImportError:
    _STAGES_AVAILABLE = False

    # ── Fallback definitions so this file works without tree_stages.py ──

    REW_STRESS_THRESHOLD = 0.40  # Granier et al. (1999), Shore (2025)

    def calculate_rew(theta: float, layer: dict) -> float:
        """Fallback REW calculation when tree_stages.py is unavailable."""
        theta_fc  = layer.get("theta_fc",  layer["theta_s"] * 0.85)
        theta_pwp = layer.get("theta_pwp", layer["theta_r"] * 3.0)
        theta_pwp = max(theta_pwp, layer["theta_r"] + 0.01)
        available = max(0.0, theta_fc - theta_pwp)
        if available < 1e-6:
            return 0.0
        return max(0.0, min(1.0, (theta - theta_pwp) / available))

    NORWEGIAN_SOIL_LAYERS = [
        {
            "label": "OLayer",  "z_top":   0.0, "z_bottom":  -8.0,
            "theta_r": 0.03, "theta_s": 0.74, "alpha": 15.0, "n": 1.27,
            "Ks": 120.0, "z_nappe": -55.0, "theta_fc": 0.60, "theta_pwp": 0.10,
        },
        {
            "label": "Mineral", "z_top":  -8.0, "z_bottom": -35.0,
            "theta_r": 0.04, "theta_s": 0.57, "alpha":  8.0, "n": 1.30,
            "Ks":  55.0, "z_nappe": -90.0, "theta_fc": 0.44, "theta_pwp": 0.08,
        },
        {
            "label": "Till",    "z_top": -35.0, "z_bottom": -80.0,
            "theta_r": 0.02, "theta_s": 0.18, "alpha":  1.0, "n": 1.15,
            "Ks":   1.5, "z_nappe": -120.0, "theta_fc": 0.14, "theta_pwp": 0.04,
        },
    ]

    # Minimal stubs so the __main__ block runs without tree_stages
    TREE_STAGES    = []
    SITE_MODIFIERS = {}

    def get_sim_config_for_stage(stage_name, soil_depth_cm=25.0):
        raise ImportError("tree_stages.py not found — cannot load stage config.")

    def validate_depth_distribution(tips, stage, site):
        print("[WARNING] tree_stages.py not found — skipping validation.")

import plantbox as pb
import plantbox.visualisation.vtk_plot as vp
import plantbox.functional.van_genuchten as vg


_HERE    = pathlib.Path(__file__).resolve().parent
XML_PATH = str(
    (_HERE / "../../modelparameter_TreeCity/structural/Picea_Abies_simple.xml")
    .resolve()
)


# ======================================================================
# SOIL LAYERS — Norwegian boreal forest profile (Shore 2025)
# Replaces generic loam/sand/clay from v1-v3.
# ======================================================================

LAYERS = NORWEGIAN_SOIL_LAYERS if _STAGES_AVAILABLE else [
    # Fallback if tree_stages.py not available
    {
        "label": "OLayer",   "z_top":   0.0, "z_bottom":  -8.0,
        "theta_r": 0.03, "theta_s": 0.74, "alpha": 15.0, "n": 1.27,
        "Ks": 120.0, "z_nappe": -55.0, "theta_fc": 0.60, "theta_pwp": 0.10,
    },
    {
        "label": "Mineral",  "z_top":  -8.0, "z_bottom": -35.0,
        "theta_r": 0.04, "theta_s": 0.57, "alpha":  8.0, "n": 1.30,
        "Ks":  55.0, "z_nappe": -90.0, "theta_fc": 0.44, "theta_pwp": 0.08,
    },
    {
        "label": "Till",     "z_top": -35.0, "z_bottom": -80.0,
        "theta_r": 0.02, "theta_s": 0.18, "alpha":  1.0, "n": 1.15,
        "Ks":   1.5, "z_nappe": -120.0, "theta_fc": 0.14, "theta_pwp": 0.04,
    },
]


def get_layer(z: float) -> dict:
    for layer in LAYERS:
        if layer["z_bottom"] <= z <= layer["z_top"]:
            return layer
    return LAYERS[-1]


# ======================================================================
# SIMULATION CONFIGURATION
# v4.0: Defaults reflect Norwegian shallow boreal forest soils
# ======================================================================

@dataclass
class SimConfig:
    sim_time : int   = 2500    # 28 years = Pole_MT default
    dt       : int   = 20

    res_x : int   = 50
    res_y : int   = 50
    res_z : int   = 40          # Fewer z-cells needed for -80 cm domain

    xmin  : float = -300.0;  xmax : float =  300.0
    ymin  : float = -300.0;  ymax : float =  300.0
    # v4.0 FIX: -80 cm not -200 cm (Norwegian soil depth from Shore 2025)
    zmin  : float =  -80.0;  zmax : float =    5.0

    n_pockets   : int   = 70
    random_seed : int   = 48

    absorption_rate  : float = 10e-6
    # v4.0: REW threshold from Shore (2025) / Granier et al. (1999)
    # REW < 0.40 = water limiting for Picea abies
    rew_stress_threshold : float = 0.40
    anoxia_threshold     : float = 0.02

    # Taproot (subType 1)
    w_grav_tap   : float = 1.2
    w_water_tap  : float = 0.4   # Shore (2025): weak soil-moisture coupling
    n_trials_tap : int   = 6
    sigma_tap    : float = 0.2

    # Lateral roots (subType 2)
    w_plagi_lat   : float = 3.0
    w_water_lat   : float = 0.5   # Shore (2025): soil moisture weak predictor
    w_radial_lat  : float = 0.5
    n_trials_lat  : int   = 8
    sigma_lat     : float = 0.3

    # Sinker roots (subType 4) — LOOP FIX v4.1:
    # Previous fix (sigma=0.50 + min_down_sin filter) still looped because:
    #   sigma=0.50 allows candidates only ~29° from current direction.
    #   From horizontal, the accepted candidate is ~29° below horizontal.
    #   This is mostly sideways, so the sinker curves gradually = loops.
    #
    # CORRECT FIX: sigma=0.80 (46° cone) + pure gravitropism (no filter).
    #   From horizontal, best candidate is ~46° below horizontal.
    #   From 46° below, best candidate is ~92° below = nearly vertical.
    #   → sinker reaches vertical in 2 STEPS (~0.6 cm arc) = no visible loop.
    #   n_trials=50 ensures the widest search always finds a steep candidate.
    w_grav_sink   : float = 10.0  # dominant gravity, no competition
    w_water_sink  : float = 0.1   # minimal — water must not fight gravity
    n_trials_sink : int   = 50    # many candidates for reliable steep downward
    sigma_sink    : float = 0.80  # wide cone → reaches steep angle in 2 steps

    # Standard fine roots (subType 3) — mineral soil
    w_water_fine  : float = 0.3   # Shore (2025): weak soil moisture coupling
    n_trials_fine : int   = 4
    sigma_fine    : float = 0.5

    # v4.0 NEW: Surface fine roots (subType 5) — O-layer
    # Shore (2025): O-layer holds 26-47% of total root length
    w_water_surf   : float = 0.3
    n_trials_surf  : int   = 4
    sigma_surf     : float = 0.6   # Wide spread for O-layer colonisation

    rock_min          : Optional[tuple] = None
    rock_max          : Optional[tuple] = None
    rock_rotation_deg : float           = -45.0

    # v4.0 FIX: -75 cm (Norwegian max soil depth)
    depth_limit : float = -75.0

    output_dir  : str  = "results/LayeredSoil"
    export_vtp  : bool = True
    export_vtk  : bool = False


# ======================================================================
# LAYERED SOIL — Norwegian profile with REW calculation
# ======================================================================

class LayeredSoil(pb.SoilLookUp):

    def __init__(self, cfg: SimConfig, obstacle=None,
                 z_nappe_offset: float = 0.0):
        super().__init__()
        self.cfg            = cfg
        self.obstacle       = obstacle
        self.z_nappe_offset = z_nappe_offset

        self._vg_params = {
            layer["label"]: vg.Parameters([
                layer["theta_r"], layer["theta_s"],
                layer["alpha"],   layer["n"], layer["Ks"],
            ])
            for layer in LAYERS
        }

        random.seed(cfg.random_seed)
        # Moisture pockets preferentially in O-layer and upper mineral soil
        self._pockets = []
        for _ in range(cfg.n_pockets):
            # 60% of pockets in top 20 cm (O-layer + upper mineral)
            if random.random() < 0.6:
                z_pk = random.uniform(cfg.zmax - 20, cfg.zmax - 2)
            else:
                z_pk = random.uniform(cfg.zmin + 20, cfg.zmax - 2)
            self._pockets.append((
                random.uniform(cfg.xmin, cfg.xmax),
                random.uniform(cfg.ymin, cfg.ymax),
                z_pk,
                random.uniform(5, 15),
                random.uniform(20, 60),
            ))

        print("Initialising Norwegian boreal soil grid …")
        t0 = time.time()
        self.grid = np.zeros((cfg.res_x, cfg.res_y, cfg.res_z))
        for i in range(cfg.res_x):
            x = cfg.xmin + i * (cfg.xmax - cfg.xmin) / (cfg.res_x - 1)
            for j in range(cfg.res_y):
                y = cfg.ymin + j * (cfg.ymax - cfg.ymin) / (cfg.res_y - 1)
                for k in range(cfg.res_z):
                    z = cfg.zmin + k * (cfg.zmax - cfg.zmin) / (cfg.res_z - 1)
                    self.grid[i, j, k] = self._compute_theta(x, y, z)

        self.grid_initial = np.copy(self.grid)
        # Track previous step moisture for MPdiff (Shore 2025)
        self.grid_previous = np.copy(self.grid)
        print(f"Norwegian soil grid ready in {time.time()-t0:.1f}s.")
        self._log_layer_stats()

    def _compute_theta(self, x, y, z):
        if self.obstacle is not None:
            if self.obstacle.getDist(pb.Vector3d(x, y, z)) < 0:
                return get_layer(z)["theta_r"]
        layer   = get_layer(z)
        params  = self._vg_params[layer["label"]]
        z_nappe = layer["z_nappe"] + self.z_nappe_offset
        h       = z_nappe - z - 5   # smaller offset for shallow soils
        bonus = max(
            (f * np.exp(-np.sqrt((x-px)**2 + (y-py)**2 + (z-pz)**2) / 30)
             for px, py, pz, _, f in self._pockets),
            default=0,
        )
        h = min(0.0, h + bonus - max(0, z + 5) * 3)
        return vg.water_content(h, params)

    def _idx(self, pos):
        cfg = self.cfg
        i = int((pos.x - cfg.xmin) / (cfg.xmax - cfg.xmin) * (cfg.res_x - 1))
        j = int((pos.y - cfg.ymin) / (cfg.ymax - cfg.ymin) * (cfg.res_y - 1))
        k = int((pos.z - cfg.zmin) / (cfg.zmax - cfg.zmin) * (cfg.res_z - 1))
        return (max(0, min(cfg.res_x-1, i)),
                max(0, min(cfg.res_y-1, j)),
                max(0, min(cfg.res_z-1, k)))

    def getWaterContent(self, pos):
        layer   = get_layer(pos.z)
        params  = self._vg_params[layer["label"]]
        theta_r = layer["theta_r"]
        z_nappe = layer["z_nappe"] + self.z_nappe_offset
        h = z_nappe - pos.z - 5
        bonus = max(
            (f * np.exp(-np.sqrt((pos.x-px)**2 + (pos.y-py)**2 + (pos.z-pz)**2) / 30)
             for px, py, pz, _, f in self._pockets),
            default=0,
        )
        h = min(0.0, h + bonus - max(0, pos.z + 5) * 3)
        theta_th = vg.water_content(h, params)
        i, j, k  = self._idx(pos)
        pumped   = self.grid_initial[i, j, k] - self.grid[i, j, k]
        return max(theta_r, theta_th - pumped)

    def getOxygen(self, pos):
        return get_layer(pos.z)["theta_s"] - self.getWaterContent(pos)

    def getREW(self, pos) -> float:
        """
        Relative Extractable Water at position.
        Source: Shore (2025), Granier et al. (1999).
        REW < 0.40 = water-limiting conditions for Picea abies.
        """
        layer = get_layer(pos.z)
        theta = self.getWaterContent(pos)
        return calculate_rew(theta, layer) if _STAGES_AVAILABLE else (
            (theta - layer["theta_r"]) / max(layer["theta_s"] - layer["theta_r"], 1e-6)
        )

    def getMPdiff(self, pos) -> float:
        """
        Rate of change of moisture (previous → current step).
        Shore (2025): MPdiff is strongest soil-based growth predictor
        (correlation -0.39 with dGROpc). Positive = soil rewetting.
        Used by hydrotropism: roots grow toward rewetting zones.
        """
        i, j, k = self._idx(pos)
        return self.grid[i, j, k] - self.grid_previous[i, j, k]

    def getValue(self, pos, organ=None):
        return self.getWaterContent(pos)

    def absorb(self, plant, dt):
        cfg    = self.cfg
        counts = {}
        for node in plant.getNodes():
            idx = self._idx(node)
            counts[idx] = counts.get(idx, 0) + 1

        # Save previous for MPdiff calculation
        self.grid_previous = np.copy(self.grid)

        for (i, j, k), n_nodes in counts.items():
            z       = cfg.zmin + k * (cfg.zmax - cfg.zmin) / (cfg.res_z - 1)
            layer   = get_layer(z)
            theta_r = layer["theta_r"]
            theta_s = layer["theta_s"]
            theta   = self.grid[i, j, k]
            air     = theta_s - theta
            if theta > theta_r and air >= cfg.anoxia_threshold:
                drop = min(cfg.absorption_rate * n_nodes * dt, 0.015)
                self.grid[i, j, k] = max(theta_r, theta - drop)

    def _log_layer_stats(self):
        cfg = self.cfg
        print(f"\n  Norwegian boreal soil moisture (Shore 2025 parameters):")
        print(f"  {'Layer':<10} {'z range':>16}  {'z_nappe':>8}  "
              f"{'mean θ':>7}  {'REW mean':>9}")
        print(f"  {'-'*60}")
        for layer in LAYERS:
            vals = []
            for k in range(cfg.res_z):
                z = cfg.zmin + k * (cfg.zmax - cfg.zmin) / (cfg.res_z - 1)
                if layer["z_bottom"] <= z <= layer["z_top"]:
                    vals.extend(self.grid[:, :, k].flatten())
            if vals:
                mean_theta = np.mean(vals)
                mean_rew   = calculate_rew(mean_theta, layer) if _STAGES_AVAILABLE else (
                    (mean_theta - layer["theta_r"]) /
                    max(layer["theta_s"] - layer["theta_r"], 1e-6))
                stress = " ← STRESS" if mean_rew < REW_STRESS_THRESHOLD else ""
                print(f"  {layer['label']:<10} "
                      f"({layer['z_top']:+5.0f}→{layer['z_bottom']:+5.0f} cm)  "
                      f"z_nappe={layer['z_nappe']:+5.0f}  "
                      f"θ={mean_theta:.3f}  "
                      f"REW={mean_rew:.2f}{stress}")
        print()


# ======================================================================
# TROPISM CLASSES
# ======================================================================

class GravHydroTropism(pb.Tropism):
    """
    Gravity + hydrotropism for taproot (subType 1).
    v4.0: REW-aware — when REW < 0.40, hydrotropism weight doubled.
    Shore (2025): strong hydrotropism when soil is water-limiting.
    """

    def __init__(self, plant, n_trials, sigma, soil,
                 w_grav=2.0, w_water=1.0, anoxia_thr=0.02,
                 rew_threshold=0.40):
        super().__init__(plant, n_trials, sigma)
        self.t_grav  = pb.Gravitropism(plant, n_trials, sigma)
        self.t_water = pb.Hydrotropism(plant, n_trials, sigma, soil)
        self.soil    = soil
        self.w_grav  = w_grav
        self.w_water = w_water
        self.anoxia  = anoxia_thr
        self.rew_thr = rew_threshold

    def tropismObjective(self, pos, old, a, b, dx, organ=None):
        pos_future = pb.Tropism.getPosition(pos, old, a, b, dx)
        sg = self.t_grav.tropismObjective(pos, old, a, b, dx, organ)

        # Anoxia escape (less relevant for Norwegian shallow soils, but retained)
        if self.soil.getOxygen(pos_future) < self.anoxia:
            return abs(sg - 0.5) * 2.0

        sw = self.t_water.tropismObjective(pos, old, a, b, dx, organ)

        # v4.1: REW-aware weighting (Shore 2025)
        # Multiplier reduced 2.0→1.3: paper shows WEAK soil moisture correlation
        # (Shore 2025 Fig 6.5: max correlation only 0.15-0.21 for soil indicators)
        # Internal tree water storage buffers soil moisture effects.
        rew = self.soil.getREW(pos)
        w_w = self.w_water * (1.3 if rew < self.rew_thr else 1.0)

        return (sg * self.w_grav + sw * w_w) / (self.w_grav + w_w)


class LateralPlagiotropism(pb.Tropism):
    """
    Active plagiotropism for lateral roots (subType 2).
    v4.0: REW-aware — hydrotropism strengthened under water stress.
    Shore (2025): Norwegian Picea laterals must stay in shallow O-layer
    + mineral soil (total 20-30 cm). Sag limited for shallow soils.
    """

    def __init__(self, plant, n_trials, sigma, soil,
                 w_plagi=3.0, w_water=1.0, w_radial=0.5,
                 anoxia_thr=0.02, rew_threshold=0.40, r_max=200.0):
        super().__init__(plant, n_trials, sigma)
        self.t_water  = pb.Hydrotropism(plant, n_trials, sigma, soil)
        self.t_grav   = pb.Gravitropism(plant, n_trials, sigma)
        self.soil     = soil
        self.w_plagi  = w_plagi
        self.w_water  = w_water
        self.w_radial = w_radial
        self.anoxia   = anoxia_thr
        self.rew_thr  = rew_threshold
        self.r_max    = r_max

    def tropismObjective(self, pos, old, a, b, dx, organ=None):
        pos_future = pb.Tropism.getPosition(pos, old, a, b, dx)

        if self.soil.getOxygen(pos_future) < self.anoxia:
            sg = self.t_grav.tropismObjective(pos, old, a, b, dx, organ)
            return abs(sg - 0.5) * 2.0

        # Plagiotropism: penalise downward movement
        # v4.0: reduced sag tolerance for Norwegian shallow soils
        # (laterals must stay in top 20-30 cm)
        r_now   = (pos.x**2 + pos.y**2)**0.5
        sag_tol = min(0.052 + 0.002 * r_now, 0.174)  # 3°→10° (tighter than v3)

        dz   = pos_future.z - pos.z
        dist = max(dx, 1e-9)
        sin_inc = dz / dist
        s_plagi = max(0.0, min(1.0, (sin_inc - sag_tol) / (1.0 - sag_tol)))

        sw = self.t_water.tropismObjective(pos, old, a, b, dx, organ)

        # Radial outward reward
        r_future = (pos_future.x**2 + pos_future.y**2)**0.5
        dr       = r_future - r_now
        s_radial = max(0.0, min(1.0, 1.0 - dr / self.r_max))

        # v4.1: REW-aware (Shore 2025 — weak soil moisture, internal buffer)
        rew = self.soil.getREW(pos)
        w_w = self.w_water * (1.3 if rew < self.rew_thr else 1.0)

        total_w = self.w_plagi + w_w + self.w_radial
        return (self.w_plagi  * s_plagi
              + w_w           * sw
              + self.w_radial * s_radial) / total_w


class SinkerTropism(pb.Tropism):
    """
    Sinker roots (subType 4): STRAIGHT DOWN from lateral plate.

    ═══════════════════════════════════════════════════════
    WHY ALL PREVIOUS VERSIONS STILL LOOPED:
    ═══════════════════════════════════════════════════════
    The insertion phi of a sinker is random (theta=pi/2 from horizontal
    lateral means the sinker starts in ANY perpendicular direction:
    up, down, sideways). This is correct plantbox behavior.

    v2: sigma=0.30, hard block dz>0  → sideways start → 13 steps to vertical → loop
    v3: sigma=0.12, hard block dz>0  → sideways start → 13 steps to vertical → loop
    v4: sigma=0.50, min_down=0.40    → accepted at 29° below horizontal →
                                       mostly sideways → gradual curve → loop

    THE ACTUAL MECHANISM OF LOOPS:
      With a hard filter (min_down_sin), the first accepted candidate is at
      the minimum angle (e.g. 24°–29° below horizontal). This IS "downward"
      but is mostly HORIZONTAL. The sinker then travels sideways for many
      steps creating the visible loop before curving to vertical.

    ═══════════════════════════════════════════════════════
    CORRECT FIX — Large sigma + pure gravitropism:
    ═══════════════════════════════════════════════════════
    sigma=0.80 means candidates span ±46° from the current direction.
    Pure gravitropism scores: 0 = straight down (best), 1 = straight up.

    From horizontal parent (sinker starts sideways at random phi):
      Step 1: candidates span ±46° from horizontal.
              Best candidate ≈ 46° below horizontal.
              Score improvement: from 0.5 (horizontal) to ≈ 0.24 (46° down).
      Step 2: candidates span ±46° from 46°-below.
              Best candidate ≈ 92° below horizontal = nearly vertical.
              Score ≈ 0.01 (almost straight down).
      → Sinker vertical in 2 STEPS = 0.6 cm of gentle arc. No visible loop.

    No hard filter needed — gravity competition alone is sufficient.
    n_trials=50 ensures the steepest candidate is always found.

    Anoxia escape: horizontal spread when waterlogged (clay / till).
    Shore (2025): w_water=0.1 because internal tree water storage buffers
                  soil moisture effects (weak soil moisture–growth correlation).
    """

    def __init__(self, plant, n_trials, sigma, soil,
                 w_grav=10.0, w_water=0.1, anoxia_thr=0.02):
        super().__init__(plant, n_trials, sigma)
        self.t_grav  = pb.Gravitropism(plant, n_trials, sigma)
        self.t_water = pb.Hydrotropism(plant, n_trials, sigma, soil)
        self.soil    = soil
        self.w_grav  = w_grav
        self.w_water = w_water
        self.anoxia  = anoxia_thr

    def tropismObjective(self, pos, old, a, b, dx, organ=None):
        pos_future = pb.Tropism.getPosition(pos, old, a, b, dx)

        # Anoxia / till boundary: spread horizontally
        if self.soil.getOxygen(pos_future) < self.anoxia:
            sg = self.t_grav.tropismObjective(pos, old, a, b, dx, organ)
            return abs(sg - 0.5) * 2.0   # 0 = horizontal, 1 = vertical

        # Pure gravitropism — NO hard filter.
        # With sigma=0.80 and n_trials=50, the steepest-downward candidate
        # among 50 is always selected. Reaches vertical in ~2 steps.
        sg = self.t_grav.tropismObjective(pos, old, a, b, dx, organ)
        sw = self.t_water.tropismObjective(pos, old, a, b, dx, organ)
        return (sg * self.w_grav + sw * self.w_water) / (self.w_grav + self.w_water)


class FineTropism(pb.Tropism):
    """
    Standard fine roots (subType 3): horizontal feeder roots in mineral soil.
    Penalises vertical movement. REW-aware hydrotropism.
    """

    def __init__(self, plant, n_trials, sigma, soil,
                 w_water=0.5, rew_threshold=0.40):
        super().__init__(plant, n_trials, sigma)
        self.t_water = pb.Hydrotropism(plant, n_trials, sigma, soil)
        self.soil    = soil
        self.w_water = w_water
        self.rew_thr = rew_threshold

    def tropismObjective(self, pos, old, a, b, dx, organ=None):
        pos_future = pb.Tropism.getPosition(pos, old, a, b, dx)
        dz     = abs(pos_future.z - pos.z)
        dist   = max(dx, 1e-9)
        s_vert = dz / dist
        sw     = self.t_water.tropismObjective(pos, old, a, b, dx, organ)
        rew    = self.soil.getREW(pos)
        w_w    = self.w_water * (1.3 if rew < self.rew_thr else 1.0)
        return (s_vert + w_w * sw) / (1.0 + w_w)


class SurfaceFineTropism(pb.Tropism):
    """
    Surface fine roots (subType 5): O-layer colonisation.
    NEW in v4.0 — based on Shore (2025) finding that 26-47% of root
    length is in the organic O-layer (0 to -8 cm).

    Behaviour:
      - Penalises movement BELOW -8 cm (out of O-layer)
      - Mild upward bias (to keep roots in O-layer)
      - Wide sigma for 3D spread through thin organic horizon
      - Weak hydrotropism (O-layer moisture is generally high)

    This creates a fuzzy cloud of fine roots in the organic layer,
    matching the field observation that most root length is there.
    """

    def __init__(self, plant, n_trials, sigma, soil,
                 w_water=0.3, o_layer_bottom=-8.0):
        super().__init__(plant, n_trials, sigma)
        self.t_water       = pb.Hydrotropism(plant, n_trials, sigma, soil)
        self.soil          = soil
        self.w_water       = w_water
        self.o_layer_bottom = o_layer_bottom

    def tropismObjective(self, pos, old, a, b, dx, organ=None):
        pos_future = pb.Tropism.getPosition(pos, old, a, b, dx)

        # Penalise movement deeper than O-layer bottom
        if pos_future.z < self.o_layer_bottom:
            # Strong penalty for leaving O-layer
            depth_penalty = min(1.0, (self.o_layer_bottom - pos_future.z) / 5.0)
            return 0.5 + 0.5 * depth_penalty

        # Within O-layer: mild horizontal + upward bias
        dz   = pos_future.z - pos.z
        dist = max(dx, 1e-9)
        # Penalise downward, reward horizontal/upward
        s_dir = max(0.0, min(1.0, 0.5 + dz / dist * 0.5))

        sw = self.t_water.tropismObjective(pos, old, a, b, dx, organ)
        return (s_dir + self.w_water * sw) / (1.0 + self.w_water)


# ======================================================================
# TIP TRACKER
# ======================================================================

class TipTracker:

    def __init__(self):
        self.tip_trajectories: dict = {}
        self.tip_positions:    dict = {}
        self.root_subtypes:    dict = {}

    def record(self, plant, t):
        snapshot = []
        for root in plant.getOrgans(pb.OrganTypes.root):
            rid   = root.getId()
            nodes = root.getNodes()
            if not nodes:
                continue
            tip = nodes[-1]
            if rid not in self.root_subtypes:
                try:
                    self.root_subtypes[rid] = root.param().subType
                except Exception:
                    self.root_subtypes[rid] = -1
            self.tip_trajectories.setdefault(rid, []).append(
                (tip.x, tip.y, tip.z, t))
            snapshot.append((rid, tip.x, tip.y, tip.z))
        self.tip_positions[t] = snapshot

    @property
    def final_tips(self):
        if not self.tip_positions:
            return []
        return self.tip_positions[max(self.tip_positions)]

    def summary(self):
        tips = self.final_tips
        n    = len(tips)
        print(f"\n{'='*66}")
        print(f"  Root Tip Summary — Norwegian Boreal Forest (Shore 2025 profile)")
        print(f"  Unique roots tracked : {len(self.tip_trajectories)}")
        print(f"  Timesteps recorded   : {len(self.tip_positions)}")
        if not tips:
            print(f"{'='*66}\n")
            return

        depths = [z for _, _, _, z in tips]
        radii  = [np.sqrt(x**2+y**2) for _, x, y, _ in tips]
        print(f"  Max depth            : {min(depths):+.1f} cm")
        print(f"  Mean depth           : {np.mean(depths):+.1f} cm")
        print(f"  Max lateral spread   : {max(radii):+.1f} cm")
        print(f"  Mean lateral spread  : {np.mean(radii):+.1f} cm")

        subtypes = {1:"Taproot  ", 2:"Lateral  ", 3:"Fine     ",
                    4:"Sinker   ", 5:"Surface  "}
        print(f"\n  Tips by root type:")
        print(f"  {'Type':<12} {'Count':>6} {'MaxDepth':>10} {'MaxRadius':>10}")
        print(f"  {'-'*42}")
        for st, label in subtypes.items():
            st_tips = [(x,y,z) for rid,x,y,z in tips
                       if self.root_subtypes.get(rid,-1) == st]
            if not st_tips:
                continue
            print(f"  {label:<12} {len(st_tips):>6} "
                  f"{min(z for _,_,z in st_tips):>+10.1f} cm  "
                  f"{max((x**2+y**2)**0.5 for x,y,_ in st_tips):>9.1f} cm")

        # Norwegian layer breakdown (Shore 2025)
        print(f"\n  Norwegian soil layers (Shore 2025):")
        norwegian = {
            "O-layer  ( 0→ -8cm)": sum(1 for _,_,_,z in tips if  -8.0 <= z <= 0),
            "Mineral  (-8→-35cm)": sum(1 for _,_,_,z in tips if -35.0 <= z < -8.0),
            "Till     (-35→-80cm)":sum(1 for _,_,_,z in tips if -80.0 <= z < -35.0),
        }
        for label, cnt in norwegian.items():
            pct = 100 * cnt / n if n else 0
            obs = "26–47%" if "O-layer" in label else ""
            print(f"    {label}: {cnt:4d} tips ({pct:.1f}%) {obs}")
        print(f"{'='*66}\n")


# ======================================================================
# AI DATASET EXPORTER
# ======================================================================

class AIDatasetExporter:

    def __init__(self):
        self.rows = []

    def build(self, tracker: TipTracker, soil: LayeredSoil):
        for rid, traj in tracker.tip_trajectories.items():
            if len(traj) < 2:
                continue
            sub_type = tracker.root_subtypes.get(rid, -1)
            for i in range(len(traj) - 1):
                x,  y,  z,  t  = traj[i]
                xn, yn, zn, _  = traj[i + 1]
                pos      = pb.Vector3d(x, y, z)
                pos_next = pb.Vector3d(xn, yn, zn)
                moisture = soil.getWaterContent(pos)
                oxygen   = soil.getOxygen(pos)
                rew      = soil.getREW(pos)
                mp_diff  = soil.getMPdiff(pos)    # Shore (2025) key predictor
                layer    = get_layer(z)
                self.rows.append({
                    "root_id"      : rid,
                    "sub_type"     : sub_type,
                    "t"            : t,
                    "x"            : round(x,  3),
                    "y"            : round(y,  3),
                    "z"            : round(z,  3),
                    "radius_cm"    : round((x**2+y**2)**0.5, 3),
                    "moisture"     : round(moisture, 4),
                    "oxygen"       : round(oxygen,   4),
                    # v4.0: REW and MPdiff from Shore (2025)
                    "REW"          : round(rew,      3),
                    "MP_diff"      : round(mp_diff,  5),
                    "rew_stress"   : int(rew < REW_STRESS_THRESHOLD),
                    "layer"        : layer["label"],
                    "dx"           : round(xn-x, 4),
                    "dy"           : round(yn-y, 4),
                    "dz"           : round(zn-z, 4),
                    "next_moisture": round(soil.getWaterContent(pos_next), 4),
                    "next_REW"     : round(soil.getREW(pos_next), 3),
                })

    def save(self, path="results/LayeredSoil/root_ai_dataset.csv"):
        if not self.rows:
            print("[WARNING] No data — run build() first.")
            return None
        os.makedirs(os.path.dirname(path), exist_ok=True)
        df = pd.DataFrame(self.rows)
        df.to_csv(path, index=False)
        print(f"\n  AI dataset → {path}")
        print(f"  Rows: {len(df)}  |  Cols: {list(df.columns)}")
        print(f"\n  Feature stats (Shore 2025 indicators highlighted):")
        show_cols = ["moisture", "REW", "MP_diff", "dx", "dy", "dz"]
        print(df[[c for c in show_cols if c in df]].describe().round(4).to_string())
        pct_stress = 100 * df["rew_stress"].mean()
        print(f"\n  REW stress (< {REW_STRESS_THRESHOLD}): {pct_stress:.1f}% of records")
        print(f"  (Shore 2025: REW < 0.40 = water-limiting for Picea abies)")
        return df


# ======================================================================
# EXPORT HELPERS
# ======================================================================

def export_soil_vtk(soil: LayeredSoil, step: int, out_dir: str):
    cfg   = soil.cfg
    fname = os.path.join(out_dir, f"soil_moisture_{step:04d}.vtk")
    dx    = (cfg.xmax - cfg.xmin) / (cfg.res_x - 1)
    dy    = (cfg.ymax - cfg.ymin) / (cfg.res_y - 1)
    dz    = (cfg.zmax - cfg.zmin) / (cfg.res_z - 1)
    with open(fname, "w") as f:
        f.write("# vtk DataFile Version 3.0\n")
        f.write(f"Norwegian boreal soil moisture — step {step}\n")
        f.write("ASCII\nDATASET STRUCTURED_POINTS\n")
        f.write(f"DIMENSIONS {cfg.res_x} {cfg.res_y} {cfg.res_z}\n")
        f.write(f"ORIGIN {cfg.xmin} {cfg.ymin} {cfg.zmin}\n")
        f.write(f"SPACING {dx:.4f} {dy:.4f} {dz:.4f}\n")
        f.write(f"POINT_DATA {cfg.res_x * cfg.res_y * cfg.res_z}\n")
        f.write("SCALARS theta float 1\nLOOKUP_TABLE default\n")
        for k in range(cfg.res_z):
            for j in range(cfg.res_y):
                for i in range(cfg.res_x):
                    f.write(f"{soil.grid[i,j,k]:.4f}\n")


def export_tips_csv(tracker: TipTracker, out_dir: str):
    path = os.path.join(out_dir, "tip_trajectories.csv")
    with open(path, "w") as f:
        f.write("root_id,sub_type,x,y,z,t,radius_cm,layer\n")
        for rid, entries in tracker.tip_trajectories.items():
            st = tracker.root_subtypes.get(rid, -1)
            for x, y, z, t in entries:
                r = (x**2+y**2)**0.5
                f.write(f"{rid},{st},{x:.3f},{y:.3f},{z:.3f},"
                        f"{t:.1f},{r:.3f},{get_layer(z)['label']}\n")
    print(f"  Tip trajectories → {path}")


def export_subtype_csv(tracker: TipTracker, out_dir: str):
    """Per-subtype CSVs for easy ParaView colouring by root type."""
    names = {1:"taproot", 2:"lateral", 3:"fine", 4:"sinker", 5:"surface"}
    for st, name in names.items():
        rows = [
            {"root_id": rid, "x": round(x,3), "y": round(y,3),
             "z": round(z,3), "t": t,
             "radius_cm": round((x**2+y**2)**0.5, 3),
             "layer": get_layer(z)["label"]}
            for rid, entries in tracker.tip_trajectories.items()
            if tracker.root_subtypes.get(rid,-1) == st
            for x, y, z, t in entries
        ]
        if rows:
            path = os.path.join(out_dir, f"tips_{name}.csv")
            pd.DataFrame(rows).to_csv(path, index=False)
            print(f"  Subtype CSV → {path} ({len(rows)} rows)")


def _progress(i, n, t_start, n_tips):
    pct  = 100 * (i+1) / n
    bar  = "█" * int(pct/5) + "░" * (20 - int(pct/5))
    elapsed = time.time() - t_start
    eta  = elapsed / (i+1) * (n-i-1) if i > 0 else 0
    sys.stdout.write(
        f"\r  [{bar}] {pct:5.1f}%  step {i+1}/{n}  "
        f"tips={n_tips:4d}  ETA {eta:5.0f}s  ")
    sys.stdout.flush()


def check_allometry(results: dict, stage_name: Optional[str] = None,
                    soil_depth_cm: float = 25.0):
    if not (_STAGES_AVAILABLE and stage_name):
        return
    from tree_stages import d0_to_root_extent, TREE_STAGES
    stage   = next((s for s in TREE_STAGES if s["stage"] == stage_name), None)
    if not stage:
        return
    extents = d0_to_root_extent(stage["d0_mean"], soil_depth_cm)
    tips    = results["final_tips"]
    if not tips:
        return
    depths = [z for _,_,_,z in tips]
    radii  = [(x**2+y**2)**0.5 for _,x,y,_ in tips]
    lat_r  = max(radii) / extents["lateral_spread_cm"]
    dep_r  = abs(min(depths)) / max(extents["max_depth_cm"], 1.0)
    def rating(r): return "GOOD ✓" if 0.6<=r<=1.4 else ("FAIR ~" if 0.4<=r<=1.6 else "POOR ✗")
    print(f"\n  Allometry Check — {stage_name} (D0={stage['d0_mean']} cm, "
          f"soil={soil_depth_cm} cm)")
    if extents["soil_constrained"]:
        print(f"  [NOTE] Depth is SOIL-CONSTRAINED (Norwegian median {soil_depth_cm} cm)")
    print(f"  {'Lateral spread':<20} {max(radii):>8.1f}cm  "
          f"expected {extents['lateral_spread_cm']:>6.0f}cm  "
          f"ratio={lat_r:.2f}  {rating(lat_r)}")
    print(f"  {'Max depth':<20} {abs(min(depths)):>8.1f}cm  "
          f"expected {extents['max_depth_cm']:>6.0f}cm  "
          f"ratio={dep_r:.2f}  {rating(dep_r)}")


# ======================================================================
# MAIN PREDICTOR
# ======================================================================

class RootGrowthPredictor:

    def __init__(self, xml_path: str, cfg: Optional[SimConfig] = None,
                 stage_name: Optional[str] = None,
                 z_nappe_offset: float = 0.0,
                 soil_depth_cm: float = 25.0):
        self.xml_path      = xml_path
        self.cfg           = cfg or SimConfig()
        self.stage_name    = stage_name
        self.z_nappe_offset = z_nappe_offset
        self.soil_depth_cm = soil_depth_cm
        os.makedirs(self.cfg.output_dir, exist_ok=True)

    def _validate_xml(self):
        if not os.path.isfile(self.xml_path):
            raise FileNotFoundError(
                f"\n[ERROR] XML not found: {self.xml_path}\n"
                f"  CWD: {os.getcwd()}")
        print(f"  XML: {os.path.abspath(self.xml_path)}")

    def _build_obstacles(self):
        cfg = self.cfg; parts = []; has_rock = False
        if cfg.rock_min and cfg.rock_max:
            rock = pb.SDF_Cuboid(pb.Vector3d(*cfg.rock_min),
                                 pb.Vector3d(*cfg.rock_max))
            rock = pb.SDF_RotateTranslate(rock, cfg.rock_rotation_deg, 2,
                                          pb.Vector3d(0,0,0))
            parts.append(rock); has_rock = True
        if cfg.depth_limit is not None:
            z = cfg.depth_limit
            floor = pb.SDF_HalfPlane(
                pb.Vector3d(-1000,-1000,z), pb.Vector3d(1000,-1000,z),
                pb.Vector3d(-1000,1000,z))
            parts.append(floor)
        if not parts:
            return None, False
        obstacle = parts[0]
        for p in parts[1:]:
            obstacle = pb.SDF_Union(obstacle, p)
        return obstacle, has_rock

    def run(self) -> dict:
        self._validate_xml()
        cfg = self.cfg

        print("\n" + "="*68)
        print("  PICEA ABIES ROOT GROWTH — v4.0 (Norwegian Boreal, Shore 2025)")
        print("="*68)
        print(f"  Stage      : {self.stage_name or 'custom'}")
        print(f"  sim_time   : {cfg.sim_time} days ({cfg.sim_time/365:.1f} yr)")
        print(f"  Domain     : z=[{cfg.zmin:.0f},{cfg.zmax:.0f}] cm  "
              f"(Norwegian soil: max 80 cm)")
        print(f"  depth_limit: {cfg.depth_limit:.0f} cm")
        print(f"  REW stress : < {cfg.rew_stress_threshold} (Shore 2025)")
        print(f"\n  Soil profile (Shore 2025, Table 3.1 averages):")
        for L in LAYERS:
            print(f"    {L['label']:<10} ({L['z_top']:+5.0f}→{L['z_bottom']:+5.0f} cm)"
                  f"  θs={L['theta_s']:.2f}  α={L['alpha']:5.1f}")
        print(f"\n  Tropism assignments:")
        print(f"    subType 1 (Taproot)       → GravHydroTropism    "
              f"w_grav={cfg.w_grav_tap}")
        print(f"    subType 2 (Lateral)       → LateralPlagiotropism "
              f"w_plagi={cfg.w_plagi_lat}")
        print(f"    subType 3 (Fine)          → FineTropism          "
              f"w_water={cfg.w_water_fine}")
        print(f"    subType 4 (Sinker)        → SinkerTropism        "
              f"sigma={cfg.sigma_sink}")
        print(f"    subType 5 (Surface/O-lyr) → SurfaceFineTropism   "
              f"[NEW Shore 2025]")
        print("="*68)

        obstacle, has_rock = self._build_obstacles()
        if has_rock:
            vp.write_container(obstacle,
                               os.path.join(cfg.output_dir, "obstacles.vtp"))

        domain    = pb.SDF_PlantBox(2000, 2000, 2000)
        navigable = pb.SDF_Difference(domain, obstacle) if obstacle else domain

        soil  = LayeredSoil(cfg, obstacle=obstacle,
                            z_nappe_offset=self.z_nappe_offset)
        plant = pb.Plant()
        plant.readParameters(self.xml_path)

        # Apply sim_time from Python config
        seed_params = plant.getOrganRandomParameter(pb.OrganTypes.seed)
        if seed_params:
            seed_params[0].simtime = cfg.sim_time
            print(f"\n  [OK] seed_params[0].simtime = {cfg.sim_time} days applied.")

        plant.setSoil(soil)
        plant.setGeometry(navigable)
        plant.initialize()

        # Build tropism objects
        t_tap = GravHydroTropism(
            plant, cfg.n_trials_tap, cfg.sigma_tap, soil,
            cfg.w_grav_tap, cfg.w_water_tap,
            cfg.anoxia_threshold, cfg.rew_stress_threshold)

        t_lat = LateralPlagiotropism(
            plant, cfg.n_trials_lat, cfg.sigma_lat, soil,
            cfg.w_plagi_lat, cfg.w_water_lat, cfg.w_radial_lat,
            cfg.anoxia_threshold, cfg.rew_stress_threshold)

        t_fine = FineTropism(
            plant, cfg.n_trials_fine, cfg.sigma_fine, soil,
            cfg.w_water_fine, cfg.rew_stress_threshold)

        t_sink = SinkerTropism(
            plant, cfg.n_trials_sink, cfg.sigma_sink, soil,
            cfg.w_grav_sink, cfg.w_water_sink,
            cfg.anoxia_threshold)

        # v4.0 NEW: Surface fine roots for O-layer
        t_surf = SurfaceFineTropism(
            plant, cfg.n_trials_surf, cfg.sigma_surf, soil,
            cfg.w_water_surf, o_layer_bottom=-8.0)

        for t in [t_tap, t_lat, t_fine, t_sink, t_surf]:
            t.setGeometry(navigable)

        plant.setTropism(t_tap,  pb.OrganTypes.root, 1)
        plant.setTropism(t_lat,  pb.OrganTypes.root, 2)
        plant.setTropism(t_fine, pb.OrganTypes.root, 3)
        plant.setTropism(t_sink, pb.OrganTypes.root, 4)

        # subType 5 (O-layer surface fine roots) only exists if the new XML
        # is loaded. Check before assigning to stay backwards-compatible
        # with older XML files that only define subTypes 1-4.
        _loaded_subtypes = {
            p.subType
            for p in plant.getOrganRandomParameter(pb.OrganTypes.root)
        }
        if 5 in _loaded_subtypes:
            plant.setTropism(t_surf, pb.OrganTypes.root, 5)
            print("  [OK] SurfaceFineTropism assigned to subType 5 (O-layer roots).")
        else:
            print("  [INFO] subType 5 not found in XML — O-layer surface roots disabled.")
            print("         Using updated Picea_Abies_simple.xml (v4.0) to enable them.")

        tracker = TipTracker()
        n_steps = round(cfg.sim_time / cfg.dt)
        t_start = time.time()
        print(f"\n  Running {n_steps} steps × {cfg.dt} d = {cfg.sim_time} d …\n")

        for i in range(n_steps):
            t = (i + 1) * cfg.dt
            plant.simulate(cfg.dt)
            soil.absorb(plant, cfg.dt)
            tracker.record(plant, t)

            if cfg.export_vtp:
                plant.write(os.path.join(cfg.output_dir, f"roots_{i:04d}.vtp"))
            if cfg.export_vtk:
                export_soil_vtk(soil, i, cfg.output_dir)

            n_tips = len(tracker.tip_positions.get(t, []))
            _progress(i, n_steps, t_start, n_tips)

        print(f"\n\n  Simulation complete in {time.time()-t_start:.1f}s.")

        tracker.summary()
        export_tips_csv(tracker, cfg.output_dir)
        export_subtype_csv(tracker, cfg.output_dir)
        check_allometry({"final_tips": tracker.final_tips},
                        self.stage_name, self.soil_depth_cm)

        return {
            "tip_trajectories": tracker.tip_trajectories,
            "tip_positions"   : tracker.tip_positions,
            "final_tips"      : tracker.final_tips,
            "tracker"         : tracker,
            "soil"            : soil,
        }


# ======================================================================
# CONVENIENCE FUNCTIONS
# ======================================================================

def tips_as_array(results, timestep): ...
def trajectory_as_array(results, root_id): ...

def deepest_tips(results, n=5):
    return sorted(results["final_tips"], key=lambda r: r[3])[:n]

def widest_tips(results, n=5):
    return sorted(results["final_tips"],
                  key=lambda r: -(r[1]**2+r[2]**2)**0.5)[:n]

def tips_by_subtype(results, tracker, sub_type):
    return [(rid,x,y,z) for rid,x,y,z in results["final_tips"]
            if tracker.root_subtypes.get(rid,-1) == sub_type]


# ======================================================================
# ENTRY POINT
# ======================================================================

if __name__ == "__main__":

    # Shore (2025) monitored plots: 35-49 year old trees → Pole_OMT/Pole_VT
    # Using Pole_MT (28 yr) here as default; change to "Pole_OMT" or "Pole_VT"
    # for Norwegian field-data-equivalent simulation.
    STAGE_NAME    = "Pole_MT"
    SOIL_DEPTH_CM = 25.0   # Norwegian median from Shore (2025), Fig 3.6

    if _STAGES_AVAILABLE and STAGE_NAME is not None:
        print(f"\nLoading SimConfig for stage: {STAGE_NAME}")
        print(f"Norwegian soil depth constraint: {SOIL_DEPTH_CM} cm")
        overrides = get_sim_config_for_stage(STAGE_NAME, SOIL_DEPTH_CM)
        print(f"  sim_time  = {overrides['sim_time']} d ({overrides['sim_time']//365} yr)")
        print(f"  zmin      = {overrides['zmin']} cm (Norwegian constraint)")

        stage_info     = next(s for s in TREE_STAGES if s["stage"] == STAGE_NAME)
        site           = SITE_MODIFIERS[stage_info["site_type"]]
        z_nappe_offset = site["z_nappe_offset"]

        cfg = SimConfig(**{k: v for k, v in overrides.items()
                           if k in SimConfig.__dataclass_fields__})
    else:
        if STAGE_NAME and not _STAGES_AVAILABLE:
            print("[WARNING] tree_stages.py not found — using SimConfig defaults.")
        cfg = SimConfig()
        z_nappe_offset = 0.0

    print(f"\n  [CHECK] cfg.sim_time = {cfg.sim_time} days "
          f"({cfg.sim_time/365:.1f} yr)")
    print(f"  [CHECK] cfg.zmin = {cfg.zmin} cm "
          f"({'OK — Norwegian soil depth' if cfg.zmin >= -100 else 'WARNING: too deep for Norwegian soils'})")

    predictor = RootGrowthPredictor(
        XML_PATH, cfg=cfg,
        stage_name=STAGE_NAME,
        z_nappe_offset=z_nappe_offset,
        soil_depth_cm=SOIL_DEPTH_CM,
    )
    results = predictor.run()

    final = results["final_tips"]
    print(f"\nTotal tips: {len(final)}")

    print("\nTop 5 deepest tips:")
    for rid, x, y, z in deepest_tips(results, 5):
        print(f"  Root {rid:4d}  ({x:+7.1f},{y:+7.1f},{z:+7.1f}) cm  "
              f"[{get_layer(z)['label']}]")

    print("\nTop 5 widest laterals:")
    for rid, x, y, z in widest_tips(results, 5):
        r  = (x**2+y**2)**0.5
        st = results["tracker"].root_subtypes.get(rid,-1)
        print(f"  Root {rid:4d}  ({x:+7.1f},{y:+7.1f},{z:+7.1f}) cm  "
              f"r={r:.1f}  subType={st}  [{get_layer(z)['label']}]")

    if _STAGES_AVAILABLE and STAGE_NAME:
        stage_info = next(s for s in TREE_STAGES if s["stage"] == STAGE_NAME)
        validate_depth_distribution(final, STAGE_NAME, stage_info["site_type"])

    print("\nBuilding AI dataset (REW + MPdiff from Shore 2025) …")
    exporter = AIDatasetExporter()
    exporter.build(results["tracker"], results["soil"])
    exporter.save(os.path.join(cfg.output_dir, "root_ai_dataset.csv"))