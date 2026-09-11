"""
root_growth_predictor.py
========================
Picea Abies root growth simulation with a realistic LAYERED soil profile.

Soil layers:
    0   →  -40 cm  :  Loam   (balanced drainage, good branching)
   -40  → -100 cm  :  Sand   (fast drainage, roots chase water deeper)
  -100  → -200 cm  :  Clay   (waterlogged, anoxia escape → horizontal growth)

Each layer has its own Van Genuchten parameters AND its own water table depth.

Features:
  - Per-layer Van Genuchten hydraulics + per-layer water table
  - Dynamic water depletion by roots (absorb())
  - Mixed tropism: gravitropism + hydrotropism + anoxia escape (plagiotropism)
  - Separate tropism for structural roots (subType 2, 3) vs sinker roots (subType 5)
  - Rock obstacle + hard depth floor
  - Exports: VTP (roots), VTK (soil moisture), CSV (tip trajectories with layer label)
  - Output: results/LayeredSoil/

Usage:
    python root_growth_predictor.py
"""

import os
import random
import pathlib
from dataclasses import dataclass
from typing import Optional

import numpy as np
import plantbox as pb
import plantbox.visualisation.vtk_plot as vp
import plantbox.functional.van_genuchten as vg


# ======================================================================
# XML PATH
# Resolved relative to this script file — works regardless of where
# you call the script from.
# If it still fails, replace XML_PATH with a full absolute path:
#   XML_PATH = "/workspaces/.../Picea_Abies_hydro_v2.xml"
# ======================================================================
_HERE    = pathlib.Path(__file__).resolve().parent
XML_PATH = str(
    (_HERE / "../../modelparameter_TreeCity/structural/Picea_Abies_hydro_v2.xml")
    .resolve()
)


# ======================================================================
# LAYERED SOIL DEFINITION
# ======================================================================
# Each layer dict contains:
#   label    : name used in logs and CSV output
#   z_top    : upper boundary (cm)  — note: z=0 is soil surface
#   z_bottom : lower boundary (cm)
#   theta_r  : residual water content      (cm³/cm³)
#   theta_s  : saturated water content     (cm³/cm³) = porosity
#   alpha    : Van Genuchten α             (1/cm)    — larger = drains faster
#   n        : Van Genuchten n             (-)       — pore-size distribution
#   Ks       : saturated hydraulic cond.   (cm/day)
#   z_nappe  : water table depth for THIS layer (cm) — KEY enhancement
#
# Biological effects:
#   Loam  → moderate drainage, balanced O₂/water → good root branching
#   Sand  → fast drainage, low retention, deep water table
#           → roots grow deep chasing moisture; high O₂, no anoxia
#   Clay  → very slow drainage, shallow water table, waterlogged
#           → anoxia threshold crossed → roots escape horizontally
#
# Reference: Carsel & Parrish (1988) — standard Van Genuchten parameters

LAYERS = [
    {
        "label"   : "Loam",
        "z_top"   :    0.0,
        "z_bottom":  -40.0,
        "theta_r" : 0.078,
        "theta_s" : 0.430,
        "alpha"   : 0.036,
        "n"       : 1.560,
        "Ks"      : 24.96,
        "z_nappe" : -150.0,   # water table well below → moderate suction
    },
    {
        "label"   : "Sand",
        "z_top"   :  -40.0,
        "z_bottom": -100.0,
        "theta_r" : 0.045,
        "theta_s" : 0.430,
        "alpha"   : 0.145,
        "n"       : 2.680,
        "Ks"      : 712.8,
        "z_nappe" : -180.0,   # deep water table → strong downward pull on roots
    },
    {
        "label"   : "Clay",
        "z_top"   : -100.0,
        "z_bottom": -200.0,
        "theta_r" : 0.068,
        "theta_s" : 0.380,
        "alpha"   : 0.008,
        "n"       : 1.090,
        "Ks"      :   4.8,
        "z_nappe" :  -80.0,   # shallow water table → near-saturated, anoxia risk
    },
]


def get_layer(z: float) -> dict:
    """
    Return the layer dict for depth z (cm).
    Falls back to the deepest layer (Clay) if z is below all defined layers.
    """
    for layer in LAYERS:
        if layer["z_bottom"] <= z <= layer["z_top"]:
            return layer
    return LAYERS[-1]   # below all layers → Clay properties


# ======================================================================
# SIMULATION CONFIGURATION
# ======================================================================
@dataclass
class SimConfig:
    # Time
    sim_time : int   = 5000   # total simulation days
    dt       : int   = 50     # timestep (days)

    # Soil grid — covers full root lateral extent (horizontal roots: 600 cm)
    res_x : int   = 50
    res_y : int   = 50
    res_z : int   = 40
    xmin  : float = -700.0;  xmax : float =  700.0
    ymin  : float = -700.0;  ymax : float =  700.0
    zmin  : float = -200.0;  zmax : float =    5.0

    # Random moisture pockets (wet patches near surface)
    n_pockets   : int   = 70
    random_seed : int   = 48

    # Root water absorption
    absorption_rate  : float = 10e-6
    anoxia_threshold : float = 0.05   # air fraction below which anoxia fires

    # Tropism — structural roots (subType 2 = Fine, 3 = Sinkers)
    w_grav_main   : float = 0.7
    w_water_main  : float = 1.0
    n_trials_main : float = 2.0
    sigma_main    : float = 1.5

    # Tropism — short sinker roots (subType 5)
    w_grav_fine   : float = 0.0025
    w_water_fine  : float = 1.0
    n_trials_fine : float = 4.0
    sigma_fine    : float = 0.5

    # Rock obstacle (tilted cuboid)
    rock_min          : tuple = (50, -150, -150)
    rock_max          : tuple = (100, 100, -20)
    rock_rotation_deg : float = -45.0

    # Hard depth floor (set to None to disable)
    depth_limit : float = -150.0

    # Output
    output_dir : str  = "results/LayeredSoil"
    export_vtp : bool = True    # per-step root VTP files
    export_vtk : bool = True    # per-step soil moisture VTK files


# ======================================================================
# DYNAMIC LAYERED SOIL
# ======================================================================
class LayeredSoil(pb.SoilLookUp):
    """
    3-D soil grid where hydraulic properties AND water table depth
    change with depth according to the LAYERS definition.

    Key enhancement over single-soil version:
        Each layer has its own z_nappe (water table depth).
        This means:
          - Loam  (0→-40):   z_nappe=-150 → moderate suction, well-drained
          - Sand  (-40→-100): z_nappe=-180 → strong capillary pull downward
          - Clay  (-100→-200): z_nappe=-80 → near-saturated, low suction

        Using the correct z_nappe per layer gives a physically realistic
        moisture gradient with a sharp dry zone in sand and a wet
        waterlogged zone in clay.
    """

    def __init__(self, cfg: SimConfig, obstacle=None):
        super().__init__()
        self.cfg      = cfg
        self.obstacle = obstacle

        # Pre-build one vg.Parameters object per layer label for fast lookup
        self._vg_params = {
            layer["label"]: vg.Parameters([
                layer["theta_r"],
                layer["theta_s"],
                layer["alpha"],
                layer["n"],
                layer["Ks"],
            ])
            for layer in LAYERS
        }

        # Random moisture pockets
        random.seed(cfg.random_seed)
        self._pockets = [
            (
                random.uniform(cfg.xmin, cfg.xmax),
                random.uniform(cfg.ymin, cfg.ymax),
                random.uniform(cfg.zmax - 100, cfg.zmax - 10),
                random.uniform(25, 50),    # radius (not used in formula but kept)
                random.uniform(40, 120),   # force  (pressure-head bonus)
            )
            for _ in range(cfg.n_pockets)
        ]

        # Fill 3-D moisture grid
        print("Initialising layered soil grid …")
        self.grid = np.zeros((cfg.res_x, cfg.res_y, cfg.res_z))
        for i in range(cfg.res_x):
            x = cfg.xmin + i * (cfg.xmax - cfg.xmin) / (cfg.res_x - 1)
            for j in range(cfg.res_y):
                y = cfg.ymin + j * (cfg.ymax - cfg.ymin) / (cfg.res_y - 1)
                for k in range(cfg.res_z):
                    z = cfg.zmin + k * (cfg.zmax - cfg.zmin) / (cfg.res_z - 1)
                    self.grid[i, j, k] = self._compute_theta(x, y, z)

        self.grid_initial = np.copy(self.grid)
        print("Layered soil grid ready.")
        self._log_layer_stats()

    # ------------------------------------------------------------------
    def _compute_theta(self, x: float, y: float, z: float) -> float:
        """
        Compute initial Van Genuchten water content at (x, y, z).
        Uses the z_nappe of the layer that (x,y,z) belongs to.
        """
        cfg = self.cfg

        # Inside an obstacle → dry (residual moisture only)
        if self.obstacle is not None:
            if self.obstacle.getDist(pb.Vector3d(x, y, z)) < 0:
                return get_layer(z)["theta_r"]

        layer   = get_layer(z)
        params  = self._vg_params[layer["label"]]
        z_nappe = layer["z_nappe"]          # ← per-layer water table

        # Base pressure head: distance below the local water table
        h = z_nappe - z - 20

        # Bonus from moisture pockets (Gaussian decay with distance)
        bonus = max(
            (f * np.exp(-np.sqrt((x-px)**2 + (y-py)**2 + (z-pz)**2) / 50)
             for px, py, pz, _, f in self._pockets),
            default=0,
        )

        # Surface evaporation: top 30 cm dries out progressively
        h = min(0.0, h + bonus - max(0, z + 30) * 5)

        return vg.water_content(h, params)

    # ------------------------------------------------------------------
    def _idx(self, pos) -> tuple:
        """Convert continuous 3-D position to clamped grid indices."""
        cfg = self.cfg
        i = int((pos.x - cfg.xmin) / (cfg.xmax - cfg.xmin) * (cfg.res_x - 1))
        j = int((pos.y - cfg.ymin) / (cfg.ymax - cfg.ymin) * (cfg.res_y - 1))
        k = int((pos.z - cfg.zmin) / (cfg.zmax - cfg.zmin) * (cfg.res_z - 1))
        return (
            max(0, min(cfg.res_x - 1, i)),
            max(0, min(cfg.res_y - 1, j)),
            max(0, min(cfg.res_z - 1, k)),
        )

    # ------------------------------------------------------------------
    def getWaterContent(self, pos) -> float:
        """
        Real-time water content at pos.
        Uses the per-layer z_nappe and subtracts water already pumped
        by roots from that voxel.
        """
        layer   = get_layer(pos.z)
        params  = self._vg_params[layer["label"]]
        theta_r = layer["theta_r"]
        z_nappe = layer["z_nappe"]          # ← per-layer water table

        h = z_nappe - pos.z - 20

        bonus = max(
            (f * np.exp(-np.sqrt((pos.x-px)**2 + (pos.y-py)**2 + (pos.z-pz)**2) / 50)
             for px, py, pz, _, f in self._pockets),
            default=0,
        )

        h = min(0.0, h + bonus - max(0, pos.z + 30) * 5)
        theta_th = vg.water_content(h, params)

        # Subtract water already pumped from this voxel
        i, j, k = self._idx(pos)
        pumped   = self.grid_initial[i, j, k] - self.grid[i, j, k]

        return max(theta_r, theta_th - pumped)

    # ------------------------------------------------------------------
    def getOxygen(self, pos) -> float:
        """Air-filled porosity = theta_s(layer) − theta_current."""
        theta_s = get_layer(pos.z)["theta_s"]
        return theta_s - self.getWaterContent(pos)

    # ------------------------------------------------------------------
    def getValue(self, pos, organ=None) -> float:
        """Required by pb.SoilLookUp — returns moisture [0, 1]."""
        return self.getWaterContent(pos)

    # ------------------------------------------------------------------
    def absorb(self, plant, dt: float):
        """
        Drain water from every voxel occupied by root nodes.
        Uses the per-layer theta_r and theta_s.
        Anoxia rule: skip waterlogged voxels (not enough air).
        """
        cfg    = self.cfg
        counts = {}
        for node in plant.getNodes():
            idx = self._idx(node)
            counts[idx] = counts.get(idx, 0) + 1

        for (i, j, k), n_nodes in counts.items():
            # Determine which layer this voxel belongs to
            z       = cfg.zmin + k * (cfg.zmax - cfg.zmin) / (cfg.res_z - 1)
            layer   = get_layer(z)
            theta_r = layer["theta_r"]
            theta_s = layer["theta_s"]
            theta   = self.grid[i, j, k]
            air     = theta_s - theta

            if theta > theta_r and air >= cfg.anoxia_threshold:
                drop = min(cfg.absorption_rate * n_nodes * dt, 0.02)
                self.grid[i, j, k] = max(theta_r, theta - drop)

    # ------------------------------------------------------------------
    def _log_layer_stats(self):
        """Print mean/min/max moisture per layer after initialisation."""
        cfg = self.cfg
        print("\n  Initial moisture statistics by layer:")
        print(f"  {'Layer':<6} {'z range':>18}  "
              f"{'z_nappe':>8}  {'mean θ':>7}  {'min θ':>7}  {'max θ':>7}")
        print(f"  {'-'*60}")
        for layer in LAYERS:
            vals = []
            for k in range(cfg.res_z):
                z = cfg.zmin + k * (cfg.zmax - cfg.zmin) / (cfg.res_z - 1)
                if layer["z_bottom"] <= z <= layer["z_top"]:
                    vals.extend(self.grid[:, :, k].flatten())
            if vals:
                print(f"  {layer['label']:<6} "
                      f"({layer['z_top']:+5.0f} → {layer['z_bottom']:+5.0f} cm)  "
                      f"z_nappe={layer['z_nappe']:+6.0f}  "
                      f"θ={np.mean(vals):.3f}  "
                      f"min={np.min(vals):.3f}  "
                      f"max={np.max(vals):.3f}")
        print()


# ======================================================================
# MIXED TROPISM
# ======================================================================
class TropismeMixte(pb.Tropism):
    """
    Scores N candidate growth directions and selects the best.

    Normal:
        score = (w_grav × gravity + w_water × hydro) / (w_grav + w_water)
        Lower score = better direction (0 = perfect).

    Anoxia escape (fires in waterlogged clay):
        score = |gravity_score − 0.5| × 2
        → horizontal direction (gravity_score ≈ 0.5) scores 0 = perfect
        → both up and down score badly → root escapes sideways
    """

    def __init__(self, plant, n_trials, sigma, soil,
                 w_grav=0.7, w_water=1.0, anoxia_thr=0.05):
        super().__init__(plant, n_trials, sigma)
        self.t_grav  = pb.Gravitropism(plant, n_trials, sigma)
        self.t_water = pb.Hydrotropism(plant, n_trials, sigma, soil)
        self.soil    = soil
        self.w_grav  = w_grav
        self.w_water = w_water
        self.anoxia  = anoxia_thr

    def tropismObjective(self, pos, old, a, b, dx, organ=None):
        sg         = self.t_grav.tropismObjective(pos, old, a, b, dx, organ)
        sw         = self.t_water.tropismObjective(pos, old, a, b, dx, organ)
        pos_future = pb.Tropism.getPosition(pos, old, a, b, dx)

        if self.soil.getOxygen(pos_future) < self.anoxia:
            return abs(sg - 0.5) * 2.0     # plagiotropism

        return (sg * self.w_grav + sw * self.w_water) / (self.w_grav + self.w_water)


# ======================================================================
# TIP TRACKER
# ======================================================================
class TipTracker:
    """
    Records every root apex position at each timestep.

    tip_trajectories : {root_id: [(x, y, z, t), ...]}
    tip_positions    : {t:       [(root_id, x, y, z), ...]}
    """

    def __init__(self):
        self.tip_trajectories: dict = {}
        self.tip_positions:    dict = {}

    def record(self, plant, t: float):
        snapshot = []
        for root in plant.getOrgans(pb.OrganTypes.root):
            rid   = root.getId()
            nodes = root.getNodes()
            if not nodes:
                continue
            tip = nodes[-1]
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
        print(f"\n{'='*56}")
        print(f"  Layered Soil — Root Tip Summary")
        print(f"  Unique roots tracked : {len(self.tip_trajectories)}")
        print(f"  Timesteps recorded   : {len(self.tip_positions)}")
        if tips:
            depths = [z for _, _, _, z in tips]
            print(f"  Max depth reached    : {min(depths):+.1f} cm")
            print(f"  Mean final depth     : {np.mean(depths):+.1f} cm")
            print(f"\n  Tips per layer at end of simulation:")
            for layer in LAYERS:
                n = sum(1 for _, _, _, z in tips
                        if layer["z_bottom"] <= z <= layer["z_top"])
                pct = 100 * n / len(tips) if tips else 0
                print(f"    {layer['label']:<6} "
                      f"({layer['z_top']:+5.0f}→{layer['z_bottom']:+5.0f} cm): "
                      f"{n:4d} tips  ({pct:.1f}%)")
        print(f"{'='*56}\n")


# ======================================================================
# VTK SOIL EXPORT
# ======================================================================
def export_soil_vtk(soil: LayeredSoil, step: int, out_dir: str):
    cfg   = soil.cfg
    fname = os.path.join(out_dir, f"soil_moisture_{step:03d}.vtk")
    dx    = (cfg.xmax - cfg.xmin) / (cfg.res_x - 1)
    dy    = (cfg.ymax - cfg.ymin) / (cfg.res_y - 1)
    dz    = (cfg.zmax - cfg.zmin) / (cfg.res_z - 1)
    with open(fname, "w") as f:
        f.write("# vtk DataFile Version 3.0\n")
        f.write(f"Layered soil moisture — step {step}\n")
        f.write("ASCII\nDATASET STRUCTURED_POINTS\n")
        f.write(f"DIMENSIONS {cfg.res_x} {cfg.res_y} {cfg.res_z}\n")
        f.write(f"ORIGIN {cfg.xmin} {cfg.ymin} {cfg.zmin}\n")
        f.write(f"SPACING {dx:.4f} {dy:.4f} {dz:.4f}\n")
        f.write(f"POINT_DATA {cfg.res_x * cfg.res_y * cfg.res_z}\n")
        f.write("SCALARS theta float 1\nLOOKUP_TABLE default\n")
        for k in range(cfg.res_z):
            for j in range(cfg.res_y):
                for i in range(cfg.res_x):
                    f.write(f"{soil.grid[i, j, k]:.4f}\n")


# ======================================================================
# CSV TIP EXPORT
# ======================================================================
def export_tips_csv(tracker: TipTracker, out_dir: str):
    """
    Save all tip trajectories to CSV.
    Columns: root_id, x, y, z, t, layer
    The 'layer' column (Loam/Sand/Clay) shows which soil layer
    each tip was in at each timestep.
    """
    path = os.path.join(out_dir, "tip_trajectories.csv")
    with open(path, "w") as f:
        f.write("root_id,x,y,z,t,layer\n")
        for rid, entries in tracker.tip_trajectories.items():
            for x, y, z, t in entries:
                label = get_layer(z)["label"]
                f.write(f"{rid},{x:.3f},{y:.3f},{z:.3f},{t:.1f},{label}\n")
    print(f"Tip trajectories saved → {path}")


# ======================================================================
# MAIN PREDICTOR
# ======================================================================
class RootGrowthPredictor:
    """
    Runs the full Picea Abies root growth simulation in layered soil.

    Parameters
    ----------
    xml_path : str        — path to Picea Abies XML parameter file
    cfg      : SimConfig  — optional config (uses defaults if not given)
    """

    def __init__(self, xml_path: str, cfg: Optional[SimConfig] = None):
        self.xml_path = xml_path
        self.cfg      = cfg or SimConfig()
        os.makedirs(self.cfg.output_dir, exist_ok=True)

    # ------------------------------------------------------------------
    def _validate_xml(self):
        if not os.path.isfile(self.xml_path):
            raise FileNotFoundError(
                f"\n[ERROR] XML parameter file not found.\n"
                f"  Path tried : {self.xml_path}\n"
                f"  CWD        : {os.getcwd()}\n"
                f"  Fix        : set XML_PATH at the top of this script "
                f"to the correct absolute path."
            )
        print(f"XML found : {os.path.abspath(self.xml_path)}")

    # ------------------------------------------------------------------
    def _build_obstacles(self):
        cfg   = self.cfg
        parts = []

        if cfg.rock_min and cfg.rock_max:
            rock = pb.SDF_Cuboid(pb.Vector3d(*cfg.rock_min),
                                 pb.Vector3d(*cfg.rock_max))
            rock = pb.SDF_RotateTranslate(
                rock, cfg.rock_rotation_deg, 2, pb.Vector3d(0, 0, 0))
            parts.append(rock)

        if cfg.depth_limit is not None:
            z     = cfg.depth_limit
            floor = pb.SDF_HalfPlane(
                pb.Vector3d(-1000, -1000, z),
                pb.Vector3d( 1000, -1000, z),
                pb.Vector3d(-1000,  1000, z))
            parts.append(floor)

        if not parts:
            return None
        obstacle = parts[0]
        for p in parts[1:]:
            obstacle = pb.SDF_Union(obstacle, p)
        return obstacle

    # ------------------------------------------------------------------
    def run(self) -> dict:
        """
        Run the simulation.

        Returns
        -------
        dict:
            tip_trajectories : {root_id: [(x,y,z,t), ...]}
            tip_positions    : {t:       [(root_id, x, y, z), ...]}
            final_tips       : [(root_id, x, y, z)] at last timestep
            tracker          : TipTracker instance
            soil             : LayeredSoil instance
        """
        self._validate_xml()
        cfg = self.cfg

        print("\n" + "="*56)
        print("  LAYERED SOIL ROOT GROWTH SIMULATION")
        print("="*56)
        for layer in LAYERS:
            print(f"  {layer['z_top']:+5.0f} → {layer['z_bottom']:+5.0f} cm : "
                  f"{layer['label']:<5}  "
                  f"θs={layer['theta_s']}  "
                  f"α={layer['alpha']}  "
                  f"n={layer['n']}  "
                  f"z_nappe={layer['z_nappe']:+.0f} cm")
        print("="*56)

        # --- Geometry ---
        obstacle  = self._build_obstacles()
        if obstacle:
            vp.write_container(
                obstacle,
                os.path.join(cfg.output_dir, "obstacles.vtp"))

        domain    = pb.SDF_PlantBox(2000, 2000, 2000)
        navigable = pb.SDF_Difference(domain, obstacle) if obstacle else domain

        # --- Layered soil ---
        soil = LayeredSoil(cfg, obstacle=obstacle)

        # --- Plant ---
        plant = pb.Plant()
        plant.readParameters(self.xml_path)
        plant.setSoil(soil)
        plant.setGeometry(navigable)
        plant.initialize()

        # --- Tropisms ---
        # t_main: structural roots (Sinkers subType 3, Fine Roots subType 2)
        #         moderate gravity pull + hydrotropism
        t_main = TropismeMixte(
            plant, cfg.n_trials_main, cfg.sigma_main, soil,
            cfg.w_grav_main, cfg.w_water_main, cfg.anoxia_threshold)

        # t_fine: short sinker roots (subType 5)
        #         very weak gravity, strong hydrotropism → follows moisture closely
        t_fine = TropismeMixte(
            plant, cfg.n_trials_fine, cfg.sigma_fine, soil,
            cfg.w_grav_fine, cfg.w_water_fine, cfg.anoxia_threshold)

        t_main.setGeometry(navigable)
        t_fine.setGeometry(navigable)

        # subType 2 = Fine Roots, 3 = Sinkers, 5 = Short Sinkers
        plant.setTropism(t_main, pb.OrganTypes.root, 2)
        plant.setTropism(t_main, pb.OrganTypes.root, 3)
        plant.setTropism(t_fine, pb.OrganTypes.root, 5)

        # --- Simulation loop ---
        tracker = TipTracker()
        n_steps = round(cfg.sim_time / cfg.dt)
        print(f"\nRunning {n_steps} steps × {cfg.dt} d = {cfg.sim_time} d …\n")

        for i in range(n_steps):
            t = (i + 1) * cfg.dt
            plant.simulate(cfg.dt)
            soil.absorb(plant, cfg.dt)
            tracker.record(plant, t)

            if cfg.export_vtp:
                plant.write(os.path.join(cfg.output_dir, f"roots_{i:03d}.vtp"))
            if cfg.export_vtk:
                export_soil_vtk(soil, i, cfg.output_dir)

            if (i + 1) % 10 == 0:
                n_tips = len(tracker.tip_positions.get(t, []))
                print(f"  step {i+1:4d}/{n_steps} | t = {t:5d} d | "
                      f"active tips = {n_tips}")

        tracker.summary()
        export_tips_csv(tracker, cfg.output_dir)

        return {
            "tip_trajectories" : tracker.tip_trajectories,
            "tip_positions"    : tracker.tip_positions,
            "final_tips"       : tracker.final_tips,
            "tracker"          : tracker,
            "soil"             : soil,
        }


# ======================================================================
# CONVENIENCE FUNCTIONS
# ======================================================================
def tips_as_array(results: dict, timestep: float) -> np.ndarray:
    """Return tip positions at a given timestep as (N, 3) numpy array [x,y,z]."""
    entries = results["tip_positions"].get(timestep, [])
    if not entries:
        raise KeyError(f"No data for timestep={timestep}.")
    return np.array([[x, y, z] for _, x, y, z in entries])


def trajectory_as_array(results: dict, root_id: int) -> np.ndarray:
    """Return full trajectory of one root as (T, 4) array [x, y, z, t]."""
    traj = results["tip_trajectories"].get(root_id)
    if traj is None:
        raise KeyError(f"Root id {root_id} not found.")
    return np.array(traj)


def deepest_tips(results: dict, n: int = 10) -> list:
    """Return the n deepest final tip positions as (root_id, x, y, z)."""
    return sorted(results["final_tips"], key=lambda r: r[3])[:n]


# ======================================================================
# ENTRY POINT
# ======================================================================
if __name__ == "__main__":

    predictor = RootGrowthPredictor(XML_PATH)
    results   = predictor.run()

    # --- Quick Python inspection ---
    final = results["final_tips"]
    print(f"Total tips at end of simulation : {len(final)}")

    print("\nTop 5 deepest root tips:")
    for rid, x, y, z in deepest_tips(results, n=5):
        label = get_layer(z)["label"]
        print(f"  Root {rid:4d}  ({x:+7.1f}, {y:+7.1f}, {z:+7.1f}) cm  [{label}]")


