import plantbox as pb
import plantbox.visualisation.vtk_plot as vp
import plantbox.functional.van_genuchten as vg
import numpy as np
import random

#gestion des chemins
import os
from pathlib import Path

# --- 1. LE SOL DYNAMIQUE 3D ---
class SolDynamique(pb.SoilLookUp):
    def __init__(self, obstacle=None):
        super().__init__()
        self.obstacle = obstacle
        
        # Paramètres d'un sol Limoneux (Loam)
        # theta_residuel (0.08) = sol sec de chez sec. theta_saturé (0.43) = boue.
        self.loam = vg.Parameters([0.08, 0.43, 0.04, 1.6, 5.0])
        self.z_nappe = -150.0

        # --- Création de la grille 3D (Voxels) ---
        self.res_x, self.res_y, self.res_z = 40, 40, 40
        self.xmin, self.xmax = -300.0, 300.0
        self.ymin, self.ymax = -300.0, 300.0
        self.zmin, self.zmax = -200.0, 5.0
        
        # Matrice qui va stocker l'humidité de chaque petit cube
        self.grid = np.zeros((self.res_x, self.res_y, self.res_z))

        random.seed(44) # On fixe la graine pour avoir les mêmes poches à chaque test !
        self.pics_aleatoires = []
        for _ in range(70):  # On sème 40 poches d'eau dans la terre
            px = random.uniform(self.xmin, self.xmax)
            py = random.uniform(self.ymin, self.ymax)
            pz = random.uniform(self.zmax - 100, self.zmax - 10) # Entre la surface et la nappe
            rayon = random.uniform(10, 35)   # Taille de la poche (entre 10 et 35 cm)
            force = random.uniform(40, 120)  # Puissance d'humidification (hPa)
            self.pics_aleatoires.append((px, py, pz, rayon, force))

        print("Remplissage initial hétérogène du sol...")
        for i in range(self.res_x):
            x = self.xmin + i * (self.xmax - self.xmin) / (self.res_x - 1)
            for j in range(self.res_y):
                y = self.ymin + j * (self.ymax - self.ymin) / (self.res_y - 1)
                for k in range(self.res_z):
                    z = self.zmin + k * (self.zmax - self.zmin) / (self.res_z - 1)

                    # 1. Gestion de l'obstacle géométrique
                    if self.obstacle is not None and self.obstacle.getDist(pb.Vector3d(x, y, z)) < 0:
                        self.grid[i, j, k] = 0.08 
                    else:
                        # 2. Le Gradient Vertical (On l'assèche un tout petit peu pour faire ressortir les flaques)
                        h_base = self.z_nappe - z - 20 
                        
                        # 3. Calcul de l'influence des pics d'eau
                        bonus_pics = 0
                        for px, py, pz, rayon, force in self.pics_aleatoires:
                            dist = np.sqrt((x - px)**2 + (y - py)**2 + (z - pz)**2)
                            #if dist < rayon:
                                # Plus on est près du centre de la poche, plus le bonus d'eau est fort
                            #    influence = force * (1.0 - (dist / rayon))
                            #    if influence > bonus_pics:
                            #        bonus_pics = influence

                            influence = force * np.exp(-dist/50)
                            if influence > bonus_pics:
                                bonus_pics = influence

                        # 4. On combine le sol de base et le bonus des pics d'eau
                        h_final = h_base + bonus_pics
                        
                        # 5. Évaporation de surface
                        if z > -30:
                            h_final -= (z + 30) * 5 
                            
                        # 6. Sécurité mathématique (Plafond de saturation)
                        h_final = min(0.0, h_final)
                            
                        # Conversion en vraie Teneur en Eau physique
                        self.grid[i, j, k] = vg.water_content(h_final, self.loam)
        self.grid_initial = np.copy(self.grid)



    def _get_indices(self, pos):
        """Trouve dans quel 'cube' (voxel) se trouve une position 3D"""
        i = int((pos.x - self.xmin) / (self.xmax - self.xmin) * (self.res_x - 1))
        j = int((pos.y - self.ymin) / (self.ymax - self.ymin) * (self.res_y - 1))
        k = int((pos.z - self.zmin) / (self.zmax - self.zmin) * (self.res_z - 1))
        
        # Sécurité pour ne pas sortir de la carte
        i = max(0, min(self.res_x - 1, i))
        j = max(0, min(self.res_y - 1, j))
        k = max(0, min(self.res_z - 1, k))
        return i, j, k

    def getValue(self, pos, organ=None):
        """Le Radar Haute Résolution de l'arbre"""
        if self.obstacle is not None and self.obstacle.getDist(pos) < 0:
            return 0.0
            
        # 1. Calcul mathématique CONTINU (Précision infinie)
        h_base = self.z_nappe - pos.z - 20 
        bonus_pics = 0
        for px, py, pz, rayon, force in self.pics_aleatoires:
            dist = np.sqrt((pos.x - px)**2 + (pos.y - py)**2 + (pos.z - pz)**2)
            #if dist < rayon:
            #    influence = force * (1.0 - (dist / rayon))
            #    if influence > bonus_pics:
            #        bonus_pics = influence
            influence = force * np.exp(-dist/50)
            if influence > bonus_pics:
                bonus_pics = influence

        h_final = h_base + bonus_pics
        if pos.z > -30:
            h_final -= (pos.z + 30) * 5 
        h_final = min(0.0, h_final)
        
        theta_theorique = vg.water_content(h_final, self.loam)

        # 2. Prise en compte de la SOIF (Lecture du Voxel)
        i, j, k = self._get_indices(pos)
        eau_pompee = self.grid_initial[i, j, k] - self.grid[i, j, k]
        
        # La réalité ressentie par la racine : La théorie MOINS ce qu'elle a déjà bu
        return max(0.08, theta_theorique - eau_pompee) * 1000000

    def pomper_eau(self, plant, dt, taux_absorption=0.0001):
        """LES RACINES BOIVENT ! On assèche le sol autour d'elles."""
        nodes = plant.getNodes()
        voxels_occupes = {}
        for node in nodes:
            idx = self._get_indices(node)
            voxels_occupes[idx] = voxels_occupes.get(idx, 0) + 1
            
        # 2. On retire l'eau de manière régulée
        for (i, j, k), nombre_noeuds in voxels_occupes.items():
            if 0.08 < self.grid[i, j, k] : # Si le sol n'est pas déjà sec 
                
                # La quantité bue dépend du nombre de racines et du temps
                baisse = (taux_absorption * nombre_noeuds * dt)
                
                # SÉCURITÉ : Une racine ne peut pas boire plus de 2% d'un gros cube en une seule étape !
                baisse = min(baisse, 0.02) 
                
                # On applique la baisse
                self.grid[i, j, k] -= baisse
                self.grid[i, j, k] = max(0.08, self.grid[i, j, k])

# --- FONCTION D'EXPORT DU SOL ---
def exporter_sol_vtk(sol, etape):
    """Génère un fichier 3D du sol numéroté pour l'animation ParaView"""
    filename = f"{dossier_simulation}/Picea_Abies_Humidite_{etape:03d}.vtk"
    with open(filename, "w") as f:
        f.write("# vtk DataFile Version 3.0\n")
        f.write(f"Teneur en Eau Etape {etape}\n")
        f.write("ASCII\n")
        f.write("DATASET STRUCTURED_POINTS\n")
        f.write(f"DIMENSIONS {sol.res_x} {sol.res_y} {sol.res_z}\n")
        f.write(f"ORIGIN {sol.xmin} {sol.ymin} {sol.zmin}\n")

        dx = (sol.xmax - sol.xmin) / (sol.res_x - 1)
        dy = (sol.ymax - sol.ymin) / (sol.res_y - 1)
        dz = (sol.zmax - sol.zmin) / (sol.res_z - 1)
        f.write(f"SPACING {dx} {dy} {dz}\n")
        
        f.write(f"POINT_DATA {sol.res_x * sol.res_y * sol.res_z}\n")
        f.write("SCALARS theta float 1\n")
        f.write("LOOKUP_TABLE default\n")
        
        for k in range(sol.res_z):
            for j in range(sol.res_y):
                for i in range(sol.res_x):
                    f.write(f"{sol.grid[i, j, k]:.4f}\n")

# --- 2. LE SCRIPT PRINCIPAL ---

#---------------préparations pour l'export des résultats et l'import du xml de l'arbre-------------------
# 1. Définir le répertoire parent du script
script_dir = Path(__file__).resolve().parent
# 2. Définir le chemin vers results/BrusselSoil
dossier_simulation = script_dir / "results" / "LoamSoil"
os.makedirs(dossier_simulation, exist_ok=True)
parameters_path = f"{script_dir.parent.parent.parent.parent}/modelparameter_TreeCity/structural/Picea_Abies_hydro_v2.xml"


plant = pb.Plant()
plant.readParameters(parameters_path)

print("Création de l'obstacle...")
roche = pb.SDF_Cuboid(pb.Vector3d(50, -150, -150), pb.Vector3d(100, 100, -20))
roche = pb.SDF_RotateTranslate(roche, -np.pi/4, 2, pb.Vector3d(0, 0, 0)) 
domaine = pb.SDF_PlantBox(2000, 2000, 2000)

limite = pb.SDF_HalfPlane(
    pb.Vector3d(-1000, -1000, -120), # Origine (coin du plan)
    pb.Vector3d(1000, -1000, -120),  # Point 1 (direction X)
    pb.Vector3d(-1000, 1000, -120)   # Point 2 (direction Y)
)

tous_obstacles = pb.SDF_Union(roche, limite)
vp.write_container(tous_obstacles, f"{dossier_simulation}/Picea_Abies_Obstacle_LoamSoil.vtp")
espace_navigable = pb.SDF_Difference(domaine, tous_obstacles)

# Initialisation
sol_dynamique = SolDynamique()
plant.setSoil(sol_dynamique)
plant.setGeometry(espace_navigable)
plant.initialize()

# Simulation
sim_time = 5000
dt = 50
n_steps = round(sim_time / dt)

print("Début de la croissance ET du pompage...")
for i in range(0, n_steps):
    # 1. L'arbre pousse en suivant le sol actuel
    plant.simulate(dt)
    
    # 2. LE COUPLAGE MAGIQUE : L'arbre assèche le sol !
    sol_dynamique.pomper_eau(plant, dt)
    
    # 3. Exportation des fichiers numérotés
    plant.write(f"{dossier_simulation}/Picea_Abies_{i:03d}.vtp")
    
    # On exporte le sol seulement 1 fois sur 5 pour ne pas saturer le disque dur
    
    exporter_sol_vtk(sol_dynamique, i)

print("Simulation biophysique terminée !")