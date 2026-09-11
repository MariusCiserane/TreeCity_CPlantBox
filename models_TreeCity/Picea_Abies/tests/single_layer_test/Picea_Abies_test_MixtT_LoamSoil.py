import plantbox as pb
import plantbox.visualisation.vtk_plot as vp
import plantbox.functional.van_genuchten as vg
import numpy as np
import random

#gestion des chemins
import os
from pathlib import Path

# ====================================================================
# 1. LE SOL DYNAMIQUE (Extrait de votre script LoamSoil)
# ====================================================================
class SolDynamique(pb.SoilLookUp):
    """
    Représente l'environnement physique et hydrologique du sol.\n
    Gère la diffusion de l'eau (modèle de Van Genuchten), l'assèchement par les racines,
    et la disponibilité en oxygène des pores du sol.
    """
    def __init__(self, obstacle=None):
        """
        Initialise la grille 3D du sol, les paramètres de Van Genuchten pour un sol limoneux (Loam),
        et génère des poches d'humidité aléatoires au-dessus de la nappe phréatique.
        
        Args:
            obstacle (SDF, optionnel): Objet géométrique représentant la roche/limite. Défaut à None.
        """
        super().__init__()
        self.obstacle = obstacle
        self.theta_s = 0.43
        self.loam = vg.Parameters([0.08, self.theta_s, 0.04, 1.6, 5.0])
        self.z_nappe = -120

        self.res_x, self.res_y, self.res_z = 40, 40, 40
        self.xmin, self.xmax = -300.0, 300.0
        self.ymin, self.ymax = -300.0, 300.0
        self.zmin, self.zmax = -200.0, 5.0
        
        self.grid = np.zeros((self.res_x, self.res_y, self.res_z))

        random.seed(48)
        self.pics_aleatoires = []
        for _ in range(30): 
            px = random.uniform(self.xmin, self.xmax)
            py = random.uniform(self.ymin, self.ymax)
            pz = random.uniform(self.zmax - 100, self.zmax - 10) 
            rayon = random.uniform(25, 50)   
            force = random.uniform(40, 120)  
            self.pics_aleatoires.append((px, py, pz, rayon, force))

        print("Remplissage initial du sol dynamique...")
        for i in range(self.res_x):
            x = self.xmin + i * (self.xmax - self.xmin) / (self.res_x - 1)
            for j in range(self.res_y):
                y = self.ymin + j * (self.ymax - self.ymin) / (self.res_y - 1)
                for k in range(self.res_z):
                    z = self.zmin + k * (self.zmax - self.zmin) / (self.res_z - 1)

                    if self.obstacle is not None and self.obstacle.getDist(pb.Vector3d(x, y, z)) < 0:
                        self.grid[i, j, k] = 0.08 
                    else:
                        h_base = self.z_nappe - z - 20 
                        bonus_pics = 0
                        for px, py, pz, rayon, force in self.pics_aleatoires:
                            dist = np.sqrt((x - px)**2 + (y - py)**2 + (z - pz)**2)
                            influence = force * np.exp(-dist/50)
                            if influence > bonus_pics:
                                bonus_pics = influence

                        h_final = h_base + bonus_pics
                        if z > -30:
                            h_final -= (z + 30) * 5 
                        h_final = min(0.0, h_final)
                        self.grid[i, j, k] = vg.water_content(h_final, self.loam)
        self.grid_initial = np.copy(self.grid)

    def _get_indices(self, pos):
        """
        Convertit une position spatiale 3D continue en indices discrets pour la grille matricielle (i, j, k).
        
        Args:
            pos (pb.Vector3d): Coordonnées (x, y, z) du point à analyser.
            
        Returns:
            tuple: Les indices (i, j, k) correspondants dans le tableau self.grid, bornés aux limites du domaine.
        """
        i = int((pos.x - self.xmin) / (self.xmax - self.xmin) * (self.res_x - 1))
        j = int((pos.y - self.ymin) / (self.ymax - self.ymin) * (self.res_y - 1))
        k = int((pos.z - self.zmin) / (self.zmax - self.zmin) * (self.res_z - 1))
        return max(0, min(self.res_x - 1, i)), max(0, min(self.res_y - 1, j)), max(0, min(self.res_z - 1, k))

    def getWaterContent(self, pos):
        """
        Calcule la teneur en eau volumique (theta) réelle et physique à un point précis.
        Prend en compte la nappe phréatique, les poches aléatoires, l'évaporation de surface (Z > -30)
        et l'eau préalablement pompée par le système racinaire.
        
        Args:
            pos (pb.Vector3d): La position spatiale à évaluer.
            
        Returns:
            float: L'humidité réelle, strictement bornée entre l'humidité résiduelle (0.08) et la saturation (0.43).
        """
        h_base = self.z_nappe - pos.z - 20 
        bonus_pics = 0
        for px, py, pz, rayon, force in self.pics_aleatoires:
            dist = np.sqrt((pos.x - px)**2 + (pos.y - py)**2 + (pos.z - pz)**2)
            influence = force * np.exp(-dist/50)
            if influence > bonus_pics:
                bonus_pics = influence

        h_final = h_base + bonus_pics

        #assèchement de la surface
        if pos.z > -30:
            h_final -= (pos.z + 30) * 5
        h_final = min(0.0, h_final)

        theta_theorique = vg.water_content(h_final, self.loam)

        i, j, k = self._get_indices(pos)
        eau_pompee = self.grid_initial[i, j, k] - self.grid[i, j, k]
        
        return max(0.08, theta_theorique - eau_pompee)
    
    def getOxygen(self, pos):
        """
        Calcule la fraction de volume d'air disponible dans les pores du sol pour la respiration.
        
        Args:
            pos (pb.Vector3d): La position spatiale à évaluer.
            
        Returns:
            float: Le pourcentage d'air disponible (Saturation totale - Humidité actuelle).
        """
        return self.theta_s - self.getWaterContent(pos) 

    def getValue(self, pos, organ=None):
        """
        Méthode native requise par pb.SoilLookUp. Sert de "radar" pour le moteur C++ (Hydrotropism).
        Doit impérativement renvoyer une valeur comprise dans l'intervalle [0, 1].
        
        Args:
            pos (pb.Vector3d): La position sondée par le moteur.
            organ (Organ, optionnel): L'organe effectuant la requête.
            
        Returns:
            float: La teneur en eau à cette position (comprise entre 0.08 et 0.43).
        """        
        return self.getWaterContent(pos)

    def pomper_eau(self, plant, dt, taux_absorption=0.0001, seuil_anoxie=0.05):
        """
        Simule l'absorption active de l'eau par le système racinaire au fil du temps.
        Applique la règle d'asphyxie : une racine ne pompe que si elle dispose de suffisamment d'oxygène.
        
        Args:
            plant (pb.RootSystem / pb.Plant): L'arbre complet dont on récupère les nœuds.
            dt (float): Le pas de temps de la simulation (en jours ou heures selon votre modèle).
            taux_absorption (float): Quantité d'eau pompée par nœud et par unité de temps.
            seuil_anoxie (float): Pourcentage d'air minimum requis pour que la racine reste fonctionnelle.
        """
        nodes = plant.getNodes()
        voxels_occupes = {}
        for node in nodes:
            idx = self._get_indices(node)
            voxels_occupes[idx] = voxels_occupes.get(idx, 0) + 1
            
        for (i, j, k), nombre_noeuds in voxels_occupes.items():
            humidite_actuelle = self.grid[i, j, k]
            volume_air = self.theta_s - humidite_actuelle 
            if 0.08 < humidite_actuelle and volume_air >= seuil_anoxie: 
                baisse = min(taux_absorption * nombre_noeuds * dt, 0.02) 
                self.grid[i, j, k] = max(0.08, humidite_actuelle - baisse)

# ====================================================================
# 1.B FONCTION D'EXPORT DU SOL EN VTK
# ====================================================================
def exporter_sol_vtk(sol, etape):
    """
    Exporte l'état hydrologique de la grille du sol au format VTK (Structured Points).
    Permet la visualisation de la dynamique de l'eau en 3D dans le logiciel ParaView.
    
    Args:
        sol (SolDynamique): L'instance du sol contenant la matrice self.grid.
        etape (int): Le numéro de l'itération temporelle (utilisé pour nommer le fichier de sortie).
    """
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


# ====================================================================
# 2. LE CERVEAU MIXTE ADAPTÉ AU SOL RÉALISTE
# ====================================================================
class TropismeMixte(pb.Tropism):
    """
    Le "Cerveau" comportemental de la racine. 
    Combine le gravitropisme (pousse vers le bas) et l'hydrotropisme (recherche de l'eau).
    Intègre les comportements de survie biologique, notamment le plagiotropisme induit par l'hypoxie.
    """

    def __init__(self, plant, n_trials, sigma, sol_dynamique, poids_grav=0.8, poids_eau=0.2, seuil_anoxie=0.05):
        """
        Initialise le gestionnaire de tropismes avec ses sous-moteurs natifs.
        
        Args:
            plant (pb.Plant): La plante à laquelle ce tropisme s'applique.
            n_trials (float): Le nombre de directions (dés) évaluées par le moteur C++ à chaque pas.
            sigma (float): La variance angulaire des essais de trajectoire.
            sol_dynamique (SolDynamique): Le sol fournissant les données hydrologiques et d'oxygénation.
            poids_grav (float): L'importance accordée à la gravité dans le comportement normal.
            poids_eau (float): L'importance accordée à la recherche d'eau.
            seuil_anoxie (float): Le taux d'air sous lequel le comportement de fuite (plagiotropisme) s'active.
        """
        super().__init__(plant, n_trials, sigma)
        self.t_gravite = pb.Gravitropism(plant, n_trials, sigma)
        self.t_eau = pb.Hydrotropism(plant, n_trials, sigma, sol_dynamique)
        self.sol = sol_dynamique
        self.w_grav = poids_grav
        self.w_eau = poids_eau
        self.seuil_anoxie = seuil_anoxie

    def tropismObjective(self, pos, old, a, b, dx, organ=None): 
        """
        Fonction d'évaluation appelée N fois par getHeading() pour chaque direction possible testée.
        Donne une "note" au chemin futur testé (0.0 = Choix parfait, 1.0 = Pire choix absolu).
        
        Args:
            pos (pb.Vector3d): La position actuelle de l'apex de la racine.
            old (pb.Matrix3d): La matrice de rotation actuelle de la racine.
            a (float): L'angle de rotation (alpha) de la direction testée.
            b (float): L'angle de rotation (beta) de la direction testée.
            dx (float): La distance du segment de croissance (bond futur).
            organ (Organ, optionnel): L'organe en cours de croissance.
            
        Returns:
            float: Le score normalisé de cette trajectoire, strictement borné dans [0, 1].
        """
        #position future
        pos_future = pb.Tropism.getPosition(pos, old, a, b, dx)

        #capteurs normaux
        score_gravite = self.t_gravite.tropismObjective(pos, old, a, b, dx, organ)
        score_eau = self.t_eau.tropismObjective(pos, old, a, b, dx, organ)

        #sécurité pour éviter les "noeuds"
        #humidite_future = self.sol.getWaterContent(pos_future)
        # Si le sol contient plus de 35% d'eau, la racine n'a plus soif.
        # On passe temporairement le poids de l'eau à 0.0 pour éviter qu'elle ne boucle.
        #poids_eau_actuel = 0.0 if humidite_future > 0.35 else self.w_eau
        poids_eau_actuel = self.w_eau
        #score normalisé pour entrer dans l'intervalle [0,1]
        somme_poids = self.w_grav + poids_eau_actuel
        score_base = (score_gravite * self.w_grav + score_eau * poids_eau_actuel) / somme_poids

        #gestion de l'anoxie
        volume_air_futur = self.sol.getOxygen(pos_future)

        #s'il n'y a pas assez d'air, on passe à du plagiotropisme (exploration horizontale)
        if volume_air_futur < self.seuil_anoxie:
            # L'horizontale correspond à un score de gravité de 0.5.
            # On transforme ce 0.5 en 0 (Score Parfait) avec une valeur absolue !
            # Ainsi, descendre (0) ou monter (1) donneront tous les deux un mauvais score (1).
            score_base = abs(score_gravite - 0.5) * 2
        return score_base

# ====================================================================
# 3. SCRIPT PRINCIPAL : ARCHITECTURE ET SIMULATION
# ====================================================================

#---------------préparations pour l'export des résultats et l'import du xml de l'arbre-------------------
# 1. Définir le répertoire parent du script
script_dir = Path(__file__).resolve().parent
# 2. Définir le chemin vers results/BrusselSoil
dossier_simulation = script_dir / "results" / "LoamSoil_MixtTropism"
os.makedirs(dossier_simulation, exist_ok=True)
parameters_path = f"{script_dir.parent.parent.parent.parent}/modelparameter_TreeCity/structural/Picea_Abies_hydro_v2.xml"


plant = pb.Plant()
plant.readParameters(parameters_path)

print("Création de la roche et des limites...")
min_roche = pb.Vector3d(50, -150, -150) 
max_roche = pb.Vector3d(100, 100, -20)   
roche = pb.SDF_Cuboid(min_roche, max_roche)
roche = pb.SDF_RotateTranslate(roche, -45, 2, pb.Vector3d(0, 0, 0)) 

limite = pb.SDF_HalfPlane(
    pb.Vector3d(-1000, -1000, -150), 
    pb.Vector3d(1000, -1000, -150),  
    pb.Vector3d(-1000, 1000, -150)   
)

tous_les_obstacles = pb.SDF_Union(roche, limite)
domaine = pb.SDF_PlantBox(2000, 2000, 2000)
espace_navigable = pb.SDF_Difference(domaine, tous_les_obstacles)

vp.write_container(tous_les_obstacles, f"{dossier_simulation}/Picea_Abies_Roche.vtp")


# Instanciation du Sol
#sol_realiste = SolDynamique(obstacle=tous_les_obstacles)
sol_realiste = SolDynamique() 
plant.setSoil(sol_realiste)
#plant.setGeometry(espace_navigable)

plant.initialize()

# Application du Tropisme Mixte
tropisme_mixte = TropismeMixte(plant, 2.0, 1.5, sol_realiste, 0.35, 1)
tropisme_mixte2 = TropismeMixte(plant, 4.0, 0.5, sol_realiste, 0.0012, 1)

#on ajoute les obstacles au cerveau pour que les racines les évitent
#tropisme_mixte.setGeometry(espace_navigable)
#tropisme_mixte2.setGeometry(espace_navigable)

plant.setTropism(tropisme_mixte, pb.OrganTypes.root, 3)
plant.setTropism(tropisme_mixte, pb.OrganTypes.root, 5)
plant.setTropism(tropisme_mixte, pb.OrganTypes.root, 2)
plant.setTropism(tropisme_mixte2, pb.OrganTypes.root, 7)

# Boucle de Simulation
sim_time = 5000
dt = 50
n_steps = round(sim_time / dt)

print("Début de la croissance biophysique (Tropisme + Pompage)...")
for i in range(0, n_steps):
    plant.simulate(dt)
    
    # LA CLÉ EST ICI : Les racines boivent et modifient le sol en temps réel
    sol_realiste.pomper_eau(plant, dt, 10e-5)
    
    plant.write(f"{dossier_simulation}/Picea_Abies_{i:03d}.vtp")
    exporter_sol_vtk(sol_realiste, i)

print("Simulation terminée avec succès !")