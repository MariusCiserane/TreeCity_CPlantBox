#gestion des chemins
import os
from pathlib import Path
import sys

sys.path.append(str(Path(__file__).resolve().parent.parent.parent.parent))

import plantbox as pb
from Soils import BrusselsSoil, NorwegianSoil, exporter_sol_vtk
import plantbox.visualisation.vtk_plot as vp
import numpy as np

class TropismeMixte(pb.Tropism):
    """
    Le "Cerveau" comportemental de la racine. 
    Combine le gravitropisme (pousse vers le bas) et l'hydrotropisme (recherche de l'eau).
    Intègre les comportements de survie biologique, notamment le plagiotropisme induit par l'hypoxie.
    """

    def __init__(self, plant, n_trials, sigma, sol_dynamique, poids_grav=0.8, poids_eau=0.2, seuil_anoxie=0.05):
        super().__init__(plant, n_trials, sigma)
        self.t_gravite = pb.Gravitropism(plant, n_trials, sigma)
        self.t_eau = pb.Hydrotropism(plant, n_trials, sigma, sol_dynamique)
        self.sol = sol_dynamique
        self.w_grav = poids_grav
        self.w_eau = poids_eau
        self.seuil_anoxie = seuil_anoxie

    def tropismObjective(self, pos, old, a, b, dx, organ=None): 
        # position future
        pos_future = pb.Tropism.getPosition(pos, old, a, b, dx)

        # capteurs normaux
        score_gravite = self.t_gravite.tropismObjective(pos, old, a, b, dx, organ)
        score_eau = self.t_eau.tropismObjective(pos, old, a, b, dx, organ)

        poids_eau_actuel = self.w_eau
        somme_poids = self.w_grav + poids_eau_actuel
        
        # Si la racine est purement agéotropique (w_grav == 0), on ne divise pas par 0
        if somme_poids == 0:
            score_base = 0.5 
        else:
            score_base = (score_gravite * self.w_grav + score_eau * poids_eau_actuel) / somme_poids

        # gestion de l'anoxie : plagiotropisme de survie
        volume_air_futur = self.sol.getOxygen(pos_future)
        if volume_air_futur < self.seuil_anoxie:
            # L'horizontale correspond à un score de gravité de 0.5.
            # On transforme ce 0.5 en 0 (Score Parfait) avec une valeur absolue.
            score_base = abs(score_gravite - 0.5) * 2
            
        return score_base


# ====================================================================
# SCRIPT PRINCIPAL : ARCHITECTURE ET SIMULATION
# ====================================================================

#---------------préparations pour l'export des résultats et l'import du xml de l'arbre-------------------
# 1. Définir le répertoire parent du script
script_dir = Path(__file__).resolve().parent
# 2. Définir le chemin vers results/BrusselSoil
dossier_simulation = script_dir / "results" 
parameters_path = f"{script_dir.parent.parent.parent.parent}/modelparameter_TreeCity/structural/Picea_Abies_hydro_v3_ia.xml"

plant = pb.Plant()
plant.readParameters(parameters_path)

print("Création de la roche et des limites...")
min_roche = pb.Vector3d(50, -150, -150) 
max_roche = pb.Vector3d(100, 100, -20)   
roche = pb.SDF_Cuboid(min_roche, max_roche)
roche = pb.SDF_RotateTranslate(roche, -np.pi/4, 2, pb.Vector3d(0, 0, 0)) 

limite = pb.SDF_HalfPlane(
    pb.Vector3d(-1000, -1000, -150), 
    pb.Vector3d(1000, -1000, -150),  
    pb.Vector3d(-1000, 1000, -150)   
)

tous_les_obstacles = pb.SDF_Union(roche, limite)
domaine = pb.SDF_PlantBox(2000, 2000, 2000)
espace_navigable = pb.SDF_Difference(domaine, tous_les_obstacles)


# ---------------------------------------------------------
print("Connexion de l'arbre au sous-sol Bruxellois...")
sol = BrusselsSoil() 

# ajout du nom au chemin d'export
dossier_simulation = dossier_simulation / sol.get_model_name()
# Créer le dossier et ses parents si nécessaire
os.makedirs(dossier_simulation, exist_ok=True)
vp.write_container(tous_les_obstacles, f"{dossier_simulation}/Picea_Abies_Obstacle_Roche.vtp")

plant.setSoil(sol)
plant.setGeometry(espace_navigable)
plant.initialize()

# ---------------------------------------------------------
# CONFIGURATION DES TROPISMES BASÉE SUR LE CSV
# ---------------------------------------------------------
# 1. Tropisme pour les pivots et sinkers (Gravitropisme dominant + recherche d'eau)
tropisme_vertical = TropismeMixte(plant, n_trials=8, sigma=0.5, sol_dynamique=sol, poids_grav=0.8, poids_eau=0.2, seuil_anoxie=0.05)

# 2. Tropisme pour les racines fines (Ageotropisme + Hydro/Chemotropisme pur selon le CSV)
tropisme_fin = TropismeMixte(plant, n_trials=4, sigma=1.0, sol_dynamique=sol, poids_grav=0.0, poids_eau=1.0, seuil_anoxie=0.05)

tropisme_vertical.setGeometry(espace_navigable)
tropisme_fin.setGeometry(espace_navigable)

# Application aux différents types d'organes racinaires (subTypes de l'XML)
# Les racines descendantes (Taproot=1, Sinkers=3, Short_Sinkers=5, SubSinkers=7)
plant.setTropism(tropisme_vertical, pb.OrganTypes.root, 1)
plant.setTropism(tropisme_vertical, pb.OrganTypes.root, 3)
plant.setTropism(tropisme_vertical, pb.OrganTypes.root, 5)
plant.setTropism(tropisme_vertical, pb.OrganTypes.root, 7)

# Les racines fines (Fine_Roots=2) -> Recherchent l'eau sans subir la gravité
plant.setTropism(tropisme_fin, pb.OrganTypes.root, 2)

# ---------------------------------------------------------
# BOUCLE DE SIMULATION DYNAMIQUE
# ---------------------------------------------------------
sim_time = 5000
dt = 50
n_steps = round(sim_time / dt)

print("Début de la croissance biophysique...")
for i in range(0, n_steps):
    plant.simulate(dt)
    sol.pomper_eau(plant, dt, 10e-5) 
    
    plant.write(f"{dossier_simulation}/Picea_Abies_{sol.get_model_name()}_{i:03d}.vtp")
    
    #if i % 5 == 0:
    exporter_sol_vtk(sol, i, dossier_simulation)

print("Simulation terminée avec succès !")