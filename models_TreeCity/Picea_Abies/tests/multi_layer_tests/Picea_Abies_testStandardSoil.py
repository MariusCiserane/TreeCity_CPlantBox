#gestion des chemins
import os
from pathlib import Path
import sys

#importation des classes de sol
sys.path.append(str(Path(__file__).resolve().parent.parent.parent.parent))
from Soils import BrusselsSoil, NorwegianSoil, exporter_sol_vtk

#module cplantbox
import plantbox as pb
import plantbox.visualisation.vtk_plot as vp
import numpy as np

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
dossier_simulation = script_dir / "results" 
parameters_path = f"{script_dir.parent.parent.parent.parent}/modelparameter_TreeCity/structural/Picea_Abies_hydro_v2.xml"

#--------------------------------Architecture et simulation--------------------------------
# 1. Création de la plante et chargement des paramètres
plant = pb.Plant()
print(f"Chargement des paramètres de croissance depuis : {parameters_path}")
plant.readParameters(parameters_path)

# ---------------------------------------------------------
# 2. INTÉGRATION DU SOL BRUXELLOIS
# ---------------------------------------------------------
print("Connexion de l'arbre au sous-sol Bruxellois...")
# On instancie la classe
sol = BrusselsSoil()

print("Création de la roche et des limites...")
espace_navigable = sol.get_espace_navigable()

# (Garde tes paramètres d'obstacles tels quels)
min_roche = pb.Vector3d(50, -150, -150) 
max_roche = pb.Vector3d(100, 100, -20)   
roche = pb.SDF_Cuboid(min_roche, max_roche)
# Attention : np.pi/4 correspond à 45 degrés en radians. 
# Selon ce qu'attend ta version de CPlantBox, vérifie si c'est en degrés ou en radians.
roche = pb.SDF_RotateTranslate(roche, -45, 2, pb.Vector3d(0, 0, 0)) 

espace_navigable = pb.SDF_Difference(espace_navigable, roche)

# ajout du nom au chemin d'export
dossier_simulation = dossier_simulation / sol.get_model_name()
# Créer le dossier et ses parents si nécessaire
os.makedirs(dossier_simulation, exist_ok=True)
# Export optionnel de l'obstacle pour vérifier sous ParaView
vp.write_container(roche, f"{dossier_simulation}/Picea_Abies_Obstacle_Roche.vtp")

# On "plante" l'arbre dans ce sol spécifique
plant.setSoil(sol)
#plant.setSoil(sol_boreal)
plant.setGeometry(espace_navigable)
plant.initialize()

# ---------------------------------------------------------
# 3. CONFIGURATION DES TROPISMES (Le Cerveau)
# ---------------------------------------------------------
# On passe l'instance de 'sol' au TropismeMixte pour qu'il 
# puisse "sentir" l'humidité et l'oxygène des vraies couches géologiques.
tropisme_mixte = TropismeMixte(plant, 2.0, 1.5, sol, 0.35, 1)
tropisme_mixte2 = TropismeMixte(plant, 4.0, 0.5, sol, 0.0012, 1)

tropisme_mixte.setGeometry(espace_navigable)
tropisme_mixte2.setGeometry(espace_navigable)

# Application aux différents types d'organes racinaires
plant.setTropism(tropisme_mixte, pb.OrganTypes.root, 3)
plant.setTropism(tropisme_mixte, pb.OrganTypes.root, 5)
plant.setTropism(tropisme_mixte, pb.OrganTypes.root, 2)
plant.setTropism(tropisme_mixte2, pb.OrganTypes.root, 7)

# ---------------------------------------------------------
# 4. BOUCLE DE SIMULATION DYNAMIQUE
# ---------------------------------------------------------
sim_time = 5000
dt = 50
n_steps = round(sim_time / dt)

print("Début de la croissance biophysique...")
for i in range(0, n_steps):
    # L'arbre pousse d'un pas de temps
    plant.simulate(dt)
    
    # LA CLÉ : L'arbre boit et assèche le sol stratifié
    # J'ai repris ton taux de 10e-5
    sol.pomper_eau(plant, dt, 10e-5) 
    
    # Export du système racinaire
    chemin_vtp = dossier_simulation / f"Picea_Abies_{sol.get_model_name()}_{i:03d}.vtp"
    plant.write(str(chemin_vtp))
    
    # Export du sol (Astuce : on exporte 1 fois sur 5 pour ne pas exploser la taille du dossier)
    
    exporter_sol_vtk(sol, i, dossier_simulation)

print("Simulation terminée avec succès !")

