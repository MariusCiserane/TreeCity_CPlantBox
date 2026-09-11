#gestion des chemins
import os
from pathlib import Path
import sys

sys.path.append(str(Path(__file__).resolve().parent.parent.parent.parent))

from Soils import creer_sol_bruxellois_dumux, creer_sol_urbain_dumux
import plantbox as pb
import plantbox.visualisation.vtk_plot as vp
from plantbox.functional.PlantHydraulicModel import HydraulicModel_Doussan
from plantbox.functional.PlantHydraulicParameters import PlantHydraulicParameters
import plantbox.functional.van_genuchten as vg
import numpy as np

# ====================================================================
# 1. LE CAPTEUR POUR LE TROPISME (Beaucoup plus simple et rapide !)
# ====================================================================
class CapteurDumux(pb.SoilLookUp):
    """Sert uniquement de 'radar' pour l'hydrotropisme de l'arbre"""
    def __init__(self, dumux_solver):
        super().__init__()
        self.s = dumux_solver
        self.theta_array = np.array(self.s.getWaterContent())

    def maj_donnees(self):
        """Met à jour le radar à chaque pas de temps"""
        self.theta_array = np.array(self.s.getWaterContent())

    def getValue(self, pos, organ=None):
        # On utilise le picker officiel de Dumux !
        idx = self.s.pick([pos.x, pos.y, pos.z])
        if idx >= 0 and idx < len(self.theta_array):
            return self.theta_array[idx]
        return 0.08 # Humidité résiduelle si la racine sort de la boîte

    def getOxygen(self, pos):
        # On gère l'anoxie selon l'altitude géologique comme avant
        z_taw_m = 60.0 + (pos.z / 100.0)
        theta_s = 0.430 if z_taw_m > 59.0 else 0.390 
        return theta_s - self.getValue(pos)

# ====================================================================
# 2. LE TROPISME MIXTE
# ====================================================================
class TropismeMixte(pb.Tropism):
    def __init__(self, plant, n_trials, sigma, capteur_sol, poids_grav=0.8, poids_eau=0.2, seuil_anoxie=0.05):
        super().__init__(plant, n_trials, sigma)
        self.t_gravite = pb.Gravitropism(plant, n_trials, sigma)
        self.t_eau = pb.Hydrotropism(plant, n_trials, sigma, capteur_sol)
        self.sol = capteur_sol
        self.w_grav = poids_grav
        self.w_eau = poids_eau
        self.seuil_anoxie = seuil_anoxie

    def tropismObjective(self, pos, old, a, b, dx, organ=None): 
        pos_future = pb.Tropism.getPosition(pos, old, a, b, dx)
        score_gravite = self.t_gravite.tropismObjective(pos, old, a, b, dx, organ)
        score_eau = self.t_eau.tropismObjective(pos, old, a, b, dx, organ)
        
        score_base = (score_gravite * self.w_grav + score_eau * self.w_eau) / (self.w_grav + self.w_eau)
        
        volume_air_futur = self.sol.getOxygen(pos_future)
        if volume_air_futur < self.seuil_anoxie:
            score_base = abs(score_gravite - 0.5) * 2
        return score_base


# ====================================================================
# Fonctions auxiliaires pour la transpiration potentielle
# ====================================================================
# 1. Fonction pour le cycle du SOLEIL (Variation sur la journée)
def sinusoidal(t):
    # t est le temps en jours (ex: 0.5 = midi, 1.0 = minuit le lendemain)
    # Vaut 0 la nuit, monte à 2 à midi. Moyenne = 1.
    return max(0.0, np.sin(2.0 * np.pi * np.array(t) - 0.5 * np.pi) + 1.0)

# 2. Fonction pour l'ÂGE DE L'ARBRE (Croissance foliaire)
def soif_selon_age(age_jours):
    """
    Calcule le besoin journalier (mL) en fonction de l'âge de l'arbre.
    Exemple linéaire : 
    - À 0 jour : 1 mL (graine/germe)
    - À 365 jours : 50 mL (jeune sapin)
    """
    besoin_min = 1.0
    besoin_a_1_an = 50.0
    
    # Règle de 3 mathématique (croissance linéaire)
    # Tu peux changer cette formule pour une courbe exponentielle si tu as des données biblio !
    besoin_actuel = besoin_min + (age_jours / 365.0) * (besoin_a_1_an - besoin_min)
    return besoin_actuel
    
# ====================================================================
# 3. SCRIPT PRINCIPAL
# ====================================================================
#---------------préparations pour l'export des résultats et l'import du xml de l'arbre-------------------
# 1. Définir le répertoire parent du script
script_dir = Path(__file__).resolve().parent
# 2. Définir le chemin vers results/BrusselSoil
dossier_simulation = script_dir / "results" / "UrbanSoil_Dumux_2"
struct_parameters_path = f"{script_dir.parent.parent.parent.parent}/modelparameter_TreeCity/structural/Picea_Abies_hydro_v2.xml"
hydraulic_parameters_path = f"{script_dir.parent.parent.parent.parent}/modelparameter_TreeCity/functional/plant_hydraulics/Picea_Abies"
os.makedirs(dossier_simulation, exist_ok=True)

# A. Initialisation du sol Bruxellois via Dumux
box_min = [-300.0, -300.0, -500.0]  # [cm]
box_max = [300.0, 300.0, 0.0]
cell_number = [20, 20, 25]

s = creer_sol_urbain_dumux(box_min, box_max, cell_number)

# B. L'arbre "Cartographié" (Le secret pour le lier au sol)
plant = pb.MappedPlant()
plant.enableExtraNode()
plant.readParameters(struct_parameters_path)

# La fonction de lien native (Picker)
def picker(x, y, z):
    return s.pick([x, y, z])

plant.setSoilGrid(picker)
plant.initialize(True) # Le "True" est important ici pour l'activer !

# C. Modèle Hydraulique (Le vrai pompage de sève)
params = PlantHydraulicParameters()
# ATTENTION : Il faut vérifier les valeurs dans le fichier de paramètres
# état actuel : produit par IA, non vérifié
#ce read_parameter rajoute de lui-même l'extension .json
params.read_parameters(hydraulic_parameters_path)
hm = HydraulicModel_Doussan(plant, params)
hm.wilting_point = -15000

# D. Configuration du comportement (Tropisme)
capteur = CapteurDumux(s)
tropisme_vertical = TropismeMixte(plant, n_trials=3, sigma=0.5, capteur_sol=capteur, poids_grav=0.8, poids_eau=0.2, seuil_anoxie=0.05)
tropisme_fin = TropismeMixte(plant, n_trials=2, sigma=1.0, capteur_sol=capteur, poids_grav=0.0, poids_eau=1.0, seuil_anoxie=0.05)

# Application aux différents types d'organes racinaires (subTypes de l'XML)
# Les racines descendantes (Taproot=1, Sinkers=3, Short_Sinkers=5, SubSinkers=7)
plant.setTropism(tropisme_vertical, pb.OrganTypes.root, 1)
plant.setTropism(tropisme_vertical, pb.OrganTypes.root, 3)
plant.setTropism(tropisme_vertical, pb.OrganTypes.root, 5)
plant.setTropism(tropisme_vertical, pb.OrganTypes.root, 7)
# Les racines fines (Fine_Roots=2) -> Recherchent l'eau sans subir la gravité
plant.setTropism(tropisme_fin, pb.OrganTypes.root, 2)

# --- LE VACCIN ANTI DIV/0 : LA PRÉ-CROISSANCE ---
plant_age = 10.0 # On laisse pousser l'arbre 10 jours "à vide"
plant.simulate(plant_age, True)

# E. Boucle de Simulation Couplée
sim_time = 3000
dt = 1  #pas de temps en jours
n_steps = round(sim_time / dt)
#t_pot = 1 # Transpiration potentielle (mL / jour)

intervalle_export = 10  # Exporter tous les 10 pas de temps

print("Début de la simulation couplée avec Dumux...")
t = 0.0
for i in range(n_steps):
    
    # 1. Mise à jour du capteur pour le tropisme
    capteur.maj_donnees()
    
    # 2. L'arbre pousse et se dirige
    plant.simulate(dt)

    real_age = plant_age + t
    t_pot = soif_selon_age(real_age)
    
    # 3. Calcul biophysique des pressions (Terre vs Sève)
    h_s = s.getSolutionHead()
    h_x = hm.solve(real_age, -t_pot, h_s, cells=True) 
    
    # 4. L'arbre pompe l'eau et dit à Dumux combien il a pris
    fluxes = hm.soil_fluxes(real_age, h_x, h_s)

    # Bouclier de sécurité pour éviter les NaN
    fluxes = np.nan_to_num(fluxes, nan=0.0)

    s.setSource(fluxes)
    s.solve(dt) 
    
    # 5. Export des résultats officiels
    
    print(f"Étape {i}/{n_steps}")
    if(i % intervalle_export == 0):
        nom_sol = f"{dossier_simulation}/Sol_Urbain_{i//intervalle_export:03d}"
        nom_arbre = f"{dossier_simulation}/Arbre_Picea_{i//intervalle_export:03d}"

        # Dumux a des fonctions natives pour écrire directement les bons fichiers 3D
        print(f"  - Export du sol avec Dumux...")
        vp.write_soil(nom_sol, s, box_min, box_max, cell_number)
        print(f"  - Export de l'arbre avec PlantBox...")
        vp.write_plant(nom_arbre, hm.ms.plant())
        # --- CORRECTION POUR PARAVIEW ---
        # On renomme instantanément le sol en .vti pour que ParaView le comprenne
        import os
        if os.path.exists(f"{nom_sol}.vtu"):
            os.rename(f"{nom_sol}.vtu", f"{nom_sol}.vti")
        
    t += dt

print("Simulation terminée !")