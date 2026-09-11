import plantbox as pb
import numpy as np

#gestion des chemins
import os
from pathlib import Path

# --- 1. CRÉATION DE LA CLASSE DE SOL SUR MESURE ---
# Nous créons une classe qui hérite du moteur C++ pb.SoilLookUp
class MonSolMixte(pb.SoilLookUp):
    
    def __init__(self):
        # On initialise la classe mère C++
        super().__init__()
        
    # On "écrase" (override) la méthode virtuelle getValue de la documentation
    # pos est le pb.Vector3d contenant (x, y, z)
    def getValue(self, pos, organ=None):
        x, y, z = pos.x, pos.y, pos.z
        
        # 1. FORCE DE GRAVITÉ (On la met à 1.0 pour qu'elle reste la force dominante)
        force_gravite = abs(z) * 1.0 
        
        # 2. FORCE DE L'EAU (Adoucie)
        centre_poche_x, centre_poche_y, centre_poche_z = 300, 0, -80
        distance = np.sqrt((x - centre_poche_x)**2 + (y - centre_poche_y)**2 + (z - centre_poche_z)**2)
        
        # L'astuce : on ajoute "+ 50" au dénominateur. 
        # Ainsi, même si la racine est pile sur l'eau (distance = 0), 
        # la force maximale sera de 1000/50 = 20. 
        # L'eau va courber la racine, mais ne pourra jamais vaincre la gravité !
        force_eau = 1000 / (distance + 50) 
        
        return force_gravite + force_eau


# --- 2. CONFIGURATION DE L'ARBRE ---

#---------------préparations pour l'export des résultats et l'import du xml de l'arbre-------------------
# 1. Définir le répertoire parent du script
script_dir = Path(__file__).resolve().parent
# 2. Définir le chemin vers results/BrusselSoil
dossier_simulation = script_dir / "results" / "Hydro"
os.makedirs(dossier_simulation, exist_ok=True)
parameters_path = f"{script_dir.parent.parent.parent.parent}/modelparameter_TreeCity/structural/Picea_Abies_hydro_v2.xml"


plant = pb.Plant()
plant.readParameters(parameters_path)

# --- 3. AFFECTATION DU SOL ---
# On instancie notre nouvel objet "MonSolMixte"
mon_environnement = MonSolMixte()

# On donne ce sol à la plante (toujours AVANT l'initialisation)
plant.setSoil(mon_environnement)

# --- 4. INITIALISATION ET SIMULATION ---
plant.initialize()

sim_time = 5000
dt = 50
n_steps = round(sim_time / dt)

print("Début de la simulation avec le nouveau sol FSPM...")

for i in range(0, n_steps):
    plant.simulate(dt)
    
    # Enregistrement de toutes les étapes
    plant.write(f"{dossier_simulation}/Picea_Abies_hydro_{i}.vtp")

print("Simulation terminée !")