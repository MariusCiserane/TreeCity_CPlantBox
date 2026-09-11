import plantbox as pb

#gestion des chemins
import os
from pathlib import Path

#---------------préparations pour l'export des résultats et l'import du xml de l'arbre-------------------
# 1. Définir le répertoire parent du script
script_dir = Path(__file__).resolve().parent
# 2. Définir le chemin vers results/BrusselSoil
dossier_simulation = script_dir / "results" / "No_Soil"
os.makedirs(dossier_simulation, exist_ok=True)
parameters_path = f"{script_dir.parent.parent.parent.parent}/modelparameter_TreeCity/structural/Picea_Abies_v2.xml"

plant = pb.Plant()
plant.readParameters(parameters_path)

# Initialize
plant.initialize()

# Simulate
sim_time = 5000
dt = 50
n_steps = round(sim_time / dt)
for i in range(0, n_steps):
    plant.simulate(dt)
    plant.write(f"{dossier_simulation}/Picea_Abies_{i}.vtp")


