#gestion des chemins
import os
from pathlib import Path
import sys

import plantbox as pb
import numpy as np

script_dir = Path(__file__).resolve().parent
dossier_simulation = script_dir / "results" 
parameters_path = f"{script_dir.parent.parent}/modelparameter_TreeCity/structural/Quercus_Robur.xml"
# 1. Création de l'objet plante
pl = pb.MappedPlant()

# 2. Chargement de ton fichier XML d'amélioration
# Assure-toi que le nom correspond exactement à ton fichier
pl.readParameters(parameters_path)

# 3. Mise en place de la Nappe Phréatique (Obstacle)
# On utilise un SDF (Signed Distance Field) de type "HalfSpace" (Demi-espace).
# - Le premier vecteur [0, 0, 1] est la normale (pointe vers le ciel, axe Z positif).
# - Le deuxième vecteur [0, 0, -300] est la position de l'obstacle. 
# Ici, la nappe phréatique est fixée à -300 cm (3 mètres de profondeur).
profondeur_nappe = -250.0
# On utilise directement les vecteurs 3D de CPlantBox au lieu de numpy
nappe_phreatique = pb.SDF_HalfPlane(pb.Vector3d(0., 0., profondeur_nappe), pb.Vector3d(0., 0., -1.))

# 4. On applique cette géométrie (cet obstacle) à l'environnement de la plante
pl.setGeometry(nappe_phreatique)

# 5. Initialisation du modèle
pl.initialize()

# 6. Lancement de la simulation
# Tu peux ajuster le nombre de jours en fonction de ce que tu avais mis dans ton XML (simulationTime)
jours_de_simulation = 3000
print(f"Simulation en cours pour {jours_de_simulation} jours...")
pl.simulate(jours_de_simulation)

# 7. Exportation du résultat pour visualisation dans ParaView (.vtp)
fichier_sortie = f"{dossier_simulation}/Quercus_Robur_Avec_Obstacle.vtp"
pl.write(fichier_sortie)

print(f"Simulation terminée avec succès ! Ouvre {fichier_sortie} dans ParaView.")