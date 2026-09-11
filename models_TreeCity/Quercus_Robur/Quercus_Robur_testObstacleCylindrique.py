import plantbox as pb
import plantbox.visualisation.vtk_plot as vp


file = "../../modelparameter_TreeCity/structural/Quercus_Robur.xml"
# 1. Création de l'objet plante
pl = pb.MappedPlant()

# 2. Chargement de ton fichier XML d'amélioration
# Assure-toi que le nom correspond exactement à ton fichier
pl.readParameters(file)

# 1. Le sol global (20m x 20m x 20m)
sol = pb.SDF_PlantBox(2000., 2000., 2000.)

# On crée un cylindre : rayon haut (40), rayon bas (40), hauteur (40), carré (False)
# Cela crée un bel obstacle rond de 80 cm de diamètre
caillou_base = pb.SDF_PlantContainer(100., 100., 40., False)

# Et on le translate exactement comme avant, par exemple à -100 cm de profondeur
caillou_profond = pb.SDF_RotateTranslate(caillou_base, pb.Vector3d([0., 0., -100.]))

# 4. La soustraction géométrique
sol_avec_obstacle = pb.SDF_Difference(sol, caillou_profond)

# 5. Application à la plante
pl.setGeometry(sol_avec_obstacle)

# -------------------------------------------------

# 5. Initialisation du modèle
pl.initialize()

# 6. Lancement de la simulation
# Tu peux ajuster le nombre de jours en fonction de ce que tu avais mis dans ton XML (simulationTime)
jours_de_simulation = 3000
print(f"Simulation en cours pour {jours_de_simulation} jours...")
pl.simulate(jours_de_simulation)

# 7. Exportation du résultat pour visualisation dans ParaView (.vtp)
fichier_sortie = "results/Quercus_Robur_Avec_ObstacleCylindrique.vtp"
pl.write(fichier_sortie)
vp.write_container(sol_avec_obstacle, "results/ObstacleCylindrique.vtp")

print(f"Simulation terminée avec succès ! Ouvre {fichier_sortie} dans ParaView.")