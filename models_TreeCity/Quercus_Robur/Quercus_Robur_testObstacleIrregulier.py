import plantbox as pb
import plantbox.visualisation.vtk_plot as vp
import numpy as np


file = "../../modelparameter_TreeCity/structural/Quercus_Robur.xml"
import plantbox as pb

plant = pb.Plant()
plant.readParameters(file)

# 1. Le sol global
sol = pb.SDF_PlantBox(2000., 2000., 2000.)

# 2. Création des 3 "morceaux" de la roche (Des cylindres de tailles différentes)
# Morceau A : Le coeur du caillou (épais et large)
"""partA_base = pb.SDF_PlantContainer(40., 40., 40., False)
partA = pb.SDF_RotateTranslate(partA_base, pb.Vector3d([0., 0., -100.]))

# Morceau B : Une excroissance sur la droite (plus fine, un peu plus haute)
partB_base = pb.SDF_PlantContainer(25., 25., 30., False)
partB = pb.SDF_RotateTranslate(partB_base, pb.Vector3d([35., 15., -90.]))

# Morceau C : Une bosse asymétrique sur la gauche (plus profonde)
partC_base = pb.SDF_PlantContainer(30., 30., 20., False)
partC = pb.SDF_RotateTranslate(partC_base, pb.Vector3d([-25., -20., -110.]))


partA2 = pb.SDF_PlantContainer(40., 40., 40., False)
partB2 = pb.SDF_PlantContainer(25., 25., 30., False)
partC2 = pb.SDF_PlantContainer(30., 30., 20., False)

partB2 = pb.SDF_RotateTranslate(partB2, pb.Vector3d([35., 15., 10.]))

partC2 = pb.SDF_RotateTranslate(partC2, pb.Vector3d([-25., -20., -10.]))

# 3. LA FUSION (Union mathématique)
# Note : SDF_Union ne prend que deux objets à la fois. On les fusionne donc en chaîne !
roche_etape1 = pb.SDF_Union(partA, partB)          # Fusionne A et B
roche_complete = pb.SDF_Union(roche_etape1, partC) # Ajoute C au résultat

roche_complete = pb.SDF_RotateTranslate(roche_complete, 50, 0, pb.Vector3d(0., 0., 0.))  # Rotation finale pour plus de réalisme

roche_complete2 = pb.SDF_Union(partA2, partB2)          # Fusionne A et B
roche_complete2 = pb.SDF_Union(roche_complete2, partC2) # Ajoute C au résultat
roche_complete2 = pb.SDF_RotateTranslate(roche_complete2, -50, 0, pb.Vector3d(0., 0., -100.))"""  # Rotation finale pour plus de réalisme

caillou = pb.SDF_PlantContainer(10., 10., 5., False)
caillou = pb.SDF_RotateTranslate(caillou, pb.Vector3d([0., 0., -250.]))

# 4. La soustraction géométrique (On enlève la grosse roche fusionnée du sol)
"""sol_avec_obstacle = pb.SDF_Difference(sol, roche_complete)
sol_avec_obstacle = pb.SDF_Difference(sol_avec_obstacle, roche_complete2)
sol_avec_obstacle = pb.SDF_Difference(sol_avec_obstacle, caillou)"""

sol_avec_obstacle = pb.SDF_Difference(sol, caillou)

# --- (Optionnel) Export pour vérifier la forme de la roche dans ParaView ---
#vp.write_container(caillou, "results/ObstacleIrregulier.vtp")
#script_paraview = caillou.writePVPScript()

# On le sauvegarde dans un fichier Python
#with open("results/MonPetitCaillou.py", "w") as fichier_python:
  #  fichier_python.write(script_paraview)
    
# Si vous voulez aussi exporter vos grosses roches pour les voir :
#with open("results/MaRoche1.py", "w") as f:
    #f.write(roche_complete.writePVPScript())

#with open("results/MaRoche2.py", "w") as f:
    #f.write(roche_complete2.writePVPScript())

with open("results/TousObstacles.py", "w") as f:
    f.write(sol_avec_obstacle.writePVPScript())
# -------------------------------------------------------------------------

# 5. Application à la plante et simulation
plant.setGeometry(sol_avec_obstacle)

plant.initialize()
sim_time = 3000  
dt = 30
n_steps = round(sim_time / dt)
for i in range(0, n_steps):
    plant.simulate(dt)
    plant.write("results/Quercus_Robur_testObstacleIrregulier_" + str(i) + ".vtp")