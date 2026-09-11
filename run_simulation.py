import plantbox as pb
import plantbox.visualisation.vtk_plot as vp

# 1. Instancier le système racinaire
rs = pb.RootSystem()

# 2. Charger votre fichier de paramètres XML
rs.readParameters("mon_espece.xml")

# 3. Initialiser la graine selon les paramètres du XML
rs.initialize()

# 4. Lancer la simulation temporelle (ex. 30 jours, pas de 1 jour)
rs.simulate(30, 1)

# 5. Exporter les résultats 3D
rs.write("resultat_espece.vtp")
rs.write("resultat_espece.rsml")
print("Simulation terminée. Fichiers resultat_espece.vtp et .rsml générés.")

# 6. Afficher le rendu 3D interactif
vp.plot_roots(rs, "subType")