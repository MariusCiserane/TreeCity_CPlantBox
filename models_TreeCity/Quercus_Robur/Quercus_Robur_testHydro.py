#gestion des chemins
import os
from pathlib import Path
import sys

import math
import random

import plantbox as pb
import numpy as np

# 1. On crée notre propre classe pour définir l'humidité du sol
class MaCarteDeau(pb.SoilLookUp):
    def __init__(self):
        super().__init__()
        
        # On crée notre sol une seule fois au démarrage de la simulation !
        self.poches = []
        nombre_de_poches = 10 # Mettez-en autant que vous voulez
        
        print(f"Génération du sol avec {nombre_de_poches} poches d'eau aléatoires...")
        
        random.seed(3) # On fixe la graine pour avoir les mêmes poches à chaque test 
        for _ in range(nombre_de_poches):
            # On tire au sort les coordonnées (dans un cube de 20 mètres sur 5 mètres de fond)
            x = random.uniform(-1000.0, 1000.0)
            y = random.uniform(-1000.0, 1000.0)
            z = random.uniform(-500.0, -20.0) # L'eau ne touche pas tout à fait la surface
            
            # On tire au sort la taille et l'attractivité de la flaque
            force = random.uniform(8.0, 15.0)  # La "quantité" d'eau
            pente = random.uniform(15.0, 50.0) # Le diviseur (plus il est petit, plus la pente est raide)
            
            self.poches.append((x, y, z, force, pente))

    # Cette fonction est appelée par chaque bout de racine à chaque pas de temps
    # 'p' est la position de la racine (p.x, p.y, p.z)
    def getValue(self, p, *args):
        # Vous pouvez même récupérer le temps si vous voulez vous amuser plus tard :
        # t = args[0] si args else 0 

        # EXEMPLE A : Une "nappe d'eau" très humide en profondeur
        #if p.z < -100.: 
           #return 1.0  # Humidité maximale sous les -100 cm
            
        # Par défaut, le sol est sec
        #return 0.1

        #---------------------------------------------------------------------------------
            
        # EXEMPLE B : Une "poche d'eau" (sphère) décalée sur le côté
        # Décommentez les lignes ci-dessous pour tester :
        distance_centre1 = ((p.x - 25)**2 + (p.y - 10)**2 + (p.z + 300)**2)**0.5

        distance_centre2 = ((p.x + 25)**2 + (p.y - 10)**2 + (p.z + 300)**2)**0.5

        distance_centre3 = ((p.x + 350)**2 + (p.y - 30)**2 + (p.z + 0)**2)**0.5

        distance_centre4 = ((p.x - 350)**2 + (p.y - 30)**2 + (p.z + 0)**2)**0.5
        
        # On crée un gradient géant qui s'étale sur 1000 cm (10 mètres !)
        #humidite = 10.0 - (distance_centre / 1000.0)
        humidite1 = 10.0 - (distance_centre1 / 20.0)

        humidite2 = 10.0 - (distance_centre2 / 20.0)

        humidite3 = 10.0 - (distance_centre3 / 50.0)

        humidite4 = 30.0 - (distance_centre4 / 10.0)
        # Le sol sec de base reste à 0.1
        return max(0.1, humidite1, humidite2, humidite3, humidite4)
        #return max(0.1, humidite1, humidite2,  humidite3)

        #---------------------------------------------------------------------------------

        #EXEMPLE C : humidité qui augmente en fonction de x

        # L'eau est à droite (axe X positif).
        # On normalise pour que l'humidité passe de 0.0 (à X = -1500) à 1.0 (à X = 1500)
        #valeur_gradient = (p.x + 2000.0) / 50.0 # 50.0 correspond à la "taille" du gradient : plus elle est petite, plus l'humidité change rapidement en fonction de x 
        
        #gradient symétrique : plus on s'éloigne, plus le sol est humide, que ce soit à gauche ou à droite
        #valeur_gradient = 0
        #if p.x>0:
            #valeur_gradient = p.x / 50
        #if p.x<-0:
            #valeur_gradient = -p.x / 50
        
        
        # On s'assure juste que ce ne soit pas négatif, on ne met pas de plafond pour que les racines soient toujours attirées vers la droite, même si elles sont déjà très à droite (valeur_gradient peut dépasser 1.0, mais ce n'est pas grave)
        #return max(0.0, valeur_gradient)

        #----------------------------------------------------------------------------------

        # EXEMPLE D : Une nappe horizontale vaste et inégale (Patchwork) --> pour l'instant, pas très concluant
        
        # 1. LA PROFONDEUR : On crée une "couche" centrée à -200 cm de profondeur.
        # L'épaisseur est de 100 cm : l'eau sera présente entre -100 et -300 cm.
        #epaisseur = 500.0
        # facteur_z vaut 1.0 exactement à -200 cm, et tombe à 0.0 quand on s'éloigne
        #facteur_z = max(0.0, 1.0 - abs(p.z + 200.0) / epaisseur)
        
        # 2. L'INÉGALITÉ HORIZONTALE : On utilise des vagues croisées
        # Le '/ 200.0' et '/ 300.0' définissent la taille des "flaques" géantes.
        # Plus le chiffre est grand, plus les flaques sont vastes.
        #onde_x = math.sin(p.x / 50.0 + math.cos(p.z / 100.0)+math.pi/2)  # On ajoute une petite variation en fonction de la profondeur pour plus de réalisme
        #onde_y = math.cos(p.y / 50.0 + math.pi/2)
        
        # L'addition des ondes donne un chiffre entre -2 et +2. 
        # On le transforme pour qu'il soit entre 0.0 (très sec) et 1.0 (très humide)
        #facteur_horizontal = (onde_x + onde_y + 2.0) / 4.0
        
        # 3. LE RÉSULTAT : On multiplie la couche de profondeur par les flaques
        # L'humidité maximale sera de 1.0 (0.1 de base + 0.9 d'eau)
        #humidite = 0.1 + (facteur_z * facteur_horizontal * 0.9)
        
        #return max(0.1, min(1.0, humidite))

        #----------------------------------------------------------------------------------

        # EXEMPLE E : L'Archipel (Poches d'eau multiples et asymétriques)
        
        # Poche 1 : Assez haute, sur la gauche
        """dist1 = ((p.x + 200)**2 + (p.y - 200)**2 + (p.z + 100)**2)**0.5
        hum1 = 30.0 - (dist1 / 15.0)
        
        # Poche 2 : Plus profonde, sur la droite
        dist2 = ((p.x - 200)**2 + (p.y + 200)**2 + (p.z + 250)**2)**0.5
        hum2 = 30.0 - (dist2 / 15.0)
        
        # Poche 3 : Loin devant (axe Y), profondeur moyenne
        dist3 = ((p.x - 200)**2 + (p.y - 200)**2 + (p.z + 150)**2)**0.5
        hum3 = 30.0 - (dist3 / 15.0)

        # Poche 4 : En profondeur, légèrement à gauche
        dist4 = ((p.x + 200)**2 + (p.y + 200)**2 + (p.z + 400)**2)**0.5
        hum4 = 30.0 - (dist4 / 15.0)
        
        # On garde l'humidité de base à 0.1
        # ATTENTION : On ne plafonne PAS à 1.0 pour garder un gradient actif jusqu'au centre !
        return max(0.1, hum1, hum2, hum3, hum4)"""
    
        #----------------------------------------------------------------------------------

        # EXEMPLE G : Le véritable champ de mines (aléatoire procédural)
        
        humidite_max = 0.1 # Le sol de base reste sec
        
        # Pour chaque millimètre de racine, on calcule la force d'attraction 
        # de TOUTES les poches d'eau générées dans le __init__
        for x, y, z, force, pente in self.poches:
            distance = ((p.x - x)**2 + (p.y - y)**2 + (p.z - z)**2)**0.5
            
            # La fameuse pente continue sans plafond !
            humidite_poche = force - (distance / pente)
            
            # On ne retient que la flaque qui attire le plus la racine à cet endroit précis
            if humidite_poche > humidite_max:
                humidite_max = humidite_poche
                
        return humidite_max

# --- INITIALISATION DE LA PLANTE ---
plant = pb.Plant()

script_dir = Path(__file__).resolve().parent
dossier_simulation = script_dir / "results/hydro" 
os.makedirs(dossier_simulation, exist_ok=True)
parameters_path = f"{script_dir.parent.parent}/modelparameter_TreeCity/structural/Quercus_Robur_hydro.xml"

plant.readParameters(parameters_path)

# ... (Ici, vous pouvez garder ou enlever vos obstacles physiques) ...

# 2. On applique notre carte d'eau à l'arbre
carte_eau = MaCarteDeau()
plant.setSoil(carte_eau)

#Ajout d'un obstacle plan à 200 cm de profondeur
domaine = pb.SDF_PlantBox(3000, 3000, 3000)

limite = pb.SDF_HalfPlane(
    pb.Vector3d(-1000, -1000, -400), # Origine (coin du plan)
    pb.Vector3d(1000, -1000, -400),  # Point 1 (direction X)
    pb.Vector3d(-1000, 1000, -400)   # Point 2 (direction Y)
)

espace_navigable = pb.SDF_Difference(domaine, limite)
plant.setGeometry(espace_navigable)

# 3. L'initialisation classique (à faire APRES setSoil)
plant.initialize()

sim_time = 3000  
dt = 30
n_steps = round(sim_time / dt)
for i in range(0, n_steps):
    plant.simulate(dt)
    plant.write(f"{dossier_simulation}/Quercus_Robur_hydro_{i}.vtp")

# --- EXPORT DE LA CARTE D'EAU POUR PARAVIEW ---
print("Génération du champ d'humidité 3D...")

# 1. On définit la zone que l'on veut "scanner" (en cm)
x_min, x_max = -1000.0, 1000.0
y_min, y_max = -1000.0, 1000.0
z_min, z_max = -1000.0, 5.0

# Résolution de la grille (50 cases par axe = 125 000 points de calcul)
res = 50 

xs = np.linspace(x_min, x_max, res)
ys = np.linspace(y_min, y_max, res)
zs = np.linspace(z_min, z_max, res)

# 2. On écrit un fichier VTK structuré
with open(f"{dossier_simulation}/Quercus_Robur_HumiditeSol.vtk", "w") as f:
    # En-tête obligatoire pour ParaView
    f.write("# vtk DataFile Version 3.0\n")
    f.write("Champ d'humidite du sol\n")
    f.write("ASCII\n")
    f.write("DATASET STRUCTURED_POINTS\n")
    f.write(f"DIMENSIONS {res} {res} {res}\n")
    f.write(f"ORIGIN {x_min} {y_min} {z_min}\n")
    
    # Calcul de l'espacement entre chaque point
    dx = (x_max - x_min) / (res - 1)
    dy = (y_max - y_min) / (res - 1)
    dz = (z_max - z_min) / (res - 1)
    f.write(f"SPACING {dx} {dy} {dz}\n")
    
    f.write(f"POINT_DATA {res**3}\n")
    f.write("SCALARS humidite float 1\n")
    f.write("LOOKUP_TABLE default\n")
    
    # 3. On calcule l'humidité pour chaque coordonnée
    # VTK a besoin des données dans l'ordre X, puis Y, puis Z
    for z in zs:
        for y in ys:
            for x in xs:
                # On utilise votre propre fonction pour calculer l'humidité ici !
                p = pb.Vector3d(x, y, z)
                valeur = carte_eau.getValue(p, 0) # 0 correspond au temps simulé
                f.write(f"{valeur:.4f}\n")
                
print("Export terminé ! Fichier HumiditeSol.vtk créé.")