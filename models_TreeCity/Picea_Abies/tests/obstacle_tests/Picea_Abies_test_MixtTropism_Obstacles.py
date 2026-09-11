import plantbox as pb
import plantbox.visualisation.vtk_plot as vp
import plantbox.rsml.rsml_writer as rw
from plantbox.rsml.rsml_writer import Metadata
import numpy as np

#gestion des chemins
import os
from pathlib import Path

# --- 1. L'ENVIRONNEMENT PUR ---
class SolHumidite(pb.SoilLookUp):
    """
    Modèle de sol simplifié contenant une unique source d'eau ponctuelle (poche d'humidité).
    Permet de tester l'attraction directionnelle basique d'une racine vers un point donné.
    """

    def __init__(self, cx=300, cy=0, cz=-80):
        """
        Initialise la position de la source d'eau unique.
        
        Args:
            cx (float): Coordonnée X du centre de la poche d'eau. Défaut à 300.
            cy (float): Coordonnée Y du centre de la poche d'eau. Défaut à 0.
            cz (float): Coordonnée Z du centre de la poche d'eau. Défaut à -80.
        """
        super().__init__()
        self.cx = cx
        self.cy = cy
        self.cz = cz

    def getDistance(self, pos):
        """
        Calcule la distance euclidienne en 3D entre une position donnée et la poche d'eau.
        
        Args:
            pos (pb.Vector3d): Position actuelle à évaluer.
            
        Returns:
            float: La distance absolue (en cm) séparant la position de la source d'eau.
        """
        return np.sqrt((pos.x - self.cx)**2 + (pos.y - self.cy)**2 + (pos.z - self.cz)**2)
        
    def getValue(self, pos, organ=None):
        """
        Calcule l'attractivité (humidité relative) du sol à une position donnée.
        La fonction est normalisée de sorte qu'elle vaille exactement 1.0 au centre 
        de la source (distance = 0) et décroisse hyperboliquement avec l'éloignement.
        
        Args:
            pos (pb.Vector3d): La position sondée par le moteur C++.
            organ (Organ, optionnel): L'organe effectuant la requête.
            
        Returns:
            float: Le score d'humidité local, strictement compris dans l'intervalle (0, 1].
        """
        # Le sol ne connait que la position de l'eau
        distance = self.getDistance(pos)
        
        # Retourne uniquement une concentration/attractivité liée à l'eau
        return 50 / (distance + 50)

#un sol plus complexe avec plusieurs poches d'eau  
# --- 1. LE SOL INTELLIGENT (Le radar est caché ici !) ---
class SolHumidite2(pb.SoilLookUp):
    """
    Modèle de sol avancé contenant de multiples poches d'humidité dispersées.
    L'attractivité du sol dépend toujours de la poche la plus proche de la racine.
    """

    def __init__(self, pics_humidite=None, obstacle=None):
        """
        Initialise l'environnement avec une liste de sources d'eau et d'éventuels obstacles.
        
        Args:
            pics_humidite (list of tuples, optionnel): Liste de coordonnées (x, y, z) 
                représentant le centre de chaque poche d'eau. Si None, une configuration par défaut est utilisée.
            obstacle (SDF, optionnel): Géométrie des obstacles (non utilisée pour le calcul d'eau ici, 
                mais conservée pour l'architecture).
        """
        super().__init__()
        self.obstacle = obstacle
        if pics_humidite is None:
            self.pics = [
                (300, 0, -80), 
                (-250, 100, -120),
                (100, -200, -40), 
                (-100, 150, -180)
            ]
        else:
            self.pics = pics_humidite
            
    def getDistance(self, pos):
        """
        Trouve la poche d'eau la plus proche de la position donnée et calcule sa distance.
        
        Args:
            pos (pb.Vector3d): Position actuelle de l'apex racinaire.
            
        Returns:
            float: La distance minimale (en cm) vers la source d'eau la plus accessible.
        """
        distance_min = float('inf')
        for cx, cy, cz in self.pics:
            dist = np.sqrt((pos.x - cx)**2 + (pos.y - cy)**2 + (pos.z - cz)**2)
            if dist < distance_min:
                distance_min = dist
        return distance_min
        
    def getValue(self, pos, organ=None):
        """
        Calcule l'attractivité (humidité relative) du sol en fonction de la source la plus proche.
        
        Args:
            pos (pb.Vector3d): La position sondée par le moteur C++.
            organ (Organ, optionnel): L'organe effectuant la requête.
            
        Returns:
            float: Le score d'humidité local, normalisé dans l'intervalle (0, 1].
        """
        distance_min = self.getDistance(pos)
        return 50 / (distance_min + 50)
    
# --- 2. LE CERVEAU MIXTE ÉPURÉ ---
class TropismeMixte(pb.Tropism):
    """
    Gestionnaire comportemental qui combine le gravitropisme et l'hydrotropisme 
    avec une pondération dynamique (sécurité anti-nœuds).
    """

    def __init__(self, plant, n_trials, sigma, sol_humide, poids_grav=0.8, poids_eau=0.2):
        """
        Initialise le comportement mixte de la racine.
        
        Args:
            plant (pb.Plant): L'arbre auquel s'applique ce tropisme.
            n_trials (float): Nombre de directions générées et testées (résolution spatiale).
            sigma (float): La variance angulaire de la recherche.
            sol_humide (pb.SoilLookUp): L'objet sol servant de référence pour l'hydrotropisme.
            poids_grav (float): Poids accordé à la croissance verticale (gravité).
            poids_eau (float): Poids accordé à la recherche de la source d'eau.
        """
        super().__init__(plant, n_trials, sigma)
        self.t_gravite = pb.Gravitropism(plant, n_trials, sigma)
        self.t_eau = pb.Hydrotropism(plant, n_trials, sigma, sol_humide)

        self.sol = sol_humide
        self.w_grav = poids_grav
        self.w_eau = poids_eau

    def tropismObjective(self, pos, old, a, b, dx, organ=None):    
        """
        Évalue une trajectoire potentielle en combinant gravité et humidité.
        Intègre une sécurité comportementale : si la racine est arrivée à destination 
        (à moins de 40 cm de la source), l'hydrotropisme est désactivé pour éviter 
        les boucles infinies (pelotes racinaires).
        
        Args:
            pos (pb.Vector3d): Position actuelle de la pointe.
            old (pb.Matrix3d): Matrice de rotation actuelle.
            a (float): Angle alpha de la direction testée.
            b (float): Angle beta de la direction testée.
            dx (float): Distance du saut de croissance.
            organ (Organ, optionnel): L'organe concerné.
            
        Returns:
            float: Le score normalisé de la trajectoire, garanti dans l'intervalle [0, 1].
        """    
        # On délègue les maths lourdes au C++ !
        score_gravite = self.t_gravite.tropismObjective(pos, old, a, b, dx, organ)
        score_eau = self.t_eau.tropismObjective(pos, old, a, b, dx, organ)

        # On garde notre sécurité anti-nœuds (si on est déjà dans l'eau, on écoute la gravité)
        poids_eau = 0 if self.sol.getDistance(pos) < 40 else self.w_eau
        
        #normalisation pour que les scores soient entre 0 et 1
        somme_poids = self.w_grav + poids_eau
        return (score_gravite * self.w_grav + score_eau * poids_eau) / somme_poids

# --- 3. LE SCRIPT PRINCIPAL ---

#---------------préparations pour l'export des résultats et l'import du xml de l'arbre-------------------
# 1. Définir le répertoire parent du script
script_dir = Path(__file__).resolve().parent
# 2. Définir le chemin vers results/BrusselSoil
dossier_simulation = script_dir / "results" / "MixtTropism"
os.makedirs(dossier_simulation, exist_ok=True)
parameters_path = f"{script_dir.parent.parent.parent.parent}/modelparameter_TreeCity/structural/Picea_Abies_hydro_v2.xml"


plant = pb.Plant()
plant.readParameters(parameters_path)


# On définit notre obstacle (une grosse plaque de roche rectangulaire)
# SDF_Cuboid prend le coin minimum (x,y,z) et le coin maximum (x,y,z)
# Plaçons un mur/roche en plein sur le chemin de la poche d'eau de droite !

print("Création de l'obstacle rocheux...")
min_roche = pb.Vector3d(50, -150, -150) # Point bas-gauche
max_roche = pb.Vector3d(100, 100, -20)   # Point haut-droit
roche = pb.SDF_Cuboid(min_roche, max_roche)

roche = pb.SDF_RotateTranslate(roche, -45, 2, pb.Vector3d(0, 0, 0)) 

# (Optionnel) Vous pouvez en créer un deuxième
# roche2 = pb.SDF_Cuboid(pb.Vector3d(-150, -50, -100), pb.Vector3d(-100, 50, -50))
# obstacles_combines = pb.SDF_Union(roche, roche2)

#ajout d'un obstacle pour limiter la propagation au delà de 60cm de la surface
limite = pb.SDF_HalfPlane(
    pb.Vector3d(-1000, -1000, -60), # Origine (coin du plan)
    pb.Vector3d(1000, -1000, -60),  # Point 1 (direction X)
    pb.Vector3d(-1000, 1000, -60)   # Point 2 (direction Y)
)

tous_les_obstacles = pb.SDF_Union(roche, limite)

# On instancie notre sol et on l'assigne
#mon_sol = SolHumidite(cx=300, cy=0, cz=-80)
mon_sol2 = SolHumidite2(obstacle=tous_les_obstacles)
plant.setSoil(mon_sol2)

# On définit la limite globale de la terre (un espace immense de 2000 cm)
# Pour que les racines aient la place de s'exprimer
domaine = pb.SDF_PlantBox(2000, 2000, 2000)

# 3. On soustrait la roche du domaine navigable
espace_navigable = pb.SDF_Difference(domaine, tous_les_obstacles)

#export de l'obstacle pour visualisation dans ParaView
vp.write_container(tous_les_obstacles, "results/ObstacleRoche.vtp")

# 4. On assigne cet espace physique à l'arbre
plant.setGeometry(espace_navigable)

# Initialisation
plant.initialize()

# On instancie notre tropisme mixte
tropisme_mixte = TropismeMixte(plant, 2.0, 1.5, mon_sol2, 0.027, 1)
tropisme_mixte2 = TropismeMixte(plant, 4.0, 0.5, mon_sol2, 0.001, 1)

tropisme_mixte.setGeometry(espace_navigable)
tropisme_mixte2.setGeometry(espace_navigable)

plant.setTropism(tropisme_mixte, pb.OrganTypes.root, 3)
plant.setTropism(tropisme_mixte, pb.OrganTypes.root, 5)
plant.setTropism(tropisme_mixte, pb.OrganTypes.root, 2)
plant.setTropism(tropisme_mixte2, pb.OrganTypes.root, 7)

# Simulation
sim_time = 5000
dt = 50
n_steps = round(sim_time / dt)

#export du système
for i in range(0, n_steps):
    plant.simulate(dt)
    plant.write(f"{dossier_simulation}/Picea_Abies_Mixte_{i}.vtp")

#export de la carte d'humidité pour ParaView
print("Génération du champ d'humidité 3D...")
x_min, x_max = -1000.0, 1000.0
y_min, y_max = -1000.0, 1000.0
z_min, z_max = -1000.0, 5.0
res = 50
xs = np.linspace(x_min, x_max, res)
ys = np.linspace(y_min, y_max, res)
zs = np.linspace(z_min, z_max, res)
with open(f"{dossier_simulation}/Picea_Abies_HumiditeSol.vtk", "w") as f:
    # En-tête obligatoire pour ParaView
    f.write("# vtk DataFile Version 3.0\n")
    f.write("Champ d'humidite du sol - Tropisme Mixte\n")
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
    
    for z in zs:
        for y in ys:
            for x in xs:
                humidite = mon_sol2.getValue(pb.Vector3d(x, y, z))
                f.write(f"{humidite}\n")
print("Champ d'humidité généré !")