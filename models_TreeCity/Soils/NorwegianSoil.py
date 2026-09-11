import plantbox as pb
import plantbox.functional.van_genuchten as vg
import numpy as np

class NorwegianSoil(pb.SoilLookUp):
    """
    Modèle hydrostatique d'un sol typique de forêt boréale norvégienne pour Picea Abies.
    Basé sur les données de la thèse de Moritz Shore (Plot 1 - Site de Référence).
    Les sols boréaux sont très superficiels, souvent riches en humus en surface 
    et reposant sur la roche mère à faible profondeur.
    """
    def __init__(self):
        super().__init__()
        
        # Profondeur de la roche mère très peu profonde (ex: 50 cm max)
        # Typiquement, la médiane en Norvège est entre 9 et 30 cm de profondeur.
        self.z_bedrock = -50.0 

        # --- PARAMÈTRES VAN GENUCHTEN (Extrait du Tableau 3.1 du PDF) ---
        # Format vg.Parameters : [theta_r, theta_s, alpha [1/cm], n, Ks [cm/d]]
        
        # 1. Topsoil (Couche de surface, très poreuse, forte capacité de rétention)
        # Valeurs : theta_r = 0.00, theta_s = 0.77, alpha = 7.49 m-1 (0.075 cm-1), n = 1.19
        # Ks estimé à 150 cm/jour pour un sol de surface très organique.
        self.vg_topsoil = vg.Parameters([0.00, 0.77, 0.075, 1.19, 150.0])
        
        # 2. Subsoil (Sous-sol minéral, juste avant la roche mère)
        # Valeurs : theta_r = 0.00, theta_s = 0.74, alpha = 17.9 m-1 (0.179 cm-1), n = 1.20
        # Ks estimé à 50 cm/jour (plus dense).
        self.vg_subsoil = vg.Parameters([0.00, 0.74, 0.179, 1.20, 50.0])

        # --- CONFIGURATION DE LA GRILLE 3D ---
        # Grille restreinte à la faible profondeur du sol
        self.xmin, self.xmax = -200, 200
        self.ymin, self.ymax = -200, 200
        self.zmin, self.zmax = self.z_bedrock, 0.0
        
        self.res_x, self.res_y, self.res_z = 20, 20, 10
        self.grid = np.zeros((self.res_x, self.res_y, self.res_z))
        
        self._initialiser_humidite()

    def _initialiser_humidite(self):
        """
        Initialise la teneur en eau (thêta) en fonction de la profondeur.
        En Norvège, la fonte des neiges printanière sature souvent le sol.
        """
        for i in range(self.res_x):
            for j in range(self.res_y):
                for k in range(self.res_z):
                    z = self.zmin + k * (self.zmax - self.zmin) / (self.res_z - 1)
                    
                    if z >= -20.0:
                        # Topsoil (Saturé au début du printemps)
                        self.grid[i, j, k] = 0.65 
                    else:
                        # Subsoil (Moyennement humide)
                        self.grid[i, j, k] = 0.50

    def get_soil_parameters(self, z):
        """ Retourne les paramètres Van Genuchten selon la profondeur z """
        if z >= -20.0:
            return self.vg_topsoil
        elif z >= self.z_bedrock:
            return self.vg_subsoil
        else:
            # Dans la roche mère, perméabilité quasi nulle
            return vg.Parameters([0.0, 0.05, 0.001, 1.1, 0.01])

    def getWaterContent(self, pos):
        """
        Retourne la teneur en eau à la position pos.
        Cette méthode est appelée par les moteurs de CPlantBox (via Tropism).
        """
        x, y, z = pos.x, pos.y, pos.z
        
        if z < self.z_bedrock:
            return 0.0 

        # Interpolation spatiale basique dans la grille
        idx_x = int((x - self.xmin) / (self.xmax - self.xmin) * (self.res_x - 1))
        idx_y = int((y - self.ymin) / (self.ymax - self.ymin) * (self.res_y - 1))
        idx_z = int((z - self.zmin) / (self.zmax - self.zmin) * (self.res_z - 1))
        
        idx_x = max(0, min(self.res_x - 1, idx_x))
        idx_y = max(0, min(self.res_y - 1, idx_y))
        idx_z = max(0, min(self.res_z - 1, idx_z))
        
        return self.grid[idx_x, idx_y, idx_z]

    def getOxygen(self, pos):
        """
        Retourne la fraction d'oxygène disponible (porosité remplie d'air).
        """
        z = pos.z
        param_vg = self.get_soil_parameters(z)
        porosite_totale = param_vg.theta_S
        teneur_eau = self.getWaterContent(pos)
        
        # L'espace qui n'est pas rempli d'eau est rempli d'air
        return max(0.0, porosite_totale - teneur_eau)

    def pomper_eau(self, plant, dt, taux_absorption=1e-5):
        """
        Simule l'absorption d'eau par les racines à chaque pas de temps.
        """
        nodes = plant.getNodes()
        
        # Comptabilise le nombre de nœuds racinaires par voxel
        for pos in nodes:
            z = pos.z
            if z < self.zmin or z > self.zmax: continue
            
            idx_x = int((pos.x - self.xmin) / (self.xmax - self.xmin) * (self.res_x - 1))
            idx_y = int((pos.y - self.ymin) / (self.ymax - self.ymin) * (self.res_y - 1))
            idx_z = int((z - self.zmin) / (self.zmax - self.zmin) * (self.res_z - 1))
            
            idx_x = max(0, min(self.res_x - 1, idx_x))
            idx_y = max(0, min(self.res_y - 1, idx_y))
            idx_z = max(0, min(self.res_z - 1, idx_z))
            
            param_vg = self.get_soil_parameters(z)
            theta_r = param_vg.theta_R
            humidite_actuelle = self.grid[idx_x, idx_y, idx_z]
            
            # Réduit l'humidité locale (limité par le theta résiduel)
            baisse = taux_absorption * dt
            self.grid[idx_x, idx_y, idx_z] = max(theta_r, humidite_actuelle - baisse)

    def get_espace_navigable(self):
        """
        Définit l'espace physique (SDF) propre à ce sol.
        Génère un plan impénétrable (roche mère) à la profondeur self.z_bedrock.
        """
        import plantbox as pb
        
        # 1. Le plan représentant la roche mère (infranchissable en dessous de z_bedrock)
        # On crée un immense plan horizontal à z = -50
        limite_roche = pb.SDF_HalfPlane(
            pb.Vector3d(-2000, -2000, self.z_bedrock), 
            pb.Vector3d(2000, -2000, self.z_bedrock),  
            pb.Vector3d(-2000, 2000, self.z_bedrock)   
        )
        
        # 2. Le domaine global (la "boîte" de la simulation)
        domaine = pb.SDF_PlantBox(4000, 4000, 4000)
        
        # 3. L'espace navigable est le domaine MOINS la limite de la roche
        espace_navigable = pb.SDF_Difference(domaine, limite_roche)
        
        return espace_navigable

    def get_model_name(self):
        return "NorwegianSoil"
