import plantbox as pb
import plantbox.functional.van_genuchten as vg
import numpy as np

class BrusselsSoil(pb.SoilLookUp):
    """
    Modèle hydrostatique complet du bassin Bruxellois (9 couches géologiques).
    """
    def __init__(self):
        super().__init__()
        
        self.z_nappe = -150.0 

        # --- PARAMÈTRES VAN GENUCHTEN ---
        self.vg_quaternaire = vg.Parameters([0.078, 0.430, 0.036, 1.56, 24.96])
        self.vg_tielt       = vg.Parameters([0.089, 0.390, 0.021, 1.48, 31.40])
        self.vg_aalbeke     = vg.Parameters([0.068, 0.380, 0.008, 1.09, 4.80])
        self.vg_moen        = vg.Parameters([0.100, 0.380, 0.027, 1.23, 2.88])
        self.vg_stmaur      = vg.Parameters([0.068, 0.380, 0.008, 1.09, 4.80])
        self.vg_grandglise  = vg.Parameters([0.045, 0.430, 0.145, 2.68, 712.80])
        self.vg_lincent     = vg.Parameters([0.068, 0.380, 0.008, 1.09, 4.80])
        self.vg_cretace     = vg.Parameters([0.034, 0.460, 0.016, 1.37, 6.00])
        self.vg_paleozoique = vg.Parameters([0.010, 0.050, 0.001, 1.10, 0.01])

        # --- CRÉATION DE LA GRILLE 3D ---
        self.res_x, self.res_y, self.res_z = 40, 40, 40
        self.xmin, self.xmax = -300.0, 300.0
        self.ymin, self.ymax = -300.0, 300.0
        self.zmin, self.zmax = -300.0, 5.0 
        
        self.grid = np.zeros((self.res_x, self.res_y, self.res_z))
        
        print("Initialisation de la stratigraphie Bruxelloise...")
        self._initialiser_grille()
        self.grid_initial = np.copy(self.grid)

    def _get_proprietes_couche(self, z_cm):
        """Convertit la profondeur CPlantBox en altitude géologique."""
        z_taw_m = 60.0 + (z_cm / 100.0)

        if z_taw_m > 59.0:     
            return self.vg_quaternaire, 0.430, 0.078
        elif z_taw_m > 39.0:   
            return self.vg_tielt, 0.390, 0.089
        elif z_taw_m > 35.0:   
            return self.vg_aalbeke, 0.380, 0.068
        elif z_taw_m > 0.0:    
            return self.vg_moen, 0.380, 0.100
        elif z_taw_m > -30.0:  
            return self.vg_stmaur, 0.380, 0.068
        elif z_taw_m > -38.0:  
            return self.vg_grandglise, 0.430, 0.045
        elif z_taw_m > -58.0:  
            return self.vg_lincent, 0.380, 0.068
        elif z_taw_m > -78.0:  
            return self.vg_cretace, 0.460, 0.034
        else:                  
            return self.vg_paleozoique, 0.050, 0.010

    def _initialiser_grille(self):
        """Remplissage de la matrice basé uniquement sur l'hydrostatique."""
        for i in range(self.res_x):
            for j in range(self.res_y):
                for k in range(self.res_z):
                    z = self.zmin + k * (self.zmax - self.zmin) / (self.res_z - 1)
                    params_vg, _, _ = self._get_proprietes_couche(z)
                    
                    # Pression pure, plus besoin de vérifier l'obstacle
                    h_final = min(0.0, self.z_nappe - z)
                    self.grid[i, j, k] = vg.water_content(h_final, params_vg)

    def _get_indices(self, pos):
        i = int((pos.x - self.xmin) / (self.xmax - self.xmin) * (self.res_x - 1))
        j = int((pos.y - self.ymin) / (self.ymax - self.ymin) * (self.res_y - 1))
        k = int((pos.z - self.zmin) / (self.zmax - self.zmin) * (self.res_z - 1))
        return max(0, min(self.res_x - 1, i)), max(0, min(self.res_y - 1, j)), max(0, min(self.res_z - 1, k))

    def getWaterContent(self, pos):
        params_vg, _, theta_r = self._get_proprietes_couche(pos.z)
        h_final = min(0.0, self.z_nappe - pos.z)
        theta_theorique = vg.water_content(h_final, params_vg)

        i, j, k = self._get_indices(pos)
        eau_pompee = self.grid_initial[i, j, k] - self.grid[i, j, k]
        return max(theta_r, theta_theorique - eau_pompee)
    
    def getOxygen(self, pos):
        _, theta_s, _ = self._get_proprietes_couche(pos.z)
        return theta_s - self.getWaterContent(pos)

    def getValue(self, pos, organ=None):      
        return self.getWaterContent(pos)

    def pomper_eau(self, plant, dt, taux_absorption=0.0001, seuil_anoxie=0.05):
        nodes = plant.getNodes()
        voxels_occupes = {}
        for node in nodes:
            idx = self._get_indices(node)
            voxels_occupes[idx] = voxels_occupes.get(idx, 0) + 1
            
        for (i, j, k), nombre_noeuds in voxels_occupes.items():
            z = self.zmin + k * (self.zmax - self.zmin) / (self.res_z - 1)
            _, theta_s, theta_r = self._get_proprietes_couche(z)
            
            humidite_actuelle = self.grid[i, j, k]
            volume_air = theta_s - humidite_actuelle 
            
            if humidite_actuelle > theta_r and volume_air >= seuil_anoxie: 
                baisse = min(taux_absorption * nombre_noeuds * dt, 0.02) 
                self.grid[i, j, k] = max(theta_r, humidite_actuelle - baisse)

    def get_espace_navigable(self):
        """
        Définit l'espace physique (SDF) propre à ce sol.
        Génère un plan impénétrable (roche mère) à la profondeur self.z_bedrock.
        """
        import plantbox as pb
        
        # 1. Le plan représentant la roche mère (infranchissable en dessous de z_bedrock)
        # On crée un immense plan horizontal à z = -50
        limite_roche = pb.SDF_HalfPlane(
            pb.Vector3d(-2000, -2000, self.z_nappe), 
            pb.Vector3d(2000, -2000, self.z_nappe),  
            pb.Vector3d(-2000, 2000, self.z_nappe)   
        )
        
        # 2. Le domaine global (la "boîte" de la simulation)
        domaine = pb.SDF_PlantBox(4000, 4000, 4000)
        
        # 3. L'espace navigable est le domaine MOINS la limite de la roche
        espace_navigable = pb.SDF_Difference(domaine, limite_roche)
        
        return espace_navigable

    def get_model_name(self):
        return "BrusselsSoil"

