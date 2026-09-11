import pandas as pd
import numpy as np
import plantbox as pb

class UrbanSoilDumuxCoupled:
    """
    Classe gérant le sol urbain hétérogène. 
    Elle lit les données Excel, calcule les PTF, configure la géométrie 3D 
    dans DuMux (pour l'eau) et dans des matrices NumPy (pour la mécanique).
    """
    def __init__(self, excel_path, profil_id, nx, ny, nz, min_b, max_b, s_dumux=None):
        self.nx, self.ny, self.nz = nx, ny, nz
        self.min_b, self.max_b = min_b, max_b
        self.s = s_dumux  
        
        print(f"Initialisation du sol urbain (Profil: {profil_id})...")
        
        # --- 1. LECTURE DES DONNÉES EXCEL (Pandas) ---
        df_profils = pd.read_excel(excel_path, sheet_name="Urban - data collected")
        self.df_vg = pd.read_excel(excel_path, sheet_name="Parameters VG per soil type")
        
        self.couches = df_profils[df_profils['Profil'] == profil_id].copy()
        if self.couches.empty:
            raise ValueError(f"Aucune donnée trouvée pour le Profil {profil_id} dans l'Excel.")
        
        # Nettoyage des données
        def clean_depth(val):
            """Convertit les profondeurs en float, gère les valeurs '>X'."""
            if isinstance(val, str) and '>' in val:
                print(f"Profondeur '>X' détectée : {val}. On prend la profondeur maximale de la boîte., ie {- self.min_b[2]} cm.")
                return - self.min_b[2]
            try:
                return float(val)
            except:
                print(f"Impossible de convertir la profondeur : {val}. On prend la profondeur maximale de la boîte., ie {- self.min_b[2]} cm.")
                return - self.min_b[2]  # Valeur par défaut si conversion impossible
            
        self.couches['Depth down\n(cm)'] = self.couches['Depth down\n(cm)'].apply(clean_depth)
        self.couches['Depth up\n(cm)'] = pd.to_numeric(self.couches['Depth up\n(cm)'], errors='coerce').fillna(0)
        self.couches = self.couches.sort_values(by='Depth down\n(cm)') # Tri par profondeur croissante
        self.couches = self.couches.drop_duplicates(subset=['Depth down\n(cm)'], keep='first') #suppression des doublons

        #On ne dépasse pas la profondeur des données dans l'Excel
        if len(self.couches) > 0:
            last_idx = self.couches.index[-1]
            if self.couches.loc[last_idx, 'Depth down\n(cm)'] < - self.min_b[2]:
                self.min_b[2] = -self.couches.loc[last_idx, 'Depth down\n(cm)']  # Ajuste la profondeur max de la boîte

        print(f"profondeur max de la boîte ajustée à {self.min_b[2]} cm selon les données Excel.")
            
        # Initialisation de la grille CPlantBox pour l'export et l'arbre
        self.grid = pb.EquidistantGrid3D(self.min_b[0], self.max_b[0], self.nx,
                                         self.min_b[1], self.max_b[1], self.ny,
                                         self.min_b[2], self.max_b[2], self.nz)
        self.pr_array = np.zeros((self.nx, self.ny, self.nz))
        
        # Matrices Python stockant les paramètres locaux pour le modèle mécanique (Dexter)
        self.grid_rho = np.zeros((self.nx, self.ny, self.nz))
        self.grid_theta_s = np.zeros((self.nx, self.ny, self.nz))
        self.grid_theta_r = np.zeros((self.nx, self.ny, self.nz))
        self.grid_alpha = np.zeros((self.nx, self.ny, self.nz))
        self.grid_n = np.zeros((self.nx, self.ny, self.nz))
        self.grid_S_index = np.zeros((self.nx, self.ny, self.nz))
        
        # Dimensions physiques d'un seul voxel (en cm)
        self.dx = (self.max_b[0] - self.min_b[0]) / self.nx
        self.dy = (self.max_b[1] - self.min_b[1]) / self.ny
        self.dz = (self.max_b[2] - self.min_b[2]) / self.nz

        # --- 2. CONSTRUCTION DE LA GÉOMÉTRIE (DuMux + Python) ---
        self.build_soil_structure()


    def build_soil_structure(self):
        """
        Calcule les paramètres via PTF et génère la géométrie 3D par écrasement successif :
        1. Couches horizontales (stratigraphie de fond)
        2. Pavés étanches (en surface)
        3. Fosse meuble (au centre, écrase les pavés et les couches)
        """
        print("--- Modélisation de la géométrie 3D ---")

        all_vg_params = [] 
        
        # =========================================================
        # ÉTAPE A : LES COUCHES HORIZONTALES (Issus de l'Excel)
        # =========================================================
        # pour étendre la couche jusqu'à la suivante, en cas de trou sans données
        self.couches['Depth_down_etendu'] = self.couches['Depth up\n(cm)'].shift(-1) 
        # Pour la toute dernière couche du tableau, on l'étire jusqu'au fond de la boîte.
        self.couches['Depth_down_etendu'] = self.couches['Depth_down_etendu'].fillna(-self.min_b[2])
        for vg_id, (idx, row) in enumerate(self.couches.iterrows()):
            #vg_id = int(idx) # ID de matériau pour DuMux
            
            
            type_sol = row.get('Soil Type', 'Sand')
            rho_b = float(row.get('bulk Density\n(g/cm3)', 1.35)) # Densité apparente (g/cm³))
            if pd.isna(rho_b): rho_b = 1.35
            
            # Récupération des paramètres de base
            vg_base = self.df_vg[self.df_vg['Soil Type'] == type_sol].iloc[0]
            # --- APPLICATION DES FORMULES PTF (Compaction Urbaine) ---
            theta_s_base = vg_base['θs (cm³/cm³)']
            theta_r_base = vg_base['θr (cm³/cm³)']
            alpha_base = vg_base['α (cm⁻¹)']
            n_base = vg_base['n']
            ks_base = vg_base.get('Ks (cm/jour)', 50.0)
            rho_b_base = 1.35 # Densité naturelle de référence

            # 1. Porosité saturée modifiée (Écrasement des pores)
            theta_s_comp = 1.0 - (rho_b / 2.65)

            # 2. Paramètre d'échelle alpha modifié (Pression d'entrée d'air)
            alpha_comp = alpha_base * (rho_b_base / rho_b)

            m = 1.0 - (1.0 / n_base)

            S_idx = -n_base * (theta_s_comp - theta_r_base) * ((1.0 + 1.0/m)**(-(1.0 + m)))

            # 3. Conductivité Ks modifiée (Équation de Kozeny-Carman)
            # Baisse drastique de la vitesse d'infiltration due à la compaction
            terme_comp = (theta_s_comp**3) * ((1.0 - theta_s_base)**2)
            terme_base = (theta_s_base**3) * ((1.0 - theta_s_comp)**2)
            ks_comp = ks_base * (terme_comp / terme_base)

            all_vg_params.append([theta_r_base, theta_s_comp, alpha_comp, n_base, ks_comp])
            
            z_top = -float(row['Depth up\n(cm)'])
            z_bottom = -float(row['Depth_down_etendu'])
                
            # Remplissage de la matrice Python pour ces profondeurs
            for k in range(self.nz):
                z_voxel = self.min_b[2] + (k + 0.5) * self.dz
                if z_bottom <= z_voxel <= z_top:
                    self.grid_theta_s[:, :, k] = theta_s_comp
                    self.grid_theta_r[:, :, k] = theta_r_base
                    self.grid_alpha[:, :, k] = alpha_comp
                    self.grid_n[:, :, k] = n_base
                    self.grid_S_index[:, :, k] = S_idx
                    self.grid_rho[:, :, k] = rho_b
                    
            print(f"  > Couche {vg_id} ({type_sol}) : Z de {z_bottom} à {z_top} cm")

        # =========================================================
        # ÉTAPE B : LES PAVÉS COMPACTÉS (Surface globale)
        # =========================================================
        paves_id = len(self.couches)
        zmin_p, zmax_p = -20.0, 0.0 # Pavés sur les 10 premiers cm de la surface
        
        # --- VOS CALCULS PTF PAVÉS ---
        qr_p, qs_p, alpha_p, n_p, ks_p, S_idx_p = 0.01, 0.15, 0.005, 1.1, 0.001, 0.01
        m_p = 1.0 - (1.0 / n_p)
        S_idx_p = -n_p * (qs_p - qr_p) * ((1.0 + 1.0/m_p)**(-(1.0 + m_p)))
        
        all_vg_params.append([qr_p, qs_p, alpha_p, n_p, ks_p])
            
        # Remplissage Python (écrase la couche du dessous)
        for k in range(self.nz):
            z_voxel = self.min_b[2] + (k + 0.5) * self.dz
            if zmin_p <= z_voxel <= zmax_p:
                self.grid_theta_s[:, :, k] = qs_p
                self.grid_theta_r[:, :, k] = qr_p
                self.grid_alpha[:, :, k] = alpha_p
                self.grid_n[:, :, k] = n_p
                self.grid_S_index[:, :, k] = S_idx_p
        print("  > Pavés appliqués en surface.")

        # =========================================================
        # ÉTAPE C : LA FOSSE MEUBLE (Centre de l'arbre)
        # =========================================================
        fosse_id = len(self.couches) + 1
        # Dimensions de la fosse centrale (ex: 60x60cm, de la surface jusqu'à -100cm)
        xmin_f, xmax_f = -30.0, 30.0
        ymin_f, ymax_f = -30.0, 30.0
        zmin_f, zmax_f = -100.0, 0.0
        
        # --- VOS CALCULS PTF FOSSE ---
        qr_f, qs_f, alpha_f, n_f, ks_f, S_idx_f = 0.05, 0.45, 0.02, 1.6, 120.0, 0.25
        m_f = 1.0 - (1.0 / n_f) 
        S_idx_f = -n_f * (qs_f - qr_f) * ((1.0 + 1.0/m_f)**(-(1.0 + m_f)))
        
        all_vg_params.append([qr_f, qs_f, alpha_f, n_f, ks_f])
            
        # Remplissage Python (écrase les pavés et les couches au centre)
        for i in range(self.nx):
            x_voxel = self.min_b[0] + (i + 0.5) * self.dx
            if xmin_f <= x_voxel <= xmax_f:
                for j in range(self.ny):
                    y_voxel = self.min_b[1] + (j + 0.5) * self.dy
                    if ymin_f <= y_voxel <= ymax_f:
                        for k in range(self.nz):
                            z_voxel = self.min_b[2] + (k + 0.5) * self.dz
                            if zmin_f <= z_voxel <= zmax_f:
                                self.grid_theta_s[i, j, k] = qs_f
                                self.grid_theta_r[i, j, k] = qr_f
                                self.grid_alpha[i, j, k] = alpha_f
                                self.grid_n[i, j, k] = n_f
                                self.grid_S_index[i, j, k] = S_idx_f
        print("  > Fosse centrale insérée (écrase le centre).")

        # =========================================================
        # INJECTION DANS DUMUX
        # =========================================================
        if self.s is not None:
            print("  > Transmission des paramètres VG à DuMux...")
            # 1. On envoie TOUS les paramètres d'un coup pour allouer le tableau C++
            # La liste plate attendue par dumux-rosi est souvent : [qr0, qs0, a0, n0, ks0, qr1, qs1, a1, n1, ks1...]
            #flat_vg_params = [item for sublist in all_vg_params for item in sublist]
            self.s.setVGParameters(all_vg_params)
            self.s.initializeProblem()
            
            print("  > Mapping spatial dans DuMux...")
            # 2. On applique les domaines géométriques
            for vg_id, (idx, row) in enumerate(self.couches.iterrows()):
                # vg_id = int(idx)
                z_top = - float(row['Depth up\n(cm)'])
                z_bottom = - float(row['Depth_down_etendu'])
                self.s.addVanGenuchtenDomain(
                    [self.min_b[0], self.min_b[1], z_bottom], 
                    [self.max_b[0], self.max_b[1], z_top], 
                    vg_id
                )
            
            # Pavés
            self.s.addVanGenuchtenDomain([self.min_b[0], self.min_b[1], zmin_p], [self.max_b[0], self.max_b[1], zmax_p], paves_id)
            # Fosse
            self.s.addVanGenuchtenDomain([xmin_f, ymin_f, zmin_f], [xmax_f, ymax_f, zmax_f], fosse_id)



    def update_hydromechanics(self):
        """ 
        Mise à jour purement mécanique. 
        Récupère l'humidité exacte de DuMux (qui gère la 3D) et applique le modèle de Dexter.
        """
        if self.s is not None:
            # Récupération directe de l'humidité volumique (theta) calculée par DuMux
            theta_flat = self.s.getWaterContent()
            theta_array = np.array(theta_flat).reshape((self.nx, self.ny, self.nz), order='F')
        else:
            # Sécurité si DuMux n'est pas utilisé
            theta_array = np.full((self.nx, self.ny, self.nz), 0.3) 
            
        a, b, c = 0.5, 0.05, 0.1 # Paramètres génériques du modèle de Dexter
        
        for i in range(self.nx):
            for j in range(self.ny):
                for k in range(self.nz):
                    theta = theta_array[i, j, k]
                    
                    # Paramètres locaux (récupérés depuis nos matrices Python)
                    theta_s = self.grid_theta_s[i, j, k]
                    theta_r = self.grid_theta_r[i, j, k]
                    alpha = self.grid_alpha[i, j, k]
                    n = self.grid_n[i, j, k]
                    S_idx = self.grid_S_index[i, j, k]
                    
                    m = 1.0 - (1.0 / n)
                    
                    # Sécurité pour éviter les divisions par zéro
                    chi = (theta - theta_r) / max(1e-5, (theta_s - theta_r))
                    chi = max(0.001, min(0.999, chi)) 

                    # Calcul de la contrainte effective (sans recalculer theta !)
                    h_eff_cm = (1.0 / alpha) * ( (chi**(-1.0/m)) - 1.0 )**(1.0/n)
                    # 3. CONVERSION D'UNITÉ : cm H2O vers MPa (1 cm H2O = 0.0000980665 MPa)
                    h_eff_MPa = h_eff_cm * 0.0000980665

                    sigma_prime = chi * h_eff_MPa
                    
                    # Résistance finale
                    PR = a + b * (1.0 / abs(S_idx)) + c * sigma_prime
                    # On empêche la formule empirique d'exploser dans les milieux ultra-secs.
                    # 10 MPa est amplement suffisant pour bloquer n'importe quelle racine d'arbre.
                    PR = min(PR, 10.0)
                    
                    self.pr_array[i, j, k] = PR
                    self.grid.setData(i, j, k, PR)
                    

    def export_paraview(self, filename="results/sol_urbain_resistance.vtk"):
        """ Exporte en VTK pour visualiser la stratigraphie et les pavés. """

        # Récupération des données DuMux
        if self.s is not None:
            h_flat = self.s.getSolutionHead()
            theta_flat = self.s.getWaterContent()
        else:
            # Sécurité au cas où DuMux n'est pas encore attaché
            h_flat = np.zeros(self.nx * self.ny * self.nz)
            theta_flat = np.zeros(self.nx * self.ny * self.nz)

        # Transformation en matrices 3D
        h_array = np.array(h_flat).reshape((self.nx, self.ny, self.nz), order='F')
        theta_array = np.array(theta_flat).reshape((self.nx, self.ny, self.nz), order='F')


        # (Laissez votre code d'export intact ici, il fonctionnait très bien)
        with open(filename, 'w') as f:
            f.write("# vtk DataFile Version 3.0\nUrban Soil\nASCII\nDATASET STRUCTURED_POINTS\n")
            f.write(f"DIMENSIONS {self.nx + 1} {self.ny + 1} {self.nz + 1}\n")
            f.write(f"ORIGIN {self.min_b[0]} {self.min_b[1]} {self.min_b[2]}\n")
            f.write(f"SPACING {self.dx} {self.dy} {self.dz}\n")
            f.write(f"CELL_DATA {self.nx * self.ny * self.nz}\n")

            # --- 1ère Variable : Résistance à la Pénétration ---
            f.write("SCALARS Penetration_Resistance float 1\n")
            f.write("LOOKUP_TABLE default\n")
            # Écriture des données (attention à l'ordre 'F' de Fortran utilisé par CPlantBox)
            for k in range(self.nz):
                for j in range(self.ny):
                    for i in range(self.nx):
                        f.write(f"{self.pr_array[i, j, k]:.4f}\n")

            # --- 2ème Variable : Pression (SolutionHead) ---
            f.write("SCALARS Pressure_Head float 1\n")
            f.write("LOOKUP_TABLE default\n")
            for k in range(self.nz):
                for j in range(self.ny):
                    for i in range(self.nx):
                        f.write(f"{h_array[i, j, k]:.4f}\n")

            # --- 3ème Variable : Teneur en eau (Water Content) ---
            f.write("SCALARS Water_Content float 1\n")
            f.write("LOOKUP_TABLE default\n")
            for k in range(self.nz):
                for j in range(self.ny):
                    for i in range(self.nx):
                        f.write(f"{theta_array[i, j, k]:.4f}\n")