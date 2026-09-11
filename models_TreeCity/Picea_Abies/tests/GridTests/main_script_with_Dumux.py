# gestion des chemins
import os
from pathlib import Path
import sys

sys.path.append(str(Path(__file__).resolve().parent.parent.parent.parent))

import plantbox as pb
from Soils import UrbanSoilDumuxCoupled
import numpy as np
from rosi.richards import RichardsWrapper
from rosi.rosi_richards import RichardsSP
from plantbox.functional.PlantHydraulicModel import HydraulicModel_Doussan
from plantbox.functional.PlantHydraulicParameters import PlantHydraulicParameters

def fonction_elongation_bengough(pr_value):
    """
    Renvoie une grille (SoilLookUp) qui applique la fonction de ralentissement de l'élongation selon Bengough et al. (2011).
    """
                
    if pr_value <= 0.5:
        return 1.0 
    elif pr_value >= 2.5:
        return 0.0 
    else:
        return max(0.0, 1.0 - ((pr_value - 0.5) / 2.0))

def soif_selon_age(age_jours):
    """
    Calcule le besoin journalier (mL) en fonction de l'âge de l'arbre.
    Exemple linéaire : 
    - À 0 jour : 1 mL (graine/germe)
    - À 365 jours : 50 mL (jeune sapin)
    """
    besoin_min = 1.0
    besoin_a_1_an = 50.0
    
    # Règle de 3 mathématique (croissance linéaire)
    # Tu peux changer cette formule pour une courbe exponentielle si tu as des données biblio !
    besoin_actuel = besoin_min + (age_jours / 365.0) * (besoin_a_1_an - besoin_min)
    return besoin_actuel

def run_simulation():
    profil_id = 2
    script_dir = Path(__file__).resolve().parent
    excel_path = f"{script_dir.parent.parent.parent.parent}/modelparameter_TreeCity/soils/SoilData_en.xlsx"  # À ajuster selon votre vrai nom de fichier
    xml_path = f"{script_dir.parent.parent.parent.parent}/modelparameter_TreeCity/structural/Picea_Abies_v2.xml"# À ajuster selon votre fichier arbre
    hydraulic_parameters_path = f"{script_dir.parent.parent.parent.parent}/modelparameter_TreeCity/functional/plant_hydraulics/Picea_Abies" # À ajuster selon votre fichier hydraulique
    export_dir = script_dir / "results" / "With_Dumux" / f"UrbanSoilTuranProfil{profil_id}_epaississement"
    export_dir.mkdir(parents=True, exist_ok=True)
    
    # --- Paramètres de la Grille ---
    nx, ny, nz = 20, 20, 20
    min_b = [-50., -50., -150.]
    max_b = [50., 50., 0.]
    
    
    # --- 1. INITIALISATION DE DUMUX ---
    s = RichardsWrapper(RichardsSP())
    s.initialize()
    s.createGrid(min_b, max_b, [nx, ny, nz])

    # --- 3. SCRIPTE D'ENVIRONNEMENT (IC / BC) ---
    # Condition Initiale : Profil hydrostatique stable centré sur -100 cm de charge
    s.setHomogeneousIC(-100.0, True)
    
    # Conditions aux limites : Surface étanche, drainage libre au fond
    s.setTopBC("constantFlux", 0.0)
    s.setBotBC("noFlux")
    
    # Initialisation finale du problème C++
    #s.initializeProblem()
    
    # --- 2. CONFIGURATION GÉOMÉTRIQUE VIA NOTRE CLASSE ---
    # Lit l'excel, calcule les PTF et crée les domaines physiques dans DuMux
    sol = UrbanSoilDumuxCoupled(excel_path, profil_id, nx, ny, nz, min_b, max_b, s_dumux=s)
    
    
    # Micro-solve CRUCIAL pour forcer l'allocation de getWaterContent() au jour 0
    print("Génération de la carte d'humidité initiale...")
    s.solve(1e-6)
    
    # Première synchronisation mécanique
    sol.update_hydromechanics()
    
    # --- 4. INITIALISATION DE L'ARBRE (CPlantBox) ---
    plant = pb.MappedPlant()
    plant.enableExtraNode()
    plant.readParameters(str(xml_path))

    print("Création de la grille de perception biologique de l'arbre...")
    # On crée une grille qui appartient à la logique de la plante
    grille_biologique = pb.EquidistantGrid3D(
        sol.min_b[0], sol.max_b[0], sol.nx,  # Axe X
        sol.min_b[1], sol.max_b[1], sol.ny,  # Axe Y
        sol.min_b[2], sol.max_b[2], sol.nz   # Axe Z
    )

    # Assignation de la grille au multiplicateur C++
    for subType in range(1, 7):
        try:
            rrpm = plant.getOrganRandomParameter(pb.OrganTypes.root, subType)
            rrpm.f_se = grille_biologique
        except Exception:
            pass
        
    plant.initialize()

    def picker(x, y, z):
        return s.pick([x, y, z])
    plant.setSoilGrid(picker)
    
    # Modèle hydraulique de Doussan pour l'aspiration de l'eau

    p_flux = PlantHydraulicParameters() # Charge les conductivités kr et kx par défaut
    # ATTENTION : Il faut vérifier les valeurs dans le fichier de paramètres
    # état actuel : produit par IA, non vérifié
    #ce read_parameter rajoute de lui-même l'extension .json
    p_flux.read_parameters(hydraulic_parameters_path)
    modele_hydraulique = HydraulicModel_Doussan(plant, p_flux)
    modele_hydraulique.wilting_point = -15000.0 # Valeur par défaut, à ajuster selon l'espèce et le sol

    # --- CORRECTION : INITIALISATION BIOLOGIQUE AVANT LA CROISSANCE ---
    # 1. On calcule les facteurs de Bengough au jour 0
    facteurs_bengough_init = np.clip(1.0 - ((sol.pr_array - 0.5) / 2.0), 0.0, 1.0)
    
    # 2. On injecte ces facteurs initiaux dans la grille
    for k in range(sol.nz):
        for j in range(sol.ny):
            for i in range(sol.nx):
                grille_biologique.setData(i, j, k, facteurs_bengough_init[i, j, k])

    plant_age = 10.0 # On laisse pousser l'arbre 10 jours "à vide"
    plant.simulate(plant_age, True)
    
    # --- 5. BOUCLE TEMPORELLE SIMULATION COUPLÉE ---
    sim_time = 3000.0     # Durée totale de la simulation (jours)
    dt = 1.0            # Pas de temps (1 jour)
    export_interval = 10 # Exporter les données tous les jours
    transpiration_cible = 50.0 # cm3/jour
    taux_epaississement = 0.002 # Vitesse de croissance radiale (cm/jour)
    
    print("\nLancement de la boucle dynamique TreeCity...")
    for t in range(int(sim_time / dt)):
        print(f"\nJour {t+1}/{int(sim_time)} :")
        
        # --- ÉTAPE A : Flux Hydrauliques & Assèchement du sol ---
        h_s = s.getSolutionHead() # Récupère l'état hydrique du sol

        real_age = plant_age + t*dt
        t_pot = soif_selon_age(real_age)
        
        # Résolution du modèle Doussan (Potentiels xylémiques h_x et flux d'extraction)
        #h_x = modele_hydraulique.solve(sim_time = t, t_act=-transpiration_cible, rsx=h_s, cells=True)
        #fluxes = modele_hydraulique.soil_fluxes(sim_time=t, rx=h_x, rsx=h_s)

        h_x = modele_hydraulique.solve(sim_time=real_age, t_act = -t_pot, rsx = h_s, cells = True)
        fluxes = modele_hydraulique.soil_fluxes(sim_time = real_age, rx = h_x, rsx = h_s)

        eau_extraite = sum(fluxes.values())
        print(f"   -> L'arbre a réellement extrait : {eau_extraite:.2f} cm3 d'eau aujourd'hui")
        
        # Injection des flux de succion de l'arbre dans DuMux
        s.setSource(fluxes)
        
        # Résolution de l'équation de Richards : le sol change d'humidité !
        s.solve(dt)
        
        # --- ÉTAPE B : Traduction mécanique de la nouvelle humidité ---
        # Met à jour automatiquement self.pr_array via self.s.getWaterContent()
        sol.update_hydromechanics()

        # L'arbre traduit la physique en biologie (Vitesse d'exécution : NumPy !)
        # La formule de Bengough est : 1.0 - ((PR - 0.5) / 2.0) plafonnée entre 0 et 1.
        # NumPy calcule cela pour les 64000 voxels d'un seul coup !
        facteurs_bengough = np.clip(1.0 - ((sol.pr_array - 0.5) / 2.0), 0.0, 1.0)
    
        # On injecte ces facteurs dans la grille C++ de la plante
        # (C'est la seule petite boucle requise pour passer de Python au C++)
        for k in range(sol.nz):
            for j in range(sol.ny):
                for i in range(sol.nx):
                    grille_biologique.setData(i, j, k, facteurs_bengough[i, j, k])
            
        # Avancement de la croissance de l'arbre
        plant.simulate(dt, False)

        # --- Epaississement progressif des racines ---
        ana = pb.SegmentAnalyser(plant)
        
        creation_times = np.array(ana.getParameter("creationTime"))
        subtypes = np.array(ana.getParameter("subType"))
        
        ages_reels_du_segment = plant_age - creation_times
        
        # On prépare le tableau des rayons
        rayons_evolutifs = np.zeros(len(creation_times))
        
        # 3. Calcul de l'épaississement progressif (à ajuster en fonction de la réalité biologique)
        for i in range(len(creation_times)):
            if subtypes[i] == 1: 
                 # Pivot : Rayon de base + (âge du segment * taux)
                rayons_evolutifs[i] = 0.20 + (ages_reels_du_segment[i] * taux_epaississement)
                                
            elif subtypes[i] == 4 or subtypes[i] == 6:
                # Charpentières : S'épaississent un peu moins vite
                rayons_evolutifs[i] = 0.15 + (ages_reels_du_segment[i] * (taux_epaississement / 2))
                                
            elif subtypes[i] == 5 or subtypes[i] == 3:
                rayons_evolutifs[i] = 0.05 + (ages_reels_du_segment[i] * (taux_epaississement / 4))
            else:
                # Racines fines : Ne s'épaississent pas (gardent leur taille de naissance)
                rayons_evolutifs[i] = 0.05 # On peut ajuster le taux pour les racines fines si nécessaire
        # On écrase le rayon par défaut avec notre cône parfait
        ana.addData("radius", rayons_evolutifs.tolist())
        
        # note bonus : il devrait être possible de ne garder que les racines qui nous intéressent 
        # avec un ana.filter("subType", minSubType, maxSubType) (les segments dont le subType n'appartient pas à [minSubType, maxSubType] seront supprimés)
        # si un seul type de racine n'est gardé, il y a aussi ana.filter("subType", val)
        # filter est un filtre qui peut être utilisé sur tous les paramètres , pas que sur subType, il suffit de mettre le nom correspondant
        # le paramètre order peut aussi être intéressant pour ne garder que les racines les plus importantes 
        
        plant_age += dt
        
        # --- ÉTAPE D : Sauvegarde des résultats ---
        if t % export_interval == 0:
            print(f"  -> Export ParaView Jour {t+1}")
            sol.export_paraview(export_dir / f"UrbanSoil_Resistance_{t//export_interval}.vtk")
            #plant.write(f"{export_dir}/Tree_Architecture_{t//export_interval}.vtp", True)
            filename = f"{export_dir}/Tree_Architecture_{t//export_interval}.vtp"
            ana.write(filename, ["age", "creationTime", "id", "order", "organType", "radius", "subType"])
            
    print("\nSimulation couplée TreeCity terminée avec succès !")

if __name__ == "__main__":
    run_simulation()