# Attention : selon ton installation, cela peut être 'rosi_richards' ou 'dumux_rosi'
from rosi.richards import RichardsWrapper 
from rosi.rosi_richards import RichardsSP

def creer_sol_urbain_dumux(box_min=[-300.0, -300.0, -300.0], box_max=[300.0, 300.0, -10.0], cell_number=[20, 20, 20]):
    """
    Initialise le solveur Dumux-ROSI (Équation de Richards) pour le sol Bruxellois.
    Utilise strictement les fonctions setLayersZ, setVGParameters et setICZ.
    """
    print("Initialisation de Dumux-ROSI (Sol Bruxellois)...")
    
    # 1. Instance du solveur
    s = RichardsWrapper(RichardsSP())
    s.initialize()

    # 2. Maillage (Grid) : Un cube de 6x6 mètres, sur 3 mètres de profondeur
    s.createGrid(box_min, box_max, cell_number) 

    # 3. Paramètres de Van Genuchten [theta_r, theta_s, alpha (1/cm), n, Ks (cm/j)]
    vg_base_granulaire = [0.02, 0.3, 0.13, 3.3, 864]
    vg_remblai_t_comp  = [0.089, 0.3, 0.0055, 1.25, 0.95]
    vg_remblai_comp    = [0.09, 0.36, 0.02, 1.38, 18]
    vg_ss_sol_nat      = [0.08, 0.380, 0.38, 1.53, 40]


    lst_vg = [vg_base_granulaire, vg_remblai_t_comp, vg_remblai_comp, vg_ss_sol_nat]

    # 4. ATTRIBUTION DES COUCHES 
    # setLayersZ définit la limite *inférieure* de chaque couche.
    # Ici : Couche 1 va jusqu'à -100 cm. Couche 2 va jusqu'à -300 cm.
    s.setLayersZ(number=[i for i in range(len(lst_vg))], z=[-35.0, -100, -200, box_min[2]])
    
    # On passe la liste complète des paramètres (dans le même ordre)
    s.setVGParameters(lst_vg)

    # 5. CONDITIONS AUX LIMITES
    s.setTopBC("noFlux")
    s.setBotBC("noFlux") 

    s.setHomogeneousIC(-10.0, True)

    s.setParameter("Soil.SourceSlope", "100")
    # Finalisation
    s.initializeProblem()

    # 7. PRESSION CRITIQUE
    # On la garde, c'est une excellente sécurité numérique pour Richards !
    s.setCriticalPressure(-15000.0)
    
    return s