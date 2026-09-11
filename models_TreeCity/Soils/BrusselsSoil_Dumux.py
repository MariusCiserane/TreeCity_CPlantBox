# Attention : selon ton installation, cela peut être 'rosi_richards' ou 'dumux_rosi'
from rosi.richards import RichardsWrapper 
from rosi.rosi_richards import RichardsSP

def creer_sol_bruxellois_dumux(box_min=[-300.0, -300.0, -300.0], box_max=[300.0, 300.0, 0.0], cell_number=[20, 20, 20]):
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
    vg_quaternaire = [0.078, 0.430, 0.036, 1.56, 24.96]
    vg_tielt       = [0.089, 0.390, 0.021, 1.48, 31.40]
    vg_aalbeke     = [0.068, 0.380, 0.008, 1.09, 4.80]
    vg_moen        = [0.100, 0.380, 0.027, 1.23, 2.88]
    vg_stmaur      = [0.068, 0.380, 0.008, 1.09, 4.80]
    vg_grandglise  = [0.045, 0.430, 0.145, 2.68, 712.80]
    vg_lincent     = [0.068, 0.380, 0.008, 1.09, 4.80]
    vg_cretace     = [0.034, 0.460, 0.016, 1.37, 6.00]
    vg_paleozoique = [0.010, 0.050, 0.001, 1.10, 0.01]

    lst_vg = [vg_quaternaire, vg_tielt, vg_aalbeke, vg_moen, vg_stmaur, vg_grandglise, vg_lincent, vg_cretace, vg_paleozoique]

    # 4. ATTRIBUTION DES COUCHES 
    # setLayersZ définit la limite *inférieure* de chaque couche.
    # Ici : Couche 1 va jusqu'à -100 cm. Couche 2 va jusqu'à -300 cm.
    s.setLayersZ(number=[i for i in range(len(lst_vg))], z=[-100.0, -2100.0, -2500.0, -6000, -9000, -9800, -11800, -13800, -15000])
    
    # On passe la liste complète des paramètres (dans le même ordre)
    s.setVGParameters(lst_vg)

    # 5. CONDITIONS AUX LIMITES
    s.setTopBC("noFlux")
    s.setBotBC("noFlux") 

    # 6. CONDITION INITIALE (Équilibre Hydrostatique)
    # setICZ prend la profondeur de la nappe phréatique (ex: -150 cm).
    s.setHomogeneousIC(-150.0, True)

    s.setParameter("Soil.SourceSlope", "100")
    # Finalisation
    s.initializeProblem()

    # 7. PRESSION CRITIQUE
    # On la garde, c'est une excellente sécurité numérique pour Richards !
    s.setCriticalPressure(-15000.0)
    
    return s