from .BrusselsSoil import BrusselsSoil
from .NorwegianSoil import NorwegianSoil
from .BrusselsSoil_Dumux import creer_sol_bruxellois_dumux
from .UrbanSoil import creer_sol_urbain_dumux
from .UrbanSoilBrussels_Dumux import creer_sol_urbain_bruxellois_dumux
from .UrbanSoilBrussels import UrbanSoilBrussels
from .urban_sol_module import UrbanSoilDumuxCoupled

def exporter_sol_vtk(sol, etape, dossier_destination):
    """
    Exporte l'état hydrologique de la grille du sol au format VTK (Structured Points).
    Permet la visualisation de la dynamique de l'eau en 3D dans le logiciel ParaView.
    
    Args:
        sol (SolDynamique): L'instance du sol contenant la matrice self.grid.
        etape (int): Le numéro de l'itération temporelle (utilisé pour nommer le fichier de sortie).
    """
    filename = f"{dossier_destination}/Picea_Abies_Humidite_{etape:03d}.vtk"
    with open(filename, "w") as f:
        f.write("# vtk DataFile Version 3.0\n")
        f.write(f"Teneur en Eau Etape {etape}\n")
        f.write("ASCII\n")
        f.write("DATASET STRUCTURED_POINTS\n")
        f.write(f"DIMENSIONS {sol.res_x} {sol.res_y} {sol.res_z}\n")
        f.write(f"ORIGIN {sol.xmin} {sol.ymin} {sol.zmin}\n")

        dx = (sol.xmax - sol.xmin) / (sol.res_x - 1)
        dy = (sol.ymax - sol.ymin) / (sol.res_y - 1)
        dz = (sol.zmax - sol.zmin) / (sol.res_z - 1)
        f.write(f"SPACING {dx} {dy} {dz}\n")
        
        f.write(f"POINT_DATA {sol.res_x * sol.res_y * sol.res_z}\n")
        f.write("SCALARS theta float 1\n")
        f.write("LOOKUP_TABLE default\n")
        
        for k in range(sol.res_z):
            for j in range(sol.res_y):
                for i in range(sol.res_x):
                    f.write(f"{sol.grid[i, j, k]:.4f}\n")