import plantbox as pb

plant = pb.Plant()
path = "../../modelparameter_TreeCity/structural/Quercus_Robur"
plant.readParameters(path + ".xml")

# Initialize
plant.initialize()

# Simulate
sim_time = 1500 
dt = 15
n_steps = round(sim_time / dt)
for i in range(0, n_steps):
    plant.simulate(dt)
    plant.write("results/Quercus_Robur_" + str(i) + ".vtp")


