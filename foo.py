def main_loop():
    numsteps = 10
    dt = 0.1

    for i in range(numsteps):
        time += dt
        bodies.update_kinematics(time)

def main():
    # might want to let a main loop handle everything
    bodies = make_bodies_and_specify_kinematics()
    dt = 0.1
    num_steps = 50
    logger = make_output_logger()

    main_loop(bodies, dt, num_steps, logger)

def main2():
    # or write the loop manually
    for i in range(numsteps):
        time += dt
        bodies.update_kinematics(time)
        wake.rollup(bodies, time)

def main3():
    airfoil = Body(parms)
    body_panels = VortexPanels(airfoil)  # includes wake panel to be shed
    # or other types of panels, e.g.:
    # body = SourceDoubletPanels(airfoil)
    # body = DoubletPanels(airfoil)
    wake_panels = VortexPanels()   # empty at first
    for i in range(numsteps):
        time = i * dt
        airfoil.update(time)
        body_panels.update_positions(airfoil)
        #rhs = body_panels.rhs(wake_panels)
        #Amat = body_panels.influence_matrix()
        #body_panels.solve(Amat, rhs)
        body_panels.update_strengths(wake_panels)
        wake_panels.update()
        wake_panels.add_panels(shed_panel)
