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

def main4():
    body = naca_airfoil("0012", 20)
    body = Pitching(body, 10, 2 * np.pi)
    Uinfty = (1,0)
    dt = 0.1
    flow = Simulation(body, Uinfty, dt, BoundVortexPanels, FreeVortexParticles)
    for i in range(num_steps):
        flow.advance(dt)
        # get any output

# try out staticmethods with derived classes
class A(object):
    def foo(self):
        print("My name is %s" % self.name())

    @staticmethod
    def name():
        return "Alice"

class B(A):
    @staticmethod
    def name():
        return "Bob"
