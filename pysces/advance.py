__all__ = ['time_advance']

def time_advance(body, wake, Uinfty, dt):
    # todo: how does dt come in to updating strengths of body panels, and in
    # particular wake panel to be shed?

    # 1) define the wake panel: need Uinf and dt to determine length
    # 2) to determine rhs for update_strengths, need Uinf for relative vel
    # 3) wake update needs Uinf and dt (relative vel)
    wake.advect(body, Uinfty, dt)
    body.update_strengths_unsteady(wake, Uinfty, dt)
    x_shed, gam_shed = body.get_newly_shed()
    wake.add_vortex(x_shed, gam_shed)
