__all__ = ['time_advance']

def time_advance(body, wake, Uinfty, dt):
    # todo: how does dt come in to updating strengths of body panels, and in
    # particular wake panel to be shed?

    # 1) define the wake panel: need Uinf and dt to determine length
    # 2) to determine rhs for update_strengths, need Uinf for relative vel
    # 3) wake update needs Uinf and dt (relative vel)
    body.update_strengths(wake, Uinfty, dt) # might need dt to determine length of TE panel
    wake.update(body, Uinfty, dt)
    shed_panels = body.get_wake_panels()
    wake.add_panels(shed_panels)
