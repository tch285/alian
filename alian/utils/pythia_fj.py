def psj_from_particle_with_index(particle, index):
    psj = fj.PseudoJet(particle.px(), particle.py(), particle.pz(), particle.e())
    psj.set_user_index(index)
    return psj
