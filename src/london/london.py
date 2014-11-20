from connector import GimicConnector
import numpy as np
import current

class London(GimicConnector):
    def __init__(self, xdens):
        self.basis, self.dmat = current.get_bas_dmat(xdens)

    def jvector(self, r):
        rho, grad_rho, jp = current.vdens(r, self.basis, self.dmat)
        return jp

    def rho(self, r):
        rho, grad_rho, jp = current.vdens(r, self.basis, self.dmat)
        return rho

    def grad_rho(self, r):
        rho, grad_rho, jp = current.vdens(r, self.basis, self.dmat)
        return grad_rho

    def calc(self, r):
        return vdens(r, self.basis, self.dmat)

