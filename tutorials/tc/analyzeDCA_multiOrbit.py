import numpy as np
from numpy import *
import h5py
import sys
import os
import matplotlib.pyplot as plt
import pandas as pd
from plotnine import *

# sys.path.append("/Users/t7m/CODES/DCAPP_maierta/DCA/tools/python_scripts/")
# import symmetrize_Nc4x4


class Analyze:

    def __init__(self, fileG4, channel="PARTICLE_HOLE_MAGNETIC", readG4=True, model="square", bin=0):
        self.fileG4 = fileG4
        self.readG4 = readG4
        self.channel = channel
        self.model = model
        self.bin = bin
        self.readData()
        self.build_bStar()
        self.build_bStar(cluster=True)

        # self.calcLatticeSusceptibility()

    def calcLatticeSusceptibility(self):
        self.susC = self.calcChiFromG4(self.G4)
        self.calcChi0Cluster()
        self.buildFullChi0Lattice(nkfine=32)
        self.calcG4Lattice()
        self.susL = self.calcChiFromG4(self.G4L)

    def readData(self):
        f = h5py.File(self.fileG4, 'r')
        self.cluster = stack(list(f["step_"+str(self.bin)+"/domains/CLUSTER/REAL_SPACE/super-basis"]), axis=0)
        # print("Cluster vectors:",self.cluster)

        # Reciprocal lattice vectors
        self.b0 = f["step_"+str(self.bin)+"/domains/CLUSTER/MOMENTUM_SPACE/super-basis"][0]
        self.b1 = f["step_"+str(self.bin)+"/domains/CLUSTER/MOMENTUM_SPACE/super-basis"][1]

        # Basis vectors for cluster
        self.bc0 = f["step_"+str(self.bin)+"/domains/CLUSTER/MOMENTUM_SPACE/basis"][0]
        self.bc1 = f["step_"+str(self.bin)+"/domains/CLUSTER/MOMENTUM_SPACE/basis"][1]

        self.iwm = array(f["step_"+str(self.bin)+'/parameters']['four-point']['frequency-transfer'])[0]  # transferred frequency in units of 2*pi*temp
        self.qchannel = array(f["step_"+str(self.bin)+'/parameters']['four-point']['momentum-transfer'])
        self.invT = array(f["step_"+str(self.bin)+'/parameters']['physics']['beta'])[0]
        self.temp = 1.0/self.invT
        if self.model in ["square", "triangular", "Kagome"]:
            self.U = array(f["step_"+str(self.bin)+'/parameters']['single-band-Hubbard-model']['U'])[0]
            self.tp = array(f["step_"+str(self.bin)+'/parameters']['single-band-Hubbard-model']['t-prime'])[0]
            self.t = np.array(f["step_"+str(self.bin)+'/parameters']['single-band-Hubbard-model']['t'])[0]
        elif self.model in ["La3Ni2O7"]:
            self.U = array(f["step_"+str(self.bin)+'/parameters']['La3Ni2O7-bilayer-model']['U'])[0]
            self.V = array(f["step_"+str(self.bin)+'/parameters']['La3Ni2O7-bilayer-model']['V'])[0]
            self.J = array(f["step_"+str(self.bin)+'/parameters']['La3Ni2O7-bilayer-model']['J'])[0]
            self.t11 = np.array(f["step_"+str(self.bin)+'/parameters']['La3Ni2O7-bilayer-model']['t11'])[0]
            self.t22 = np.array(f["step_"+str(self.bin)+'/parameters']['La3Ni2O7-bilayer-model']['t22'])[0]
            self.t12 = np.array(f["step_"+str(self.bin)+'/parameters']['La3Ni2O7-bilayer-model']['t12'])[0]
            self.tperp11 = np.array(f["step_"+str(self.bin)+'/parameters']['La3Ni2O7-bilayer-model']['t-perp_11'])[0]
            self.tperp22 = np.array(f["step_"+str(self.bin)+'/parameters']['La3Ni2O7-bilayer-model']['t-perp_22'])[0]
            self.Delta = np.array(f["step_"+str(self.bin)+'/parameters']['La3Ni2O7-bilayer-model']['Delta'])[0]
        elif self.model in ["two-orbital"]:
            self.U = array(f["step_"+str(self.bin)+'/parameters']['two-orbital-model']['U'])[0]
            self.V = array(f["step_"+str(self.bin)+'/parameters']['two-orbital-model']['V'])[0]
            self.J = array(f["step_"+str(self.bin)+'/parameters']['two-orbital-model']['J'])[0]
            self.t11 = np.array(f["step_"+str(self.bin)+'/parameters']['two-orbital-model']['t11'])[0]
            self.t22 = np.array(f["step_"+str(self.bin)+'/parameters']['two-orbital-model']['t22'])[0]
            self.t12 = np.array(f["step_"+str(self.bin)+'/parameters']['two-orbital-model']['t12'])[0]
        elif self.model in ["bilayer-Hubbard-model"]:
            self.U = array(f["step_"+str(self.bin)+'/parameters']['bilayer-Hubbard-model']['U'])[0]
            self.t= np.array(f["step_"+str(self.bin)+'/parameters']['bilayer-Hubbard-model']['t'])[0]
            self.tprime= np.array(f["step_"+str(self.bin)+'/parameters']['bilayer-Hubbard-model']['t-prime'])[0]
            self.tperp= np.array(f["step_"+str(self.bin)+'/parameters']['bilayer-Hubbard-model']['t-perp'])[0]
        self.mu = np.array(f["step_"+str(self.bin)+'/DCA-loop-functions/chemical-potential'])[0]
        self.fill = array(f["step_"+str(self.bin)+'/parameters']['physics']['density'])[0]
        self.dens = array(f["step_"+str(self.bin)+'/DCA-loop-functions']['density'])
        self.nk = array(f["step_"+str(self.bin)+'/DCA-loop-functions']['n_k'])

        self.Green = array(f["step_"+str(self.bin)+'/functions/cluster_greens_function_G_k_w'])[:, :, 0, :, 0, :]
        self.Nc = self.Green.shape[1]
        self.NwG = self.Green.shape[0]
        self.iwG0 = int(self.NwG/2)
        self.sigma = array(f["step_"+str(self.bin)+'/functions/Self_Energy'])[:, :, 0, :, 0, :]
        self.wn = np.array(f["step_"+str(self.bin)+'/domains']['frequency-domain']['elements'])
        self.wnSet = np.array(f["step_"+str(self.bin)+'/domains']['vertex-frequency-domain (COMPACT)']['elements'])
        self.Kvecs = stack(list(f["step_"+str(self.bin)+'/domains/CLUSTER/MOMENTUM_SPACE/elements']), axis=0)
        self.NwTP = 2*np.array(f["step_"+str(self.bin)+'/parameters']['domains']['imaginary-frequency']['four-point-fermionic-frequencies'])[0]
        self.qmcSign = list(f["step_"+str(self.bin)+'/DCA-loop-functions/sign'])
        self.wnG = np.array(f["step_"+str(self.bin)+'/domains/frequency-domain/elements'])

        if self.readG4:
            self.G4 = array(f["step_"+str(self.bin)+'/functions/G4_'+self.channel])[..., :, :, :, :]  # last 4 indices are the orbitals
            print("G4 shape:", self.G4.shape)
            if self.channel == "PARTICLE_HOLE_MAGNETIC":
                self.G4 = swapaxes(self.G4, 4, 6)
                self.G4 = swapaxes(self.G4, 5, 7)
            elif "PARTICLE_PARTICLE" in self.channel:
                # self.G4 = swapaxes(self.G4, 2, 4)
                # self.G4 = swapaxes(self.G4, 3, 5)
                # self.G4 = swapaxes(self.G4, 6, 9)
                # self.G4 = swapaxes(self.G4, 7, 8)
                # self.G4 = swapaxes(self.G4, 4, 6)
                # self.G4 = swapaxes(self.G4, 5, 7)

                if self.model != "La3Ni2O7":
                    self.G4 = self.G4.swapaxes(4, 6).swapaxes(5, 7)
                else:
                    self.G4 = self.G4.swapaxes(4, 9).swapaxes(5, 8).swapaxes(6, 9).swapaxes(7, 8)

            else:
                print("Channel not implemented; bailing out")
                exit()

            # G4 shape now: G4(wm, Q, wn, K, l1, l2, wn', K', l3, l4)
            self.Qindices = array(f["step_"+str(self.bin)+'/domains']['Momentum exchange domain.']['element_indices_'])
            self.wnSet = np.array(f["step_"+str(self.bin)+'/domains']['vertex-frequency-domain (COMPACT)']['elements'])
            self.Nwm = self.G4.shape[0]
            self.Nq = self.G4.shape[1]
            self.NwG4 = self.G4.shape[2]
            self.nOrb = self.G4.shape[4]
            self.nt = self.Nc*self.NwG4*self.nOrb*self.nOrb
            self.iwG40 = int(self.NwG4/2)
            self.G4 = self.G4.reshape(self.Nwm*self.Nq, self.nt, self.nt)  # Now G4 is a 3D array

        self.aVecs = np.zeros((self.nOrb, 2))
        if self.model == "Kagome":
            self.aVecs = np.array([[0, 0], [0.5, 0], [0.25, np.sin(np.pi / 3.0) / 2.0]])

        self.setupMomentumTables()

        # Remove vaccuum term for charge channel G4
        if self.readG4 & (self.channel == "PARTICLE_HOLE_CHARGE"):
            for iQ in range(self.G4.shape[0]):
                (iwm, iq) = np.unravel_index(iQ, (self.Nwm, self.Nq))
                if (iq == self.iK0) & (iwm == 0):

                    for iK1 in range(self.G4.shape[1]):
                        (iw1, ik1, l1, l2) = np.unravel_index(iK1, (self.NwG4, self.Nc, self.nOrb, self.nOrb))
                        for iK2 in range(self.G4.shape[2]):
                            (iw2, ik2, l3, l4) = np.unravel_index(iK2, (self.NwG4, self.Nc, self.nOrb, self.nOrb))
                            iw1Green = iw1 - self.iwG40 + self.iwG0
                            iw2Green = iw2 - self.iwG40 + self.iwG0
                            self.G4[iQ, iK1, iK2] -= 2.0 * self.Green[iw1Green, ik1, l1, l2] * self.Green[iw2Green, ik2, l4, l3]

        # Add phase factors when K+Q is outside 1. BZ
        # For PHM channel, we have
        #
        #   l1, K+Q          l3, K'+Q
        #    ---<-----------<----
        #           |   |                     x exp(-i*(K+Q-[K+Q]_1BZ)*r_l1) x exp(i*(K'+Q-[K'+Q]_1BZ)*r_l3)
        #           |   |
        #    --->----------->----
        #    l2, K           l4, K'
        #
        self.addPhaseFactors()

        f.close()

    def addPhaseFactors(self):
        if self.readG4 & (self.channel == "PARTICLE_HOLE_MAGNETIC"):
            print("Adding phase factors to G4")
            G4 = self.unfoldG4(self.G4)
            for iq in range(self.Nq):
                Q = self.Kvecs[self.Qindices[iq]]
                for ik1 in range(self.Nc):
                    K1 = self.Kvecs[ik1]
                    iK1pQ = self.iKSum[ik1, self.Qindices[iq]]  # index of K1+Q mapped to 1.BZ
                    for ik2 in range(self.Nc):
                        K2 = self.Kvecs[ik2]
                        iK2pQ = self.iKSum[ik2, self.Qindices[iq]]  # index of K2+Q mapped to 1.BZ
                        for l1 in range(self.nOrb):
                            for l3 in range(self.nOrb):
                                pf1 = np.exp(-1j*np.dot(K1+Q-self.Kvecs[iK1pQ, :], self.aVecs[l1, :]))
                                pf2 = np.exp(1j*np.dot(K2+Q-self.Kvecs[iK2pQ, :], self.aVecs[l3, :]))
                                G4[:, iq, :, ik1, l1, :, :, ik2, l3, :] *= pf1 * pf2

    def dispersion(self, kx, ky):
        ek = np.zeros((self.nOrb, self.nOrb), dtype='complex')
        if self.model == "square":  # 1-orbital model
            ek[0, 0] = -2.*self.t*(cos(kx)+cos(ky)) - 4.0*self.tp*cos(kx)*cos(ky)
        elif self.model == "triangular":  # 1-orbital model
            ek[0, 0] = -2. * self.t * cos(kx) - 4. * self.t * cos(sqrt(3.) * ky / 2.) * cos(kx / 2.)
        elif self.model == "Kagome":  # 3-orbital model
            ek[0, 1] = -2.0 * self.t * np.cos(0.5 * kx)
            ek[0, 2] = -2.0 * self.t * np.cos(0.25 * kx + 0.25 * np.sqrt(3.0) * ky)
            ek[1, 2] = -2.0 * self.t * np.cos(0.25 * kx - 0.25 * np.sqrt(3.0) * ky)

            ek[1, 0] = -2.0 * self.t * np.cos(0.5 * kx)
            ek[2, 0] = -2.0 * self.t * np.cos(0.25 * kx + 0.25 * np.sqrt(3.0) * ky)
            ek[2, 1] = -2.0 * self.t * np.cos(0.25 * kx - 0.25 * np.sqrt(3.0) * ky)
        elif self.model == "La3Ni2O7":  # two-orbital bilayer model for La3Ni2O7
            val11 = -2. * self.t11 * (cos(kx) + cos(ky))
            val22 = -2. * self.t22 * (cos(kx) + cos(ky))
            val12 = 2. * self.t12 * (cos(kx) - cos(ky))

            ek[0, 0] = val11 + self.Delta
            ek[1, 1] = val22
            ek[2, 2] = val11 + self.Delta
            ek[3, 3] = val22

            ek[0, 1] = val12
            ek[1, 0] = val12
            ek[2, 3] = val12
            ek[3, 2] = val12

            ek[0, 2] = -self.tperp11
            ek[2, 0] = -self.tperp11
            ek[1, 3] = -self.tperp22
            ek[3, 1] = -self.tperp22

        elif self.model == "two-orbital":  # generic two-orbital model
            val11 = -2. * self.t11 * (cos(kx) + cos(ky))
            val22 = -2. * self.t22 * (cos(kx) + cos(ky))
            val12 = -4. * self.t12 * sin(kx) * sin(ky)

            ek[0, 0] = val11
            ek[1, 1] = val22

            ek[0, 1] = val12
            ek[1, 0] = val12

        elif self.model == "bilayer-Hubbard-model":  # generic two-orbital model
            val = -2. * self.t * (cos(kx) + cos(ky)) - 4. * self.tprime * cos(kx)*cos(ky)

            ek[0, 0] = val
            ek[1, 1] = val

            ek[0, 1] = -self.tperp
            ek[1, 0] = -self.tperp

        else:
            print("Model not implemented")
            sys.exit(0)
        return ek

    def calcChi0Cluster(self):
        print("Now calculating chi0 on cluster")
        NwG = self.NwG

        self.chic0 = np.zeros_like(self.G4)
        # Load frequency and K domain
        wnSet = self.wnSet
        # Kset = self.shift21BZ(self.Kvecs)
        Kset = self.Kvecs

        for iQ in range(self.G4.shape[0]):
            (iwm, iq) = np.unravel_index(iQ, (self.Nwm, self.Nq))

            for iwn, wn in enumerate(wnSet):  # reduced tp frequencies !
                for iK, K in enumerate(Kset):  # cluster K
                    if "PARTICLE_HOLE" in self.channel:
                        iKG2 = self.iKSum[iK, self.Qindices[iq]]  # k+Q
                        iwG = int(iwn - self.iwG40 + self.iwG0)
                        iwPlusiwm = self.iwnPlusiwm(iwG, iwm, maxIndex=NwG)  # iwn+iwm
                        G1 = self.Green[iwG, iK, :, :]
                        G2 = self.Green[iwPlusiwm, iKG2, :, :]
                        cc = np.outer(G2.reshape(self.nOrb*self.nOrb), G1.reshape(self.nOrb*self.nOrb))

                        for ind1 in range(self.nOrb*self.nOrb):
                            (l1, l2) = np.unravel_index(ind1, (self.nOrb, self.nOrb))
                            iK1 = np.ravel_multi_index((iwn, iK, l1, l2), dims=(self.NwG4, self.Nc, self.nOrb, self.nOrb))
                            for ind2 in range(self.nOrb*self.nOrb):
                                (l3, l4) = np.unravel_index(ind2, (self.nOrb, self.nOrb))
                                iK2 = np.ravel_multi_index((iwn, iK, l3, l4), dims=(self.NwG4, self.Nc, self.nOrb, self.nOrb))
                                ccind1 = np.ravel_multi_index((l3, l1), dims=(self.nOrb, self.nOrb))
                                ccind2 = np.ravel_multi_index((l2, l4), dims=(self.nOrb, self.nOrb))
                                # Add phase factor exp((K+Q-K[iK2])*r_l1-r_l2)
                                Q = self.Kvecs[self.Qindices[iq]]
                                pf = np.exp(1j*np.dot(K+Q-Kset[iKG2], self.aVecs[l3, :]-self.aVecs[l1, :]))
                                self.chic0[iQ, iK1, iK2] = -cc[ccind1, ccind2] * pf

                            # Note that the disconnected vaccum like term has already been subtracted for the charge channel
                            # for the particle-hole-magnetic channel, it is in principle there for q=0 and wm=0, but the
                            # sum over the spins with a factor of sigma*sigma' cancels it for spin symmetric models
                            # For Rashba-Hubbard and Moire-Hubbard models, this term is finite since the spin degeneracy is lifted
                            # We therefore treat those models with a different analysis script.

                    elif "PARTICLE_PARTICLE" in self.channel:
                        iKG2 = self.iKSum[self.iKDiff[self.iK0, iK], self.Qindices[iq]]  # -k+Q
                        iwG = int(iwn - self.iwG40 + self.iwG0)
                        miwPlusiwm = self.iwnPlusiwm(NwG-iwG-1, iwm, maxIndex=NwG)  # -iwn+iwm
                        G1 = self.Green[iwG, iK, :, :]
                        G2 = self.Green[miwPlusiwm, iKG2, :, :]
                        cc = np.outer(G1.reshape(self.nOrb*self.nOrb), G2.reshape(self.nOrb*self.nOrb))

                        for ind1 in range(self.nOrb*self.nOrb):
                            (l1, l2) = np.unravel_index(ind1, (self.nOrb, self.nOrb))
                            iK1 = np.ravel_multi_index((iwn, iK, l1, l2), dims=(self.NwG4, self.Nc, self.nOrb, self.nOrb))
                            for ind2 in range(self.nOrb*self.nOrb):
                                (l3, l4) = np.unravel_index(ind2, (self.nOrb, self.nOrb))
                                iK2 = np.ravel_multi_index((iwn, iK, l3, l4), dims=(self.NwG4, self.Nc, self.nOrb, self.nOrb))

                                ccind1 = np.ravel_multi_index((l1, l3), dims=(self.nOrb, self.nOrb))
                                ccind2 = np.ravel_multi_index((l2, l4), dims=(self.nOrb, self.nOrb))
                                # ccind1 = np.ravel_multi_index((l2, l4), dims=(self.nOrb, self.nOrb))
                                # ccind2 = np.ravel_multi_index((l3, l1), dims=(self.nOrb, self.nOrb))

                                # Add phase factor exp((K+Q-K[iK2])*r_l1-r_l2)
                                Q = self.Kvecs[self.Qindices[iq]]
                                pf = np.exp(1j*np.dot(K+Q-Kset[iKG2], self.aVecs[l3, :]-self.aVecs[l1, :]))
                                self.chic0[iQ, iK1, iK2] = cc[ccind1, ccind2] * pf
                    else:
                        print("Channel not implemented! Bailing out")
                        exit()

    # def calcGammaIrr(self):
    #     # Calculate the irr. vertex GammaIrr
    #     Nq=self.Nq; Nwm=self.Nwm; Nc=self.Nc; NwG4=self.NwG4; NwG=self.NwG; nt = self.nt
    #     self.GammaM = np.zeros((Nwm, Nq, nt, nt), dtype='complex')
    #     for iwm in range(Nwm):
    #         for iq in range(Nq):
    #             G4 = self.G4[iwm,iq,...].reshape(nt,nt)
    #             chic0 = np.diag(self.chic0[iwm,iq,...].reshape(nt))
    #             G4 = linalg.inv(G4)
    #             chic0 = linalg.inv(chic0)
    #             self.GammaM[iwm, iq, ...] = chic0 - G4
    #             self.GammaM[iwm, iq, ...] *= float(Nc)*self.invT
    #     self.Gamma = self.GammaM.reshape(Nwm, Nq, NwG4, Nc, NwG4, Nc)


    def calcChiFromG4(self, G4, form_factor=None):
        #  G4 data layout: G4(Q, K, K')
        #  Without form factor: Sum over K, K'
        if form_factor is None:
            form_factor = self.swave
        G4 = self.unfoldG4(G4)

        susK = np.zeros((self.Nq, self.Nwm, self.Nc, self.nOrb, self.nOrb, self.Nc, self.nOrb, self.nOrb))
        sus = np.zeros((self.Nq, self.Nwm, self.nOrb, self.nOrb, self.nOrb, self.nOrb))
        for iwm in range(self.Nwm):
            for iQ in range(self.Nq):
                susK[iQ, iwm, ...] = np.sum(G4[iwm, iQ, :, :, :, :, :, :, :, :].real, axis=(0, 4))

        for iK1 in range(self.Nc):
            gk1 = form_factor(self.Kvecs[iK1, 0], self.Kvecs[iK1, 1])
            for iK2 in range(self.Nc):
                gk2 = form_factor(self.Kvecs[iK2, 0], self.Kvecs[iK2, 1])
                sus[:, :, ...] += gk1 * susK[:, :, iK1, :,:, iK2, :,:] * gk2

        sus /= self.invT * self.Nc
        # print("Susceptibility from G4: ", sus)
        return sus[0,0, ...], form_factor

    def calcChiFromG4_Test(self, G4):
        #  G4 data layout: G4(Q, K, K')
        #  Without form factor: Sum over K, K'
        G4 = self.unfoldG4(G4)[0, 0, ...]

        susK = np.zeros((self.Nc, self.nOrb, self.nOrb, self.Nc, self.nOrb, self.nOrb))
        sus = 0.0
        susD = 0.0
        susDM = np.zeros((self.nOrb, self.nOrb))
        susK = np.sum(G4[:, :, :, :, :, :, :, :].real, axis=(0, 4))

        gks = self.swave(self.Kvecs[:, 0], self.Kvecs[:, 1])
        gkxs = self.xswave(self.Kvecs[:, 0], self.Kvecs[:, 1])
        gkspm = self.spmwave(self.Kvecs[:, 0], self.Kvecs[:, 1])
        gkd = self.dwave(self.Kvecs[:, 0], self.Kvecs[:, 1])

        for iK1 in range(self.Nc):
            for iK2 in range(self.Nc):
                for l1 in range(self.nOrb):
                    for l2 in range(self.nOrb):
                        for l3 in range(self.nOrb):
                            for l4 in range(self.nOrb):
                                if (l1//2==l1/2) == (l2//2==l2/2): # both even or both odd --> same orbital
                                    gk1 = gks[iK1]+gkxs[iK1]+gkspm[iK1]
                                else:
                                    gk1 = gkd[iK1]
                                if (l3//2==l3/2) == (l4//2==l4/2): # both even or both odd --> same orbital
                                    gk2 = gks[iK2]+gkxs[iK2]+gkspm[iK2]
                                else:
                                    gk2 = gkd[iK2]
                                sus += gk1 * susK[iK1, l1, l2, iK2, l3, l4] * gk2
                                if (l3==l1) & (l4==l2):
                                    susD += gk1 * susK[iK1, l1, l2, iK2, l1, l2] * gk2
                                    susDM[l1, l2] += gk1 * susK[iK1, l1, l2, iK2, l1, l2] * gk2
        return susDM

        sus /= self.invT * self.Nc
        print("Susceptibility from G4 Test: ", sus, susD)
        # return sus[0,0, ...], form_factor

    def swave(self, kx, ky):
        return np.ones_like(self.Kvecs[:,0])

    def xswave(self, kx, ky):
        return np.cos(kx) + np.cos(ky)

    def spmwave(self, kx, ky):
        return np.cos(kx) * np.cos(ky)

    def dwave(self, kx, ky):
        return np.cos(kx) - np.cos(ky)

    def buildFullChi0Lattice(self, nkfine=128):
        Nq = self.Nq
        Nwm = self.Nwm
        Nc = self.Nc
        NwG4 = self.NwG4
        nOrb = self.nOrb
        nt = NwG4*Nc*nOrb*nOrb
        self.chi0 = np.zeros_like(self.G4)

        kPatch = self.build_kGrid(cluster=True, nkfine=nkfine)

        for iQ in range(self.G4.shape[0]):
            self.chi0[iQ, :, :] = self.buildChi0Lattice(iQ, kPatch)

    def buildChi0Lattice(self, iQ, kPatch):
        (iwm, iq) = np.unravel_index(iQ, (self.Nwm, self.Nq))
        Qvec = self.Kvecs[self.Qindices[iq], :]
        print("Now calculating chi0 on lattice for Q=", Qvec," and iwm=", iwm)

        NwG = self.NwG

        # Load frequency and K domain
        wnSet = self.wnSet
        # Kset = self.shift21BZ(self.Kvecs)
        Kset = self.Kvecs

        # Now coarse-grain G*G to build chi0(K) = Nc/N sum_k Gc(K+k')Gc(-K-k')
        nOrb = self.nOrb
        chi0 = np.zeros_like(self.G4[iQ, :, :])

        # Pre-calculate energies
        ek1 = np.zeros((self.Nc, kPatch.shape[0], self.nOrb, self.nOrb), dtype='complex')
        ek2 = np.zeros((self.Nc, kPatch.shape[0], self.nOrb, self.nOrb), dtype='complex')
        for iK, K in enumerate(Kset):  # cluster K
            for ik, k in enumerate(kPatch):
                kx = K[0] + k[0]
                ky = K[1] + k[1]
                if "PARTICLE_HOLE" in self.channel:
                    ek1[iK, ik, ...] = self.dispersion(kx, ky)
                    ek2[iK, ik, ...] = self.dispersion(kx + Qvec[0], ky + Qvec[1])  # e(k+Q)
                elif "PARTICLE_PARTICLE" in self.channel:
                    ek1[iK, ik, ...] = self.dispersion(kx, ky)
                    ek2[iK, ik, ...] = self.dispersion(-kx + Qvec[0], -ky + Qvec[1])  # e(-k+Q)

        self.cG1 = np.zeros((self.NwG4, self.Nc, self.nOrb, self.nOrb), dtype='complex')
        self.cG2 = np.zeros((self.NwG4, self.Nc, self.nOrb, self.nOrb), dtype='complex')
        for iwn, wn in enumerate(wnSet):  # reduced tp frequencies !
            for iK, K in enumerate(Kset):  # cluster K
                cc = np.zeros((self.nOrb*self.nOrb, self.nOrb*self.nOrb), dtype='complex')
                cG1 = np.zeros((self.nOrb, self.nOrb), dtype='complex')
                cG2 = np.zeros((self.nOrb, self.nOrb), dtype='complex')

                if "PARTICLE_HOLE" in self.channel:
                    iK2 = self.iKSum[iK, self.Qindices[iq]]  # k+Q
                    for ik, k in enumerate(kPatch):
                        e1 = ek1[iK, ik, ...]
                        e2 = ek2[iK, ik, ...]
                        iwG = int(iwn - self.iwG40 + self.iwG0)
                        # minusiw = self.minusiwn(iwG, maxIndex=NwG) # -iwn
                        iwPlusiwm = self.iwnPlusiwm(iwG, iwm, maxIndex=NwG)  # iwn+iwm
                        sigmaK1 = self.sigma[iwG, iK, :, :]
                        sigmaK2 = self.sigma[iwPlusiwm, iK2, :, :]
                        G1 = np.linalg.inv((1j*wn + self.mu)*np.identity(nOrb) - e1 - sigmaK1)
                        G2 = np.linalg.inv((1j*self.wnG[iwPlusiwm] + self.mu)*np.identity(nOrb) - e2 - sigmaK2)
                        cc += -np.outer(G2.reshape(self.nOrb*self.nOrb), G1.reshape(self.nOrb*self.nOrb))
                        cG1 += G1
                        cG2 += G2

                    for ind1 in range(self.nOrb*self.nOrb):
                        (l1, l2) = np.unravel_index(ind1, (self.nOrb, self.nOrb))
                        iK1 = np.ravel_multi_index((iwn, iK, l1, l2), dims=(self.NwG4, self.Nc, self.nOrb, self.nOrb))
                        self.cG1[iwn, iK, :, :] = cG1
                        self.cG2[iwn, iK, :, :] = cG2
                        for ind2 in range(self.nOrb*self.nOrb):
                            (l3, l4) = np.unravel_index(ind2, (self.nOrb, self.nOrb))
                            iK2 = np.ravel_multi_index((iwn, iK, l3, l4), dims=(self.NwG4, self.Nc, self.nOrb, self.nOrb))
                            ccind1 = np.ravel_multi_index((l3, l1), dims=(self.nOrb, self.nOrb))
                            ccind2 = np.ravel_multi_index((l2, l4), dims=(self.nOrb, self.nOrb))
                            chi0[iK1, iK2] = cc[ccind1, ccind2]

                elif "PARTICLE_PARTICLE" in self.channel:
                    iK2 = self.iKSum[self.iKDiff[self.iK0, iK], self.Qindices[iq]]  # -k+Q
                    for ik, k in enumerate(kPatch):
                        e1 = ek1[iK, ik, ...]
                        e2 = ek2[iK, ik, ...]
                        iwG = int(iwn - self.iwG40 + self.iwG0)
                        miwPlusiwm = self.iwnPlusiwm(NwG-iwG-1, iwm, maxIndex=NwG)  # iwn+iwm
                        sigmaK1 = self.sigma[iwG, iK, :, :]
                        sigmaK2 = self.sigma[miwPlusiwm, iK2, :, :]
                        G1 = np.linalg.inv((1j*wn + self.mu)*np.identity(nOrb) - e1 - sigmaK1)
                        G2 = np.linalg.inv((1j*self.wnG[miwPlusiwm] + self.mu)*np.identity(nOrb) - e2 - sigmaK2)
                        cc += np.outer(G1.reshape(self.nOrb*self.nOrb), G2.reshape(self.nOrb*self.nOrb))

                    for ind1 in range(self.nOrb*self.nOrb):
                        (l1, l2) = np.unravel_index(ind1, (self.nOrb, self.nOrb))
                        iK1 = np.ravel_multi_index((iwn, iK, l1, l2), dims=(self.NwG4, self.Nc, self.nOrb, self.nOrb))
                        self.cG1[iwn, iK, :, :] = cG1
                        self.cG2[iwn, iK, :, :] = cG2
                        for ind2 in range(self.nOrb*self.nOrb):
                            (l3, l4) = np.unravel_index(ind2, (self.nOrb, self.nOrb))
                            iK2 = np.ravel_multi_index((iwn, iK, l3, l4), dims=(self.NwG4, self.Nc, self.nOrb, self.nOrb))
                            ccind1 = np.ravel_multi_index((l1, l3), dims=(self.nOrb, self.nOrb))
                            ccind2 = np.ravel_multi_index((l2, l4), dims=(self.nOrb, self.nOrb))
                            chi0[iK1, iK2] = cc[ccind1, ccind2]

        chi0 /= kPatch.shape[0]
        self.cG1 /= kPatch.shape[0]
        self.cG2 /= kPatch.shape[0]

        return chi0

    def calcG4Lattice(self):
        nt = self.nt
        self.G4L = np.zeros_like(self.G4)
        self.GammaIrr = np.zeros_like(self.G4)
        for iQ in range(self.G4.shape[0]):
            G4 = self.G4[iQ, ...]
            chic0 = self.chic0[iQ, ...]
            chi0 = self.chi0[iQ, ...]
            G4Inv = np.linalg.inv(G4)
            chic0Inv = np.linalg.inv(chic0)
            chi0Inv = np.linalg.inv(chi0)
            self.G4L[iQ, ...] = np.linalg.inv(chi0Inv - chic0Inv + G4Inv)
            self.GammaIrr[iQ, ...] = (chic0Inv - G4Inv) * self.invT * self.Nc

    def diagonalizePPKernel(self, iQ=0):
        self.pm = np.dot(self.GammaIrr[iQ, ...], self.chi0[iQ, ...]) / (self.invT * self.Nc)
        w, v = np.linalg.eig(self.pm)
        wt = abs(w - 1)
        ilead = np.argsort(wt)
        self.lambdas = w[ilead]
        self.evecs = v[:, ilead]

    def diagonalizeSymmPPKernel(self, iQ=0):
        chiD = linalg.eig(self.chi0[iQ, ...])

# Helper functions

    def unfoldG4(self, G4):
        return G4.reshape(self.Nwm, self.Nq, self.NwG4, self.Nc, self.nOrb, self.nOrb, self.NwG4, self.Nc, self.nOrb, self.nOrb)

    def unfold_K_Index(self, index):
        return np.unravel_index(index, (self.NwG4, self.Nc, self.nOrb, self.nOrb))

    def setupMomentumTables(self):
        # build tables for K+K' and K-K'
        self.iK0 = self.K_2_iK(0.0, 0.0)
        self.iKDiff = zeros((self.Nc, self.Nc), dtype='int')
        self.iKSum = zeros((self.Nc, self.Nc), dtype='int')
        Nc = self.Nc
        for iK1 in range(Nc):
            Kx1 = self.Kvecs[iK1, 0]
            Ky1 = self.Kvecs[iK1, 1]
            for iK2 in range(0, Nc):
                Kx2 = self.Kvecs[iK2, 0]
                Ky2 = self.Kvecs[iK2, 1]
                iKS = self.K_2_iK(Kx1+Kx2, Ky1+Ky2)
                iKD = self.K_2_iK(Kx1-Kx2, Ky1-Ky2)
                self.iKDiff[iK1, iK2] = iKD
                self.iKSum[iK1, iK2] = iKS

    def K_2_iK(self, Kx, Ky):
        delta = 1.0e-5
        K = np.array([Kx, Ky])
        M = np.column_stack((self.b0, self.b1))
        MTM = np.dot(M.T, M)
        MTMinv = np.linalg.inv(MTM)

        for iK in range(0, self.Nc):
            K_minus_Ktarget = K - self.Kvecs[iK]
            Y = np.dot(np.dot(MTMinv, M.T), K_minus_Ktarget)
            # print("iK, Y: ", iK, Y)
            if (abs(Y[0] - round(Y[0])) < delta) & (abs(Y[1] - round(Y[1])) < delta):
                return iK
        print("No Kvec found!!!", Kx, Ky)

    def iwnPlusiwm(self, iwn, iwm, maxIndex): # iwn + iwm (wn fermionic, wm bosonic)
        return int(min([max([iwn + iwm, 0]), maxIndex-1]))

    def minusiwn(self, iwn, maxIndex):  # -iwn (wn fermionic)
        return maxIndex-iwn-1

    def iw1Minusiw2(self, iw1, iw2, maxIndex):  # both iw1 and iw2 are fermionic
        iww = min(max(iw1 - iw2, 0), maxIndex)
        return iww

    def build_bStar(self, cluster=False):
        #  Build 8 rec. lattice vectors that surround K0=(0,0)
        if cluster is False:
            b0 = self.b0  # use rec. lattice vectors
            b1 = self.b1
        else:
            b0 = self.bc0  # use basis vectors for cluster K-points
            b1 = self.bc1
        b = np.zeros((8, 2))
        b[0] = b0
        b[1] = b1
        b[2] = -b0
        b[3] = -b1
        b[4] = b0 + b1
        b[5] = b0 - b1
        b[6] = -b0 + b1
        b[7] = -b0 - b1
        if cluster is False:
            self.bStar = b
        else:
            self.bcStar = b

    def shift21BZ(self, Kin):
        b = self.bStar
        #  Now check if K-point is closest to K0=(0,0) and, if not, subtract from it the b-vector that it is closest to
        K = Kin.copy()
        for iK, Kvec in enumerate(K):
            dist0 = np.linalg.norm(Kvec)
            distClosest = dist0
            ibClosest = -1  # corresponds to K=(0,0)
            for ib in range(8):
                distb = np.linalg.norm(Kvec - b[ib])
                if distb < distClosest:
                    distClosest = distb
                    ibClosest = ib
            if ibClosest > -1:
                Kvec -= b[ibClosest]
        return K

    def build_kGrid(self, nkfine=128, cluster=False):
        #  build fine k mesh for 1. BZ if cluster = False, for DCA patch if cluster = True
        if cluster is False:
            b = self.bStar
            b0 = self.b0
            b1 = self.b1
        else:
            b = self.bcStar
            b0 = self.bc0
            b1 = self.bc1

        k = []
        for ik1 in range(nkfine):
            for ik2 in range(nkfine):
                k.append(ik1/float(nkfine)*b0 + ik2/float(nkfine)*b1)

        #  Now check if k-point is closest to K0=(0,0) and, if not, subtract the b-vector
        for iK, Kvec in enumerate(k):
            dist0 = np.linalg.norm(Kvec)
            distClosest = dist0
            ibClosest = -1  # corresponds to K=(0,0)
            for ib in range(8):
                distb = np.linalg.norm(Kvec - b[ib])
                if distb < distClosest:
                    distClosest = distb
                    ibClosest = ib
            if ibClosest > -1:
                Kvec -= b[ibClosest]
        k = np.stack(list(k))
        return k

#  Plotting functions

    def plotKvecs(self, shiftTo1BZ=True, plotBZ=True, plotPatch=True):
        if shiftTo1BZ is True:
            Kvecs = self.shift21BZ(self.Kvecs)
        else:
            Kvecs = self.Kvecs
        if plotBZ is True:
            k = self.build_kGrid()
        if plotPatch is True:
            kP = self.build_kGrid(cluster=True)
        fig, ax = plt.subplots(figsize=(6, 6))
        if plotBZ is True:
            ax.scatter(k[:, 0], k[:, 1], alpha=0.1, c="lightgrey")
        if plotPatch is True:
            ax.scatter(kP[:, 0], kP[:, 1], alpha=0.1, c="grey")

        ax.scatter(Kvecs[:, 0], Kvecs[:, 1])
        ax.annotate("", xy=(self.b0[0], self.b0[1]), xytext=(0, 0), arrowprops=dict(arrowstyle="->"))
        ax.annotate("", xy=(self.b1[0], self.b1[1]), xytext=(0, 0), arrowprops=dict(arrowstyle="->"))

        for iK in range(self.Kvecs.shape[0]):
            ax.text(x=Kvecs[iK, 0], y=Kvecs[iK, 1], s=str(iK))
        ax.set_aspect('equal')

    def plotEvecs(self, evecs, id=0):

        fig, ax = plt.subplots(nrows=self.nOrb, ncols=self.nOrb, figsize=(4*self.nOrb, 4*self.nOrb), sharey=True, sharex=True)
        for i1 in range(self.nOrb):
            for i2 in range(self.nOrb):
                for iK in range(self.Nc):
                    ax[i1, i2].plot(self.wnSet, evecs[:, id].reshape(self.NwG4, self.Nc, self.nOrb, self.nOrb)[:, iK, i1, i2].real, label="real; K="+str(iK))
                    ax[i1, i2].plot(self.wnSet, evecs[:, id].reshape(self.NwG4, self.Nc, self.nOrb, self.nOrb)[:, iK, i1, i2].imag, '--', label="imag; K="+str(iK))
                    ax[i1, i2].set_title(str(i1)+str(i2))
        [ax[i1, -1].set(xlabel=r"$\omega_n$") for i1 in range(self.nOrb)]

        # ax[1, 0].legend(bbox_to_anchor=(1.2, 1.0), loc='right')
        # plt.subplots_adjust(left=0.1, bottom=0.1, right=0.88, top=0.9, wspace=0.5, hspace=0.2)
        fig.suptitle(r"$\lambda=$"+str(self.lambdas[id]))

    def plotEvecs3(self, id, imag=False):
        if imag:
            ev = self.evecs[:, id].imag.reshape(self.NwG4, self.Nc, self.nOrb, self.nOrb)
        else:
            ev = self.evecs[:, id].real.reshape(self.NwG4, self.Nc, self.nOrb, self.nOrb)

        dd = pd.DataFrame([[self.wnSet[j], iK,  str(s1)+str(s2), ev[j, iK, s1, s2]] for s1 in range(
            self.nOrb) for s2 in range(self.nOrb) for iK in range(self.Nc) for j in range(self.NwG4)])
        dd.columns = [r"$\omega_n$", "K", "orbital", r"$\phi'$"]
        # theme_set(theme_minimal)
        self.evplot = (ggplot(dd, aes(x=r"$\omega_n$", y=r"$\phi'$", color="factor(orbital)"))
                       + geom_line()
                       + facet_wrap('K')
                       + labs(color="orbital")
                       + theme_minimal(base_size=16, base_family="Arial")
                       + theme(panel_grid_major=element_line(color="Darkgrey"))
                       + theme(figure_size=(10, 10))
                       + theme(plot_background=element_rect(fill='white', color='white'))
                       # + theme(strip_background=element_rect(fill="Darkgrey", size=1.4, alpha=.95),)
                       )

        print(self.evplot)

if __name__ == '__main__':
    leading_eigenvalues = []
    leading_eigenvecs = []
    for filename in ['T=0.1/dca_tp.hdf5','T=0.09/dca_tp.hdf5','T=0.08/dca_tp.hdf5','T=0.07/dca_tp.hdf5']:
        print('\n=== Processing', filename, '===')
        T_analysis = Analyze(filename,channel='PARTICLE_PARTICLE_UP_DOWN')
        T_analysis.calcChi0Cluster()
        T_analysis.buildFullChi0Lattice()
        T_analysis.calcG4Lattice()
        T_analysis.diagonalizePPKernel()
        print('Leading eigenvalues (sorted by distance from 1):')
        print(T_analysis.lambdas[:5])
        print('Leading eigenvalue (real part):', T_analysis.lambdas[0].real)
        leading_eigenvalues.append(T_analysis.lambdas.real)
        leading_eigenvecs.append(T_analysis.evecs[:,0].reshape(T_analysis.NwG4, T_analysis.Nc).real)
