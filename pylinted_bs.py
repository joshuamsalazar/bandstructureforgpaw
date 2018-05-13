#!/usr/bin/python
"""Interface to BandElement and BandStructure classes for
Computational Physics I Final Project
(c) April 2018"""

from matplotlib import pyplot as plt
from gpaw import GPAW

calc = GPAW('graphene_ground_state_pw_PBEsol_G-K-M.gpw')

class BandElementInterface(object):
    """Interface for an element of a band"""
    def __init__(self, nband, kpt, evalue):    # band_structure, calc, Extra arguments
        """Constructor for the BandElement
        nband (int) band number
        kpt (int) k-point number
        band_structure l(ist) list of BandElement objects [nbands, nkpts]
        """
        self.nbands = calc.get_number_of_bands()
        self.nband = nband
        self.kpt = kpt
        self.evalue = evalue    #Eigenvalue of the Band Element
        self.overlaps = {}    # Overlaps
        self.best_overlap = 0    #Best current overlap for iterating
        self.next = None    # Pointer to next element in band
        self.prev = None    # Pointer to previous element in band

    def __repr__(self):
        return "band "+str(self.nband)+" k-point "+str(self.kpt)

    def set_overlaps(self):
        """Set the array of overlaps <psi_nk|psi_n'k+1>
        overlaps (list) array of overlaps between this band
        and the next k-point of length nbands"""
        raise NotImplementedError

    def match(self, minoverlap=0.15):
        """Match this BandElement with the next k-point
        If the maximum overlap is less than minoverlap
        the BandElement is left unmatched"""
        #self.kmatch = overlaps[c_best][0]

    def print_band(self):
        """Print the band to the passed file handle"""
        band = [[], []]
        while self.next is not None:
            #kpt = self.kpt
            band[0].append(self.kpt)
            band[1].append(self.evalue)
            self = self.next
        return band

class BandStructureInterface(object):
    """Interface for a band structure class"""
    def __init__(self, calc, refenergy=None):
        """Generate Band Structure for Calculator
        calc (GPAW calc object)
        refenergy (real) Fermi level or energy reference
        If None the Fermi level is taken from calc"""
        self.calc = calc
        self.nbands = calc.get_number_of_bands()
        self.nkpts = len(calc.get_ibz_k_points())
        if refenergy is None:
            refenergy = calc.get_fermi_level()
        self.refenergy = refenergy
        self.bands = []
        self.initialize_band_structure()
        self.initialize_overlaps()
        self.match_bands()

    def initialize_band_structure(self):
        """Initialize the band structure matrix (n,kpt) with BandElements"""
        self.bs_matrix = [[BandElementInterface(n, k, 0) for k in range(self.nkpts)]
                          for n in range(self.nbands)]
        for n in range(self.nbands):
            for k in range(self.nkpts):
                self.bs_matrix[n][k].evalue = self.calc.get_eigenvalues(k)[n]

    def initialize_overlaps(self):
        """Initialize the overlap list for each band_nk """
        for n in range(self.nbands):
            for k in range(self.nkpts - 1):
                for n_prim in range(self.nbands):
                    print("Currently overlapping n: " + str(n) + ". k: " +str(k) +". n': " +
                          str(n_prim))
                    self.bs_matrix[n][k].overlaps[n_prim] = abs(
                        (self.calc.get_pseudo_wave_function(n, k).conj() *
                         self.calc.get_pseudo_wave_function(n_prim, k+1)).sum() /
                        (self.calc.get_pseudo_wave_function(n, k).conj() *
                         self.calc.get_pseudo_wave_function(n, k)).sum())
                self.bs_matrix[n][k].overlaps = sorted(
                    self.bs_matrix[n][k].overlaps.items(), key=lambda x: -x[1])
                        #Cloning the dictionary into a overlap-ordered list [(nband, overlap),...}

    def match_bands(self):
        """Trace bands through the k-space
        calling the match method for each BandElement"""
        done = False
        while not done:
            for k in range(self.nkpts - 1):
                for n in range(self.nbands):
                    current_element = self.bs_matrix[n][k]
                        #Variable containing the current band element in the FOR loop process 
                    best_overlap_index = current_element.best_overlap
                            #Index number of the best overlap on the list
                    next_match = self.bs_matrix[current_element.overlaps[
                        best_overlap_index][0]][k+1]
                            #Next band element candidate to be matched
                    n_next = next_match.nband
                        # next_match.prev is the pointer to the currently matched item in conflict
                    print("Currently matching n: " + str(n) + ". k: " +str(k))
                    if next_match.prev is None:   #If the best overlap candidate is unmatched
                        current_element.next = next_match    #Then stablish new default connections
                        next_match.prev = current_element
                        print("Normally matched n: " + str(n) + ". k: " +str(k) + "with n': " +
                              str(n_next) + " " + str(k+1))

                    #Check whether the new matching item has conflict (same level of best overlap)
                    #One of them should be lower, then clear the lower
                    else:
                        if current_element.overlaps[
                                best_overlap_index][1] < next_match.prev.overlaps[
                                    self.bs_matrix[n_next][k].best_overlap][1] \
                            and current_element.overlaps[
                                best_overlap_index][0] == next_match.prev.overlaps[
                                    self.bs_matrix[n_next][k].best_overlap][0]:
                            print("Conflict n: " + str(n) + ". k: " +str(k))
                            best_overlap_index += 1
                            if k == self.nkpts - 2 and n == self.nbands - 1:
                                done = True
                            print("Ended loop")
                        else:
                            print("Unmatching n: " + str(n) + ". k: " +str(k) +"with . n: ")
                            next_match.prev = None
                            current_element.next = next_match
                            next_match.prev = current_element.next
                            if k == self.nkpts - 2 and n == self.nbands - 1:
                                done = True
                            print("Ended loop")

    def get_overlap(self, kpt1, kpt2):
        """Calculate the overlap between wavefunctions at
        kpt1, and kpt2 <psi_n,kpt1|psi_n',kpt2>
        returns a nbands by nbands array"""
        raise NotImplementedError

    def get_bands(self):
        """Print all bands to filename by calling the print function
        each BandElement at the first k-point"""
        for i in range(self.nbands):
            self.bands.append(self.bs_matrix[i][0].print_band())

    def print_bands(self, filename="tmp.dat"):
        fig = plt.figure()
        plt.title("Band Structure")
        for i in range(self.nbands-1):
            #print(self.bands[i][0][i])
            print(i)
            plt.plot(self.bands[i][0], self.bands[i][1])
        plt.xlabel('k-points')
        plt.ylabel('Eigenenergies [eV]')
        plt.grid()
        plt.legend()
        plt.show()
 #   def get_pw_overlap(calc, k1=0, k2=1):
 #       """Calculate PW overlaps for kpts k1 and k2"""
 #       calc.set_positions(calc.get_atoms())
 #       kpt1 = calc.wfs.kpt_u[k1]
 #       kpt2 = calc.wfs.kpt_u[k2]
 #       nbands = calc.get_number_of_bands()
 #       dtype = calc.wfs.dtype
 #       pw_overlap = zeros((nbands, nbands), dtype=dtype)
 #       pw_overlap = abs(dot(kpt1.psit_nG[:],
 #                            kpt2.psit_nG[:].transpose().conj()))
 #       # Normalization
 #       norm = abs(dot(kpt1.psit_nG[:],
 #                      kpt1.psit_nG[:].transpose().conj())).diagonal()
 #       pw_overlap /= norm
 #       return pw_overlap
class MainMenu():
    def show_options(self):
        print("#################################################################")
        print("#################### Band Structure Plotter #####################")
        print("#################################################################")
        print("\n\n\n\n\n")
        chosen_option = False
        while not chosen_option:
            print("\n>>>>\t Which method should I use for calculating your band structure?")
            options = ["\ta. Finite difference(Real Space)\n",
                       "\tb. Plane Wave Expansion\n",
                       "\tc. Linear Combination of Atomic Orbitals (LCAO)\n",
                       "\td. Projector Augmented Wave\n\n >>>"]
            option = raw_input(''.join([options[i] for i in range(len(options))])).lower()
            if option == "a":
                band_struc = BandStructureInterface(calc)
                band_struc.print_bands()
                chosen_option = True
            elif option == "b":
                print("No.")
            elif option == "c":
                print("Chao.")
            else:
                print("Please chose a valid option")
MENU = MainMenu()
MENU.show_options()
