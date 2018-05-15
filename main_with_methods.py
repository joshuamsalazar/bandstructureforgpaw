#!/usr/bin/python
"""Interface to BandElement and BandStructure classes for
Computational Physics I Final Project
(c) April 2018"""

from numpy import *
from matplotlib import pyplot as plt
from gpaw import GPAW
import sys

calc = GPAW(str(sys.argv[1]))

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
        self.counter = 0
        

    def initialize_band_structure(self):
        """Initialize the band structure matrix (n,kpt) with BandElements"""
        self.bs_matrix = [[BandElementInterface(n, k, 0) for k in range(self.nkpts)]
                          for n in range(self.nbands)]
        for n in range(self.nbands):
            for k in range(self.nkpts):
                self.bs_matrix[n][k].evalue = self.calc.get_eigenvalues(k)[n]

    def initialize_pw_overlaps(self):
        """Initialize the overlap list for each band_nk """
        for k in range(self.nkpts-1):
            overlaps = self.get_pw_overlap(k,k+1)
            for n in range(self.nbands):
                for n_prim in range(self.nbands):   
                    percentage = k/float(self.nkpts)*100 + 1/float(self.nkpts)/self.nbands*n*100   
                        #For LCAO, the program runs through the ns, given k (viceversa for RS)
                    print("Currently overlapping n: " + str(n) + ". k: " +str(k) +". n': " +
                          str(n_prim) + "\t" +str(round(percentage, 3)) + " percent completed.")
                    self.bs_matrix[n][k].overlaps[n_prim] = overlaps[n][n_prim]#mod
                self.bs_matrix[n][k].overlaps = sorted(self.bs_matrix[n][k].overlaps.items(), key=lambda x: -x[1])
        print("PW failed "+str(self.counter)+"times.")
                #Cloning the dictionary into a overlap-ordered list [(nband, overlap),...}

    def initialize_lcao_overlaps(self):
        """Initialize the overlap list for each band_nk """
        for k in range(self.nkpts-1):
            overlaps = self.get_lcao_overlap(k,k+1)
            for n in range(self.nbands):
                for n_prim in range(self.nbands):   
                    percentage = k/float(self.nkpts)*100 + 1/float(self.nkpts)/self.nbands*n*100   
                        #For LCAO, the program runs through the ns, given k (viceversa for RS)
                    print("Currently overlapping n: " + str(n) + ". k: " +str(k) +". n': " +
                          str(n_prim) + "\t" +str(round(percentage, 3)) + " percent completed.")
                    self.bs_matrix[n][k].overlaps[n_prim] = overlaps[n][n_prim]#mod
                self.bs_matrix[n][k].overlaps = sorted(self.bs_matrix[n][k].overlaps.items(), key=lambda x: -x[1])
                #Cloning the dictionary into a overlap-ordered list [(nband, overlap),...}
                
    def get_overlap(self, kpt1, kpt2):
        """Calculate the overlap between wavefunctions at
        kpt1, and kpt2 <psi_n,kpt1|psi_n',kpt2>
        returns a nbands by nbands array
        this version of get_overlap is for the real space"""
        overlaps = [[] for k in range(self.nbands)]
        for i in range(self.nbands):
                #Obtain the first wave function wf_nk1"
            wf_nk1 = calc.get_pseudo_wave_function(i,kpt1)
                #Normalize wf_nk1"
            N = (wf_nk1.conjugate()*wf_nk1).sum()
            for j in range(self.nbands):
                    #Obtain the second wave function wf_nk2"
                wf_nk2 = calc.get_pseudo_wave_function(j,kpt2)
                    #Get the overlap"
                overlp = abs(((wf_nk1.conjugate()*wf_nk2).sum())/N)
                    #"Add to the overlap list of the band element (kpt1,i)"
                overlaps[i].append(overlp)
        return overlaps

    def initialize_rs_overlaps(self):
        """Initialize the overlap list for each band_nk """
        for n in range(self.nbands):
            for k in range(self.nkpts - 1):
                for n_prim in range(self.nbands):
                    percentage = n/float(self.nbands)*100 + 1/float(self.nbands)/self.nkpts*k*100  
                    print("Currently overlapping n: " + str(n) + ". k: " +str(k) +". n': " +
                          str(n_prim) + "\t" +str(round(percentage, 3)) + " percent completed.")
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
            for k in range(self.nkpts - 2):
                for n in range(self.nbands):
                    current_element = self.bs_matrix[n][k]
                        #Variable containing the current band element in the FOR loop process 
                    best_overlap_index = current_element.best_overlap   #Index number of the best overlap on the list
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
                            if k == self.nkpts - 3 and n == self.nbands - 1:
                                done = True
                                print("Ended loop")
                        else:
                            print("Unmatching n: " + str(n) + ". k: " +str(k) +"with . n: ")
                            next_match.prev = None
                            current_element.next = next_match
                            next_match.prev = current_element.next
                            if k == self.nkpts - 3 and n == self.nbands - 1:
                                done = True
                                print("Ended loop")

    def get_bands(self):
        """Print all bands to filename by calling the print function
        each BandElement at the first k-point"""
        for i in range(self.nbands):
            self.bands.append(self.bs_matrix[i][0].print_band())

    def print_bands(self, filename="tmp.dat"):
        fig = plt.figure()
        plt.title("Band Structure")
        for i in range(self.nbands-1):
            plt.plot(self.bands[i][0], self.bands[i][1])
        plt.xlabel('k-points')
        plt.ylabel('Eigenenergies [eV]')
        plt.grid()
        plt.legend()
        plt.show()

    def get_lcao_overlap(self, k1, k2):
        """Calculate LCAO overlaps for kpts k1 and k2"""
        calc.set_positions(calc.get_atoms())
        kpt1 = calc.wfs.kpt_u[k1]
        kpt2 = calc.wfs.kpt_u[k2]
        nbands = calc.get_number_of_bands()
        dtype = calc.wfs.dtype
        lcao_overlap = zeros((nbands, nbands), dtype=dtype)
        lcao_overlap = abs(dot(kpt1.C_nM, dot(kpt1.S_MM,
                                              kpt2.C_nM.transpose().conj())))
        return lcao_overlap

    def get_proj_ani(self, k):
        """Obtain Projector Matrix"""
        # We can work in the IBZ since we have no momentum transfer
        nks = calc.wfs.kd.get_rank_and_index(s=0, k=k)[1]
        kpt = calc.wfs.kpt_u[nks]
        nsym = calc.wfs.kd.sym_k[k]
        proj_ani = {}
        for nat in range(len(calc.wfs.setups.id_a)):
            a_sa = calc.wfs.kd.symmetry.a_sa[nsym, nat]
            proj_ani[nat] = dot(kpt.P_ani[a_sa],
                                calc.wfs.setups[nat].R_sii[nsym])
        return proj_ani

    def get_paw_overlap(self, k1, k2):
        """Calculate PAW Corrections to the overlap qvnm matrix"""
        proj_ani = self.get_proj_ani(k1)
        # PAW corrections from other spin channel
        proj_ami = self.get_proj_ani(k2)
        nbands = calc.get_number_of_bands()
        dtype = calc.wfs.dtype
        paw_overlap_qvnm = zeros((nbands, nbands), dtype=dtype)
        setups = calc.wfs.setups
        for nat in range(len(calc.wfs.setups.id_a)):
            paw_overlap_qvnm += dot(proj_ani[nat],
                                    proj_ami[nat].transpose())
        return paw_overlap_qvnm
    
    def get_pw_overlap(self, k1, k2):
        """Calculate PW overlaps for kpts k1 and k2"""
        calc.set_positions(calc.get_atoms())
        kpt1 = calc.wfs.kpt_u[k1]
        kpt2 = calc.wfs.kpt_u[k2]
        nbands = calc.get_number_of_bands()
        dtype = calc.wfs.dtype
        pw_overlap = zeros((nbands, nbands), dtype=dtype)
        w1 = kpt1.psit_nG[:]
        w2 = kpt2.psit_nG[:].transpose().conj()
        if w1.shape[1] != w2.shape[0]:
            dim_min = min(w1.shape[1],w2.shape[0])
            if w1.shape[1] < w2.shape[0]:
                w2 = resize(w2,(dim_min,w2.shape[1]))
            if w1.shape[1] > w2.shape[0]:
                w1 = resize(w1,(w1.shape[0],dim_min))
            pw_overlap = abs(dot(w1,w2))
        N = abs(dot(w1,w1.conjugate().transpose())).diagonal()
        pw_overlap /= N
        for i in range(self.nbands):
            if amax(pw_overlap[i])<0.6:
                self.counter += 1
                return self.get_overlap(k1,k2)
        return pw_overlap

def display_notif(text):
    print("\n"*24)
    print("#"*60)
    print("\t\t" + str(text))
    print("#"*60)
    print("\n"*6)

class MainMenu():
    def display_menu(self):
        print("\n"*10)
        display_notif("Band Structure interpolator v0.2")
        print("\n"*6)
        chosen_option = False
        while not chosen_option:
            print("\n>>>>\t Which method should I use for calculating your band structure?\n\n")
            options = ["\ta. Finite difference(Real Space)\n",
                       "\tb. Plane Wave Expansion\n",
                       "\tc. Linear Combination of Atomic Orbitals (LCAO)\n",
                       "\td. Projector Augmented Wave\n",
                       "\te. Close the program\n\nEnter your option >>> "]
            option = raw_input(''.join([options[i] for i in range(len(options))])).lower()
            if option == "a":
                band_struc = BandStructureInterface(calc)
                band_struc.initialize_rs_overlaps()
                band_struc.match_bands()
                band_struc.get_bands()
                chosen_option = True
                display_notif("Finished Calculation")
            elif option == "b":
                band_struc = BandStructureInterface(calc)
                band_struc.initialize_pw_overlaps()
                band_struc.match_bands()
                band_struc.get_bands()
                chosen_option = True
                display_notif("Finished Calculation")
            elif option == "c":
                band_struc = BandStructureInterface(calc)
                band_struc.initialize_lcao_overlaps()
                band_struc.match_bands()
                band_struc.get_bands()
                chosen_option = True
                display_notif("Finished Calculation")
            elif option == "d":
                display_notif("Pending code")
            elif option == "e":
                exit()
            else:
                display_notif("Please choose a valid option")
        while chosen_option:
            print("\n>>>>\t What do you want to do now?\n\n")
            options = ["\ta. Display the plot of your band structure\n",
                       "\tb. Export your band structure data into a .txt file\n",
                       "\tc. Redo the calculations\n",
                       "\td. Close the program\n\nEnter your option >>> "]
            option = raw_input(''.join([options[i] for i in range(len(options))])).lower()
            if option == "a":
                band_struc.print_bands()
                chosen_option = True
                display_notif("Finished Calculation")
            elif option == "b":
                display_notif("Not yet implemented")
            elif option == "c": 
                chosen_option = False
                MainMenu().display_menu()
            elif option == "d":
                exit()
            else:
                print("Please chose a valid option")
                
    def set_options(self, op_list): #Pending method or displaying menus of varialbe options
        options = ["\ta. Finite difference(Real Space)\n",
                       "\tb. Plane Wave Expansion\n",
                       "\tc. Linear Combination of Atomic Orbitals (LCAO)\n",
                       "\td. Projector Augmented Wave\n\n>>> "]
        for i in range(len(op_list)):
            options.append("\t"+str(op_list[i])+"\n")
        print("\n>>> ")
        option = raw_input(''.join([options[i] for i in range(len(options))])).lower()
        
MENU = MainMenu()
MENU.display_menu()
