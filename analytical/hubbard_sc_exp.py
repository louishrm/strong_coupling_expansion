import numpy as np 
from more_itertools import distinct_permutations
#from scipy.linalg import expm, block_diag, logm



class HubbardED:

    """Class for numerical diagonalization of the Hubbard model on the square lattice."""

    def __init__(self, n_orbitals:int):
        
        self.n_orbitals = n_orbitals #number of orbitals
        self.dim = 4**n_orbitals #dimension of the Hilbert space

        self.basis = self.all_bases() #dictionary with the basis states for each particle number and magnetization sector

        self.HU, self.Hmu = self.full_Hamiltonian_diag() #diagonal part of the Hamiltonian in each sector



    def get_basis_fixed_N_Sz(self, Nup, Ndown):

        """Gets the basis states in a symmetry sector with a fixed number of particles
        and a fixed magnetization."""

        if Nup > self.n_orbitals or Ndown > self.n_orbitals:

            raise ValueError('Number of particles in each spin sector cannot exceed the number'
                             'of orbitals.')

        shuffle_array_up = np.array([1]*Nup + [0]*(self.n_orbitals-Nup))
        shuffle_array_down = np.array([1]*Ndown + [0]*(self.n_orbitals-Ndown))


        up_basis = np.array(list(distinct_permutations(shuffle_array_up)))
        down_basis = np.array(list(distinct_permutations(shuffle_array_down)))

        basis = {}; counter=0


        for up in up_basis:
             for down in down_basis:
                basis_state = np.concatenate((down,up))
                basis_state_int = int(''.join(map(str, basis_state)), 2)
                basis[basis_state_int] = counter; counter += 1


        return basis
    

    def all_bases(self):

        """Gets all the possible basis states, each dictionary contains the basis for a fixed particle number and magnetization sector."""

        all_bases = []

        for Ndown in range(self.n_orbitals+1):
            for Nup in range(self.n_orbitals+1):

                if Ndown == 0 and Nup == 0: #forget about the 0 particle space. 
                    continue


                all_bases.append(self.get_basis_fixed_N_Sz(Nup,Ndown))

        return all_bases



    def ns(self, state: int,i:int,s:int):
        """returns the eigenvalue for the number operator for spin sigma (0: down, 1:up) on site i
        either 0 or 1."""

        return (state>>(i+self.n_orbitals*s))&1


    def phase(self,state: int,i:int,j:int,s:int):
        """returns the phase for the creation operator c^dagger_i c_j +/-1"""

        parity = sum([(state>>(l+self.n_orbitals*s))&1 for l in range(min(i,j)+1,max(i,j))])%2
        return (-1)**parity


    def checkcdagc(self,state:int,i:int,j:int,s:int):
        """returns True if the hopping is possible, False otherwise"""

        if self.ns(state,j,s) - self.ns(state,i,s) ==1:
            return True
        else:
            return False


    def cdagc(self,state:int,i:int,j:int,s):
        """Returns a new state with a particle created on site i and destroyed on site j as well as the phase, 
        assumes that check is True."""

        p = self.phase(state,i,j,s)
        mask = (1<<(i+self.n_orbitals*s)) | (1<<(j+self.n_orbitals*s))
        newstate = state^mask
        return p, newstate




    def Hamiltonian_Diag(self, basis:dict):
        
        """Constructs the diagonal part (U) and (mu) parts of the Hubbard Hamiltonian in the basis corresponding to a fixed particle number and magnetization sector."""

        dim = len(basis)
        HU = np.zeros((dim,dim),dtype=complex)
        Hmu = np.zeros((dim,dim),dtype=complex)

        for state, col in basis.items():

            #interaction+chemical potential
            for i in range(self.n_orbitals):

                HU[col,col] += self.ns(state,i,0)*self.ns(state,i,1) 

                Hmu[col,col] += -(self.ns(state,i,0) + self.ns(state,i,1))

        return HU, Hmu
    

    def full_Hamiltonian_diag(self):

        """Get the full diagonal part of the Hamiltonian, by block diagonalizing it in each 
        particle number and magnetization sector."""

        HUs, Hmus = [], []

        for basis in self.basis:
            HU, Hmu = self.Hamiltonian_Diag(basis)
            HUs.append(HU); Hmus.append(Hmu)

        return HUs, Hmus


    def Hamiltonian_Hopping(self, neighbors_and_t: list, basis: dict):
        """Constructs the hopping part of Hubbard Hamiltonian in the basis corresponding 
        to a fixed particle number and magnetization sector."""


        dim = len(basis)
        H = np.zeros((dim,dim),dtype=complex)

        for state,col in basis.items():
            #hoppings
            for i,j,t in neighbors_and_t:
                for s in range(2):

                    if self.checkcdagc(state,i,j,s): #hopping j to i with spin s
                        p, new_state = self.cdagc(state,i,j,s)
                        row =  basis[new_state]
                        H[row,col] += -t*p
                        H[col,row] += -t*p        

        return H
    

    def full_Hamiltonian_Hopping(self, neighbors_and_t: list):

        """Get the full hopping part of the Hamiltonian, by block diagonalizing it in each 
        particle number and magnetization sector."""

        Hs = []

        for basis in self.basis:
            H = self.Hamiltonian_Hopping(neighbors_and_t,basis)
            Hs.append(H)

        return Hs
    


    def full_spectrum(self, U: float, mu: float, neighbors_and_t: list):

        """Get the full spectrum of the Hamiltonian by diagonalizing it in each particle number and magnetization sector."""

        energies = [np.array([0])]

        Hts= self.full_Hamiltonian_Hopping(neighbors_and_t)

        for i in range(len(self.basis)):
            
            H = U*self.HU[i] + mu*self.Hmu[i] + Hts[i]
            E = np.linalg.eigvals(H)
            energies.append(E)

        return np.concatenate(energies)
        
    

    def partition_function(self, rho: np.ndarray):

        """Get the partition function of the Hubbard model at inverse temperature beta."""
        Z = np.trace(rho)

        return Z
    
                
        


