import numpy as np
from astropy.table import Table


print("Programme produced by Max Howe, email: mh2003@cam.ac.uk \n \n",
      "Commands: \n",  
      "linear(n) - returns the Hückel π-energies and degeneracies for a linear polyene with n carbons. \n",
      "cyclic(n) - returns the Hückel π-energies and degeneracies for a cyclic polyene with n carbons. \n",
      "platonic(n) - returns the Hückel π-energies and degeneracies for the platonic solid with n carbons. \n \n",
      "Note: Buckminsterfullerene (n=60) is a truncated icosahedron and can be investigated like the platonic solids. \n")


value_err_msg = "Only non-negative integers are permitted in the bond chain."


## Task 1: Determine the Hückel π-energies and degeneracies for a linear polyene with n carbons

## An input of 'linear(n)', where n is the number of carbon atoms, will give the desired results
def linear(n):

    ## n must be a non-negative integer
    if n != int(n):
        print(value_err_msg)
        return
    elif n <= 0:
        print(value_err_msg)
        return

    ## Set the matrix size
    Hückel_matrix = np.zeros([n,n])
    evals = np.empty(shape=(0))
    evects = np.empty(shape=(0, 0))

    ## Here we have chosen that α=0 and β=-1. So our energies are relative to α and in terms of |β|
    ## Iterate through the molecule and set matrix elements
    i = 0
    for j in range(n):
        if  j == i+1:
             Hückel_matrix[j , i ] = -1 
             Hückel_matrix[i , j ] = -1
        i=j

    ## Compute the eigenvalues and eigenvectors of our square array
    ## The eigenvalues are the π-energies of the system and the eigenvectors give the orbital coefficients on each atom
    evals, evects = np.linalg.eig(Hückel_matrix)

    ## Define a degeneracy vector whose elements correspond to the degeneracies of the eigenvectors
    ## No need for rigorous degeneracy analysis here as none of the MOs are degenerate 
    degeneracy = np.ndarray(n)
    for i in range(n):
        degeneracy[i] = sorted(evals).count(sorted(evals)[i])
        
    ## Make a table of results showing the energies and degeneracies of the π MOs formed
    ## Rounding eigenvalues makes them easier to interpret
    ## Note this works up to about n=2000, above which the eigenvalues are too close together in energy so more significant figures should be used
    t = Table([sorted(np.round(evals,5)),degeneracy], names = ('Relative Energy /|β|','Degeneracy'))
    print(t)


## Task 2: Determine the Hückel π-energies and degeneracies for a cyclic polyene with n carbons

## An input of 'cyclic(n)', where n is the number of carbon atoms, will give the desired results
def cyclic(n):

    ## n must be a non-negative integer
    if n != int(n):
        print(value_err_msg)
        return
    elif n <= 0:
        print(value_err_msg)
        return

    ## Set the matrix size
    Hückel_matrix = np.zeros([n,n])
    evals = np.empty(shape=(0))
    evects = np.empty(shape=(0, 0))
    unique_evals = []

    ## Here we have chosen that α=0 and β=-1. So our energies are relative to α and in terms of |β|
    ## Iterate through the molecule and set matrix elements
    i = 0
    for j in range(n):
        if  j == i+1:
             Hückel_matrix[j , i ] = -1 
             Hückel_matrix[i , j ] = -1
        i=j
        ## Same matrix as before but with top right and bottom left elements also having a value of -1, due to end connections       
        Hückel_matrix[n-1 , 0 ] = -1 
        Hückel_matrix[0 , n-1 ] = -1

    ## Compute the eigenvalues and eigenvectors of our square array
    ## The eigenvalues are the π-energies of the system and the eigenvectors give the orbital coefficients on each atom
    evals, evects = np.linalg.eig(Hückel_matrix)

    ## Now we do have degeneracies so to only print each eigenvalue once in our table we need a new list of unique eigenvalues
    ## To get degeneracies, eigenvalues need to be rounded as otherwise due to computer errors the exact values of the eigenvalues won't be the same

    e = np.round(evals,5)

    ## Set up a list of unique eigenvalues
    def unique(e):
        for x in e:
            if x not in unique_evals:
                unique_evals.append(x)
    unique(e)
    
    ## Define a degeneracy vector whose elements correspond to the degeneracies of the eigenvectors
    degeneracy = np.ndarray(len(unique_evals))
    for i in range(len(unique_evals)):
        ## It is necessary to round the energies as otherwise due to computer limitations in the accuracy of the values, the degeneracies will not be shown 
        degeneracy[i] = sorted(np.real(np.round(evals,5))).count(sorted(np.real(unique_evals))[i])
        
    ## Make a table of results showing the energies and degeneracies of the π MOs formed.
    ## Rounding eigenvalues makes them easier to interpret
    ## Note this works up to about n=2000, above which the eigenvalues are too close together in energy so more significant figures should be used
    t = Table([sorted(np.real(unique_evals)),degeneracy], names = ('Relative Energy /|β|','Degeneracy'))
    print(t)


## Tasks 3 & 4: Determine the Hückel π-energies and degeneracies for the platonic solids and buckminsterfullerene

## An input of 'platonic(n)', where n is the number of carbon atoms in the molecule, will give the desired results
def platonic(n):

    ## Define the connectivities of the molecules in terms of the number of atoms they contain
    ## Tetrahedron
    if n == 4:
        m = (1, 2, 3, 4, 2, 0, 3, 1, 4)
    ## Cube    
    elif n == 8:
        m = (1, 2, 3, 4, 1, 5, 6, 7, 8, 5, 0, 2, 6, 0, 3, 7, 0, 4, 8)
    ## Dodecahedron
    elif n == 20:
        m = (1, 2, 3, 4, 5, 1, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 6, 0, 2, 8, 0, 3, 10, 0, 4, 12,
             0, 5, 14, 0, 7, 16, 17, 18, 19, 20, 16, 0, 9, 17, 0, 11, 18, 0, 13, 19, 0, 15, 20)
    ## Octahedron
    elif n == 6:
        m = (1, 2, 3, 4, 5, 2, 6, 3, 1, 4, 6, 5, 1)
    ## Icosahedron
    elif n == 12:
        m = (1, 3, 5, 1, 0, 2, 4, 6, 2, 0, 1, 8, 3, 10, 5, 12, 1, 0, 2, 9, 4, 11, 6, 7, 2, 8, 0, 3,
             9, 0, 4, 10, 0, 5, 11, 0, 6, 12, 1, 7, 8, 9, 10, 11, 12, 7)
    ## Buckminsterfullerene
    elif n == 60:
        m = (1, 2, 3, 4, 5, 1, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 6, 0, 2, 9,
             0, 3, 12, 0, 4, 15, 0, 5, 18, 0, 7, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33,
             34, 35, 36, 37, 38, 39, 40, 21, 22, 0, 8, 25, 0, 10, 26, 0, 11, 29, 0, 13, 30, 0, 14,
             33, 0, 16, 34, 0, 17, 37, 0, 19, 38, 0, 20, 21, 0, 40, 41, 42, 43, 44, 45, 46, 47, 48,
             49, 50, 51, 52, 53, 54, 55, 41, 0, 23, 42, 0, 24, 44, 0, 27, 45, 0, 28, 47, 0, 31, 48,
             0, 32, 50, 0, 35, 51, 0, 36, 53, 0, 39, 54, 0, 55, 56, 57, 58, 59, 60, 56, 0, 43, 57,
             0, 46, 58, 0, 49, 59, 0, 52, 60)

    value_err_msg_2 = "Only n corresponding to the number of atoms in a platonic solid (or buckminsterfullerene) are permitted"
    ## n must correspond to the number of atoms in one of the platonic solids
    if n not in set([4, 8, 20, 6, 12, 60]):
        print(value_err_msg_2)
        return

    ## Turns a bond chain into a Huckel matrix
    def chain(m):

        ## Set the matrix size
        Hückel_matrix = np.zeros([n,n])
        evals = np.empty(shape=(0))
        evects = np.empty(shape=(0, 0))
        unique_evals = []

        # Iterate through the chain and set matrix elements.
        i = 0
        for j in m:
            if i != 0 and j != 0 and j != i:
                Hückel_matrix[j-1 , i-1 ] = -1
                Hückel_matrix[i-1 , j-1 ] = -1
            i = j

        ## Compute the eigenvalues and eigenvectors of our square array
        ## The eigenvalues are the π-energies of the system and the eigenvectors give the orbital coefficients on each atom
        evals, evects = np.linalg.eig(Hückel_matrix)

        ## Set up a list of unique eigenvalues
        def unique(e):
            for x in e:
                if x not in unique_evals:
                    unique_evals.append(x)

        ## We do have degeneracies so to only print each eigenvalue once in our table we need a new list of unique eigenvalues
        ## To get degeneracies, eigenvalues need to be rounded as otherwise due to computer errors the exact values of the eigenvalues won't be the same
        e = np.round(evals,5)
        unique(e)
    
        ## Define a degeneracy vector whose elements correspond to the degeneracies of the eigenvectors
        degeneracy = np.ndarray(len(unique_evals))
        for i in range(len(unique_evals)):
            ## Here it is necessary to round the energies as otherwise due to computer limitations in the accuracy of the values, the degeneracies will not be shown 
            degeneracy[i] = sorted(np.real(np.round(evals,5))).count(sorted(np.real(unique_evals))[i])
        ## Make a table of results showing the energies and degeneracies of the π MOs formed.
        ## Rounding eigenvalues makes them easier to interpret
        t = Table([sorted(np.real(unique_evals)),degeneracy], names = ('Relative Energy /|β|','Degeneracy'))
        print(t)

    chain(m)
