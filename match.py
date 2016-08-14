import numpy as np
import pdb

def match(A,B,test=False):
    '''
    Return indices of input arrays that correspond to matching values

    code modified from StackOverflow
    '''

    # Get sorted unique elements of A and B and the indices based on the uniqueness
    unqA,udx1,idx1 = np.unique(A,return_inverse=True,return_index=True)
    unqB,udx2,idx2 = np.unique(B,return_inverse=True,return_index=True)

    # Create mask equivalent to np.in1d(A,B) and np.in1d(B,A) for unique elements
    # mask will be true for matching elements
    # mask1 has all of the elements in unqA
    mask1 = (np.searchsorted(unqB,unqA,'right') - np.searchsorted(unqB,unqA,'left'))==1
    # mask2 has all of the elements in unqB
    mask2 = (np.searchsorted(unqA,unqB,'right') - np.searchsorted(unqA,unqB,'left'))==1

    # find the corresponding indices for unqA objects in unqB
    m=np.searchsorted(unqB,unqA,'left')
    # only select objects that are matched
    gd=np.where(mask1)[0]
    # get indices of these objects in original unsorted input array, B
    m2=udx2[m[gd]]

    # now do the same in the other direction
    m=np.searchsorted(unqA,unqB,'left')
    gd=np.where(mask2)[0]
    m1=udx1[m[gd]]

    if test :
        if len(m1) != len(m2) : print('length of index arrays does not match!')
        for i in range(len(m1)) :
            if A[m1[i]] != B[m2[i]] : print('failed! ' )

    return m1, m2

