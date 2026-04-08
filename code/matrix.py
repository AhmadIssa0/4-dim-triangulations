
class Matrix:
    def __init__(self, mat=None):
        """A matrix is given by a list of rows (a row being a list of integers)."""
        if mat is not None:
            self.mat = mat
        else:
            self.mat = []
            
    def switch(self, i, j):
        self.mat[i], self.mat[j] = self.mat[j], self.mat[i]
        
    def colSwitch(self, i, j):
        A = self.mat
        for k in range(len(A)):
            A[k][i], A[k][j] = A[k][j], A[k][i]
            
    def colAddMultiple(self, c, i, j):
        A = self.mat
        for k in range(len(A)):
            A[k][j] += c*A[k][i]
        
    def addMultiple(self, c, i, j):
        # row j -> c*(row i) + (row j)
        A = self.mat
        for k in xrange(len(A[j])):
            if A[i][k] != 0:
                A[j][k] += c*A[i][k]
            
    def scaleRow(self, c, i):
        self.mat[i] = [x*c for x in self.mat[i]]
        
    def scaleCol(self, c, i):
        A = self.mat
        for k in range(len(A)):
            A[k][i] = c*A[k][i]
        
    def adjoinIdentity(self):
        num_rows = len(self.mat)
        num_cols = len(self.mat[0])
        n_mat = []
        for i in range(num_rows):
            n_mat.append(self.mat[i][:] + [0]*i + [1] + [0]*(num_rows-i-1))
        return Matrix(n_mat)
        
    def inverse(self):
        # inverts, only works if invertible
        # otherwise may behave badly
        return self.RRDecomposition()[1]
        
    def __str__(self):
        return str(self.mat)
        
    def multiply(self, m2):
        A = self.mat
        A_rows = len(A)
        
        B = m2.mat
        B_cols = len(B[0])
        
        res = []
        for i in range(A_rows):
            res.append([0]*B_cols)
            for j in range(B_cols):
                for k in range(len(A[0])):
                    res[i][j] += A[i][k]*B[k][j]
        return Matrix(res)
        
    def transpose(self):
        n_rows = len(self.mat)
        n_cols = len(self.mat[0])
        n_mat = []
        for i in range(n_cols):
            n_mat.append([])
            for j in range(n_rows):
                n_mat[i].append(self.mat[j][i])
        return Matrix(n_mat)
        
    def RRDecomposition(self):
        # decomposition of matrix into A = U R
        # where U is unimodular and R is row reduced
        # if self is unimodular then U is its inverse
        n_rows = len(self.mat)
        n_cols = len(self.mat[0])
        adj_mat = self.adjoinIdentity()
        adj_mat.rowReduce(n_cols)
        
        R = [x[:n_cols] for x in adj_mat.mat]
        U = [x[n_cols:] for x in adj_mat.mat]
        return (Matrix(R), Matrix(U))
        
    def dim(self):
        return str(len(self.mat)) + ' x ' + str(len(self.mat[0]))
        
    def clone(self):
        c_mat = []
        c_mat = [x[:] for x in self.mat]
        return Matrix(c_mat)
        
    @staticmethod
    def zeros(m, n):
        mat = []
        for i in range(m):
            mat.append([0]*n)
        return Matrix(mat)
        
    @staticmethod
    def identity(m):
        mat = []
        for i in range(m):
            mat.append([0]*i + [1] + [0]*(m-i-1))
        return Matrix(mat)

    # Note this isn't quite smith normal form.
    # We only put the matrix into diagonal form.
    # We do not ensure that d_i divides d_{i+1},
    # where d_j is the jth diagonal element.
    def smithNormalForm(self):
        c_mat = self.clone() # copy
        A = c_mat.mat
        n_rows = len(A)
        n_cols = len(A[0])
        
        # we want to write D = P A Q
        # where D is diagonal (as much as possible)
        P = Matrix.identity(n_rows)
        Q = Matrix.identity(n_cols)
        
        for p in range(min(n_cols, n_rows)):
            # p is the pivot we're dealing with
            
            # make sure pivot is non-zero
            if A[p][p] == 0:
                # search rows for non-zero entry
                for j in range(p+1,n_rows):
                    if A[j][p] != 0:
                        c_mat.switch(j,p)
                        P.switch(j,p)
                        break
                       
                if A[p][p] == 0:
                    for j in range(p+1,n_cols):
                        if A[p][j] != 0:
                            c_mat.colSwitch(j,p)
                            Q.colSwitch(j,p)
                            
                if A[p][p] == 0:
                    continue
                    
            if A[p][p] < 0:
                c_mat.scaleRow(-1, p)
                P.scaleRow(-1, p)
                
            finished = False
            while not finished:
                # make pivot gcd of the column
                # and eliminate all other entries in column
                finished = True
                for i in range(p+1,n_rows):
                    while A[i][p] != 0:
                        finished = False
                        if A[i][p] < 0:
                            c_mat.scaleRow(-1, i)
                            P.scaleRow(-1, i)
                        if A[i][p] < A[p][p]:
                            c_mat.switch(i,p)
                            P.switch(i,p)
                        c = A[i][p] / A[p][p]
                        c_mat.addMultiple(-1*c, p, i)
                        P.addMultiple(-1*c, p, i)
                        
                # make pivot gcd of the row
                # eliminate all other entries in row
                for i in range(p+1,n_cols):
                    while A[p][i] != 0:
                        finished = False
                        if A[p][i] < 0:
                            c_mat.scaleCol(-1, i)
                            Q.scaleCol(-1, i)
                        if A[p][i] < A[p][p]:
                            c_mat.colSwitch(p,i)
                            Q.colSwitch(p,i)
                        c = A[p][i] / A[p][p]
                        c_mat.colAddMultiple(-1*c, p, i)
                        Q.colAddMultiple(-1*c, p, i)
                
        return (c_mat, P, Q)
            
    def rowReduce(self, cols=None):
        # cols tells us how many columns to row reduce
        n_rows = len(self.mat)
        n_cols = len(self.mat[0])
        A = self.mat
        
        p_i = 0 # pivot column index
        for p_j in range(cols):
            for i in range(p_i,n_rows):
                if A[i][p_j] != 0:
                    if i > p_i:
                        self.switch(i,p_i)
                    break
            if A[p_i][p_j] == 0: # this column has no pivot
                continue
            if A[p_i][p_j] < 0:
                self.scaleRow(-1, p_i)
            for i in range(p_i+1,n_rows):
                while A[i][p_j] != 0:
                    if A[i][p_j] < 0:
                        self.scaleRow(-1, i)
                    if A[i][p_j] < A[p_i][p_j]:
                        self.switch(i,p_i)
                    c = A[i][p_j] / A[p_i][p_j]
                    self.addMultiple(-1*c, p_i, i)
            ## eliminate entries above ##
            for i in range(0, p_i):
                if A[i][p_j] % A[p_i][p_j] == 0:
                    c = A[i][p_j] / A[p_i][p_j]
                    self.addMultiple(-1*c, p_i, i)
            ####
            p_i += 1
            

#mat = Matrix([[2, -1, 0, 0, 0, 0, 0, 0, 0, -1], [-1, 2, -1, 0, 0, 0, 0, 0, 0, 0], [0, -1, 3, -2, 0, 0, 0, 0, 0, 0], [0, 0, -2, 2, 1, 0, 0, 0, 0, 0], [0, 0, 0, 1, -2, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 2, -1, 0, 0, -1], [0, 0, 0, 0, 0, -1, 2, -1, 0, 0], [0, 0, 0, 0, 1, 0, -1, 1, -1, 0], [0, 0, 0, 0, 0, 0, 0, -1, 1, 1], [-1, 0, 0, 0, 0, -1, 0, 0, 1, 1]])
#print mat.smithNormalForm()[0]
"""
mat = Matrix([[-2,1,0,0],[1,-2,1,0],[0,1,-2,2],[0,0,2,-4]])
mat = Matrix([[-2,0,1],[0,-2,1],[1,1,-2]])
mat = Matrix([[-2,1,1,1,1],[1,2,0,0,0],[1,0,3,0,0],[1,0,0,-3,0],[1,0,0,0,-4]])
print mat.smithNormalForm()[0]
"""
"""
R, U = mat.RRDecomposition()
print 'R=', R
print 'U=', U
"""
