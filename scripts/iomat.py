##########################################################################################
# This script does the following:
# - Constructs the IO matrices from processed WIOD data files
# - Writes the matrices to csv files and latex tables
##########################################################################################

##########################################################################################
# Imports, constants, etc.
##########################################################################################
import pandas as pd
import numpy as np
import itertools
import locale
import regions

locale.setlocale(locale.LC_ALL,'en_US.utf8')

nc=regions.nr
year=2011
countries=['UK','EU','ROW']

inpath = 'output/'
outpath = 'output/'

##########################################################################################
# Function definitions
##########################################################################################

def read_wiod_data(year):
    # Reads the WIOD data for a given year and returns three Pandas DataFrames:

    va = pd.read_csv(inpath + 'va' + str(year) +'.txt',sep=' ')
    inin = pd.read_csv(inpath + 'inin' + str(year) +'.txt',sep=' ')
    fin = pd.read_csv(inpath + 'fin' + str(year) +'.txt',sep=' ')
    va=va[pd.notnull(va['i'])]
    inin=inin[pd.notnull(inin['i'])]
    fin=fin[pd.notnull(fin['i'])]
    return va, inin, fin

def read_iomat(fname):
    # Reads an IO matrix from text file and returns a numpy matrix

    iomat = np.genfromtxt(fname=fname,dtype='float',delimiter=',',names=None)
    colsums=iomat[-1,:]
    rowsums=iomat[:,-1]

    return iomat, rowsums, colsums

def construct_iomat(va,inin,fin):
    # Given three Pandas DataFrames that contain value added (va), intermediate inputs (inin),
    # and final demand (fin), constructs a world input-output matrix. Returns three Numpy arrays:

    va=va.sort_values(by=['i'])
    fin=fin.sort_values(by=['i','j'])
    inin=inin.sort_values(by=['i','j'])

    rowsums = np.zeros( nc +1 )
    colsums = np.zeros( nc + nc )

    M = inin.pivot_table(values='m', index=['j'], columns=['i']).values
    V = va['va'].values.reshape((1,nc))
    F = fin.pivot_table(values=['f'], index=['j'], columns='i').values
    V = np.hstack((V,np.zeros((1,nc))))
    iomat=np.vstack( ( np.hstack((M,F)) , V ) )

    for row in range(0,nc + 1):
        rowsums[row] = np.sum(iomat[row,:])

    for col in range(0,nc + nc):
        colsums[col] = np.sum(iomat[:,col])

    return iomat, rowsums, colsums

def coeffs(iomat):
    # Given world IO matrix (iomat), calculates IO coefficients and returs them in A

    A=np.zeros(iomat.shape)
    for col in range(0,A.shape[1]):
        A[:,col] = iomat[:,col]/np.sum(iomat[:,col])
    return A

def ras(iomat0,rowsums1,colsums1):
    # Given an initial IO matrix (iomat), and desired rowsums (rowsums1) and colsums (colsums1),
    # performs the RAS balancing procedure. Returns a new IO matrix (iomat) that is consistent
    # with desired row- and colsums.

    A0 = coeffs(iomat0)
    iomat = np.dot(A0,np.diag(colsums1))

    go=True
    iter=0
    maxit=10000
    tol=1.0e-8

    while go:
        iter=iter+1
        rowsums = np.sum(iomat,axis=1)
        r = np.divide(rowsums1,rowsums)
        iomat = np.dot(np.diag(r),iomat)
        colsums = np.sum(iomat,axis=0)
        s = np.divide(colsums1,colsums)
        iomat = np.dot(iomat,np.diag(s))
        colsums = np.sum(iomat,axis=0)
        rowsums = np.sum(iomat,axis=1)

        norm1 = max(np.divide(abs(rowsums-rowsums1),rowsums1))
        norm2 = max(np.divide(abs(colsums-colsums1),colsums1))
        if((norm1 <tol and norm2 <tol) or iter == maxit):
            go=False

    if iter==maxit:
        print 'RAS iteration did not converge!'
        print 'iter = ', iter, ' diff = ', max(norm1,norm2)
    else:
        print 'RAS converged after ',str(iter),' iterations'


    return iomat


def write_iomat_csv(iomat,fname):
    # Write world IO matrix (iomat) to csv file called filename

    iomat2 = np.vstack((iomat,np.sum(iomat,axis=0).reshape((1,nc+nc))))
    iomat2 = np.hstack((iomat2,np.sum(iomat2,axis=1).reshape((nc+2,1))))
    iomat2 = 100*iomat2/iomat2[nc,0]
    np.savetxt(fname=outpath+fname,X=iomat2,fmt='%0.15f',delimiter=',')
    
def write_iomat_latex(iomat,rowsums,colsums,caption,label,fname):
    # Given a world IO matrix (iomat), rowsums, colsums, creates a latex file
    # in location filename that contains a table with given caption and label.

    ukgdp = iomat[nc,0].sum()
    iomat2 = 100*iomat[:,:]/ukgdp
    rowsums2 = 100*rowsums/ukgdp
    colsums2 = 100*colsums/ukgdp
    M=iomat2[0:nc,0:nc]
    V=iomat2[-2,0:nc]
    F=iomat2[0:nc,nc:(2*nc)]

    with open('tex/' + fname + '.tex','wb') as file:

        file.write('\\begin{table}[p]\n')
        file.write('\\begin{center}\n')
        file.write('\\caption{'+caption+'}\n')
        file.write('\\label{tab:'+label+'}\n')

        # number of columns
        file.write('\\begin{tabular}{c')
        for i in range(0,nc): # intermediate inputs columns
            file.write('c')
        for i in range(0,nc): # final demand columns
            file.write('c')
        file.write('c}\n')
        file.write('\\toprule\n')

        file.write('&\\multicolumn{'+str(nc)+'}{c}{Intermediate inputs}&\\multicolumn{'+str(nc)+'}{c}{Final demand}&\\\\\n')
        file.write('\\cmidrule(rl){2-4}\\cmidrule(rl){5-7}\n')

        # middle headers: country names (underline them)
        for c in countries:
            file.write('& '+c)
        for c in countries:
            file.write('& '+c)
        file.write('& GO\\\\\n')

        file.write('\\midrule\n')

        # write intermediate input and final use data
        for i in range(0,nc):
            file.write(countries[i])
            for j in range(0,nc):
                tmpstr = locale.format('%0.1f',M[i,j],grouping=True)
                file.write('&'+tmpstr)
            for j in range(0,nc):
                tmpstr = locale.format('%0.1f', F[i,j],grouping=True)
                file.write('&'+tmpstr)

            tmpstr=locale.format('%0.1f',rowsums2[i],grouping=True)
            file.write('&'+tmpstr)
            file.write('\\\\\n') 

        # write value added portion
        file.write('\\midrule\n')
        file.write('VA')
        for i in range(0,nc):
            tmpstr = locale.format('%0.1f',V[i],grouping=True)
            file.write('&'+tmpstr)
        for i in range(0,nc):
            file.write('& - ')
        tmpstr = locale.format('%0.1f',sum(V),grouping=True)
        file.write('&'+tmpstr)
        file.write('\\\\\n')

        # write gross output portion
        file.write('\\midrule\n')
        file.write('GO')
        for i in range(0,nc):
            tmpstr = locale.format('%0.1f',colsums2[i],grouping=True)
            file.write('&'+tmpstr)
        for i in range(0,nc):
            tmpstr = locale.format('%0.1f',colsums2[nc+i],grouping=True)
            file.write('&'+tmpstr)
        file.write('&\\\\\n')
            
        file.write('\\bottomrule\n')
        file.write('\\end{tabular}\n')
        file.write('\\end{center}\n')
        file.write('\\end{table}\n')

##########################################################################################
# Main script body
##########################################################################################

# load the preprocessed data and convert to an IO matrix
va, inin, fin = read_wiod_data(year)
iomat, rowsums, colsums = construct_iomat(va, inin, fin)

# make sure it's balanced
iomat = ras(iomat,rowsums,colsums)

# balance the matrix so that all countries have balanced trade, then write to CSV
iomatb = np.zeros((nc*2,nc*2))
iomatb[0,:]=iomat[0,:]
iomatb[1,:]=iomat[1,:]
iomatb[2,:]=iomat[2,:]
iomatb[3,0]=iomat[3,0]
iomatb[4,1]=iomat[3,1]
iomatb[5,2]=iomat[3,2]

rowsumsb=np.zeros(nc*2)
rowsumsb[0]=rowsums[0]
rowsumsb[1]=rowsums[1]
rowsumsb[2]=rowsums[2]
rowsumsb[3]=iomat[3,0]
rowsumsb[4]=iomat[3,1]
rowsumsb[5]=iomat[3,2]

colsumsb=colsums
colsumsb[nc:(2*nc)]=iomat[nc,0:nc]

iomatb = ras(iomatb,rowsumsb,colsumsb)
iomatb2 = np.zeros(iomat.shape)
iomatb2[0,:]=iomatb[0,:]
iomatb2[1,:]=iomatb[1,:]
iomatb2[2,:]=iomatb[2,:]
iomatb2[3,0]=iomatb[3,0]
iomatb2[3,1]=iomatb[4,1]
iomatb2[3,2]=iomatb[5,2]

# write the raw and balanced matrices to csv and LaTeX
write_iomat_csv(iomat,'iomat'+str(year)+'.csv')
iomat, rowsums, colsums = read_iomat(outpath+'iomat'+str(year)+'.csv')
write_iomat_latex(iomat,
                  rowsums,
                  colsums,
                  str(year)+' raw world input-output table (UK GDP = 100)',
                  'iomat'+str(year),
                  'iomat'+str(year))


write_iomat_csv(iomatb2,'iomat_balanced'+str(year)+'.csv')
iomat, rowsums, colsums = read_iomat(outpath+'iomat_balanced'+str(year)+'.csv')
write_iomat_latex(iomat,
                  rowsums,
                  colsums,
                  str(year)+' balanced world input-output table (UK GDP = 100)',
                  'iomat_balanced'+str(year),
                  'iomat_balanced'+str(year))
