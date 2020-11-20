#This is a program to generate the input file for DD-vMCG

def dot_aligned(seq):
    snums = [str(n) for n in seq]
    dots = [s.find('.') for s in snums]
    m = max(dots)
    return [' '*(m - d) + s for s, d in zip(snums, dots)]

#ask the user to give input details
filename = input('Please give a filename : ')
sourcefile = input('Please give the filepath of the transfile: ')
qcprogram = input('Please give the name of the quantum chemistry program used to produce the transfile (gaussian/molpro/molcas): ')
ndof_in= input('Please give the number of degrees of freedom : ')

ndof = int(ndof_in)
j=0
jj=0
ener='0'
nbasis='0'
atom=list()
datax=list()
datay=list()
dataz=list()
dof=list()
sym=list()
freq=list()

if qcprogram == 'gaussian':
    with open(sourcefile) as sf:
        previous = None
        current = next(sf).strip()
        for i,line in enumerate(sf):
            line=line.strip()
            if 'Standard orientation:' in line:
                j=i
            if j!=0 and i>j+4:
                if '----------' in line:
                    j=0
                else:
                    lst=line.split()
                    atom.append(lst[1])
                    datax.append(lst[3])
                    datay.append(lst[4]) 
                    dataz.append(lst[5])
        
            if  line.startswith('Frequencies ---'):
                if len(dof) < ndof:
                    dof_lst=previous.split()
                    dof.extend(dof_lst)
                    sym_lst=current.split()
                    sym.extend(sym_lst)
                    freq_lst=line.split()
                    freq.extend(freq_lst[2:])
            previous = current
            current = line

            for part in line.split():        
                if "nbasis=" in part.lower():
                    nbasis = line.split()[line.split().index(part) + 1]
                if "zero-point" in part.lower():
                    next_part=line.split()[line.split().index(part) + 1]
                    if next_part ==  'Energies=':
                        ener = line.split()[line.split().index(part) + 2]
elif qcprogram == 'molcas':
    with open(sourcefile) as sf:
        previous = None
        current = next(sf).strip()
        for i,line in enumerate(sf):
            line=line.strip()
            if 'Cartesian coordinates in ' in line:
                j=i
            if j!=0 and i>j+3:
                if '----------' in line:
                    j=0
                else:
                    lst=line.split()
                    lst2=list(lst[1])
                    atom.append(lst2[0])
                    datax.append(lst[2])
                    datay.append(lst[3]) 
                    dataz.append(lst[4])
        
            if  line.startswith('Harmonic frequencies in '):
                jj=i
            if jj!=0 and len(dof) < ndof and line.startswith('Frequency:'):
                dof_lst=previous.split()
                dof.extend(dof_lst)
                sym_lst=current.split()
                sym.extend(sym_lst)
                freq_lst=line.split()
                freq.extend(freq_lst[1:])
            previous = current
            current = line

            if  line.startswith('Basis functions'):
                nbasis_lst=line.split()
                nbasis=(nbasis_lst[2])

            if  line.startswith('RASSCF energy for state  1 '):
                ener_lst=line.split()
                ener=(ener_lst[5])

elif qcprogram == 'molpro':
    with open(sourcefile) as sf:
        previous = None
        current = next(sf).strip()
        for i,line in enumerate(sf):
            line=line.strip()
            if line.startswith('Atomic Coordinates'):
                j=i
            if j!=0 and i>j+3:
                if len(line.strip()) == 0 :
                    j=0
                else:
                    lst=line.split()
                    atom.append(lst[1])
                    datax.append(lst[3])
                    datay.append(lst[4]) 
                    dataz.append(lst[5])
        
            if  line.startswith('Normal Modes'):
                jj=i
            if jj!=0 and len(dof) < ndof and line.startswith('Wavenumbers [cm-1]'):
                dof_lst=current.split()
                for n in range(len(dof_lst)):
                    if (n % 2) == 0:
                        dof.append(dof_lst[n])
                    else:
                        sym.append(dof_lst[n])
                freq_lst=line.split()
                freq.extend(freq_lst[2:])
            previous = current
            current = line

            if  line.startswith('NUMBER OF CONTRACTIONS:'):
                nbasis_lst=line.split()
                nbasis=(nbasis_lst[3])

            if  line.startswith('!MCSCF STATE 1.1 Energy'):
                ener_lst=line.split()
                ener=(ener_lst[4])                             

f1=open(filename+'.inp','w+' )
f1.write ('#' * 40 +'\n')
f1.write ('#' * 3 +  ' '*10+ 'DD propagation'+ ' '*10+ '#' * 3 + '\n')
f1.write ('#' * 40 +'\n')
f1.write ('\n')

f1.write ('RUN-SECTION\n')
f1.write ('name = '+filename+'_dd\n')
f1.write ('propagation direct = nmodes gwp= \n')
f1.write ('tfinal =     tout =     tpsi = \n' )
f1.write ('psi  gridpop   update  steps  auto \n' )
f1.write ('end-run-section\n')
f1.write ('\n')

f1.write ('INITIAL-GEOMETRY-SECTION\n')
f1.write ('nstates =      init_state = \n' )
f1.write ('cartesian = \n' )
for (a,b,c,d) in zip(atom,dot_aligned(datax),dot_aligned(datay),dot_aligned(dataz)):
    f1.write(a+'     '+b +'   '+c+'   '+d+'\n')
f1.write ('end-cartesian\n')
f1.write ('nmode\n')

for num in range(len(dof)):
    if not sym:
        line_new="{:<15}{:^15}{:<15} {:<15}".format(dof[num]+'A','0.0000',freq[num]+',','cm-1 width = 0.7071')
        f1.write(line_new+'\n')
    else:
        line_new="{:<15}{:^15}{:<15} {:<15}".format(dof[num]+sym[num],'0.0000',freq[num]+',','cm-1 width = 0.7071')
        f1.write(line_new+'\n')

f1.write ('end-nmode\n')
f1.write ('end-initial-geometry-section\n')
f1.write ('\n')



f1.write ('DIRDYN-SECTION\n')
f1.write ('data =      dd_diab=global \n')
sr=list(sourcefile.split('/'))
f1.write ('transfile = '+ sr[-1]+'\n')
f1.write ('update = always\n')
f1.write ('ener0 = '+ ener+'\n')
f1.write ('qcprogram =     method  = \n' )
f1.write ('db =      dbmin  = \n' )
f1.write ('nbasis = '+ nbasis+'\n')
f1.write ('subcmd = \n')
f1.write ('hess_upd  ddlog dbsave\n')
f1.write ('end-dirdyn-section\n')
f1.write ('\n')
f1.write ('end-input')