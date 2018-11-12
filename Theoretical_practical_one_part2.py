#Imports
import numpy as np
import math
import matplotlib.pyplot as plt
import pandas as pd

#Function Definitions

#Takes in a list of connections of the molecule, starting at 0. Eg cyclopropane would be [0,1,2,0], a tetrahedral cage would be [0,1,2,3,1,0,2] etc

def atom_bridge_list(list1):
    matA = np.zeros((max(list1)+1,max(list1)+1))

    for i in range(len(list1)):
        if i-1>=0:
            tupl1 = tuple(sorted([list1[i-1],list1[i]]))
        
            matA[tupl1] = -1

    Huckelmat = matA + matA.T

    evals, evecs = np.linalg.eig(Huckelmat)
    evals = evals[evals.argsort()]
    
    bonds = list(matA.reshape((max(list1)+1)**2)).count(-1)
    
    return evals,bonds

#Special case for 2x diatomics

def diatomics():
    matA = np.zeros((4,4))

    matA[1,0] = -1
    matA[3,2] = -1

    Huckelmat = matA + matA.T

    evals, evecs = np.linalg.eig(Huckelmat)
    evals = evals[evals.argsort()]
    
    bonds = 2
    
    return evals,bonds

def Huckeldiatomics(bridge_list):
    
    bonds = 2
    atoms = 4
    
    
    
    
    evals = diatomics()[0]
    electronlist, energylist = Energylist(evals,bridge_list)
    moment_list = nth_moment(evals,bridge_list)
    modified_evals = modify_energies(atoms,bonds,evals)
    modified_energy_list = Energylist(modified_evals,bridge_list)[1]
    return electronlist,energylist,modified_energy_list,evals,moment_list



#MOfiller fills up the orbitals and spits out the total energy

def MOfiller(eval_set,elec):
    orbitals_occupied = elec/2
    
    tot_energy=0

    for i in range(math.ceil(orbitals_occupied)):
        tot_energy+=2*eval_set[i]
    
    if elec%2==1:
        tot_energy-=eval_set[math.ceil(orbitals_occupied)-1]

    return round(tot_energy,3)


#Energylist returns the orbital energies

def Energylist(eval_set,list):

    energylist = []
    electronlist = []
    for j in range(1,2*(max(list)+1)+1):
        energylist.append(MOfiller(eval_set,j))
        electronlist.append(j)
    
    return electronlist,energylist
        

#mewP returns the pth moment

def mewP(p,evals,list1):
    mewp = 0

    for l in range(max(list1)+1):
        mewp += evals[l]**p
    return mewp


#nth moment creates a list of the moments

def nth_moment(evals,list):
    return [round(mewP(p,evals,list),3) for p in range(max(list)+1)]
    


#modifyenergies returns the modified energies of the molecules once taking into account nuclear repulsion

def modify_energies(atoms,bonds,evals):
    
    return [i*np.sqrt(atoms/(2*bonds)) for i in evals]
    





#Huckel spits out all the important data to create the images and tables

def Huckel(bridge_list):
    
    bonds = atom_bridge_list(bridge_list)[1]
    atoms = max(bridge_list)+1
    
    
    
    
    evals = atom_bridge_list(bridge_list)[0]
    electronlist, energylist = Energylist(evals,bridge_list)
    moment_list = nth_moment(evals,bridge_list)
    modified_evals = modify_energies(atoms,bonds,evals)
    modified_energy_list = Energylist(modified_evals,bridge_list)[1]
    return electronlist,energylist,modified_energy_list,evals,moment_list
    


#Creating the data

three_atom_linear = Huckel([0,1,2])
three_atom_cyclic = Huckel([0,1,2,0])
four_atom_linear = Huckel([0,1,2,3])
four_atom_cyclic = Huckel([0,1,2,3,0])
four_atom_rhombus= Huckel([0,1,2,3,0,3,1])
four_atom_tetrahedron= Huckel([0,1,2,3,0,3,1,0,2])
four_atom_roadsign= Huckel([0,1,2,3,1])
two_diatomics= Huckeldiatomics([0,1,2,3])


#Plotting the 4 graphs

fig_adjusted = plt.figure(figsize = (8,5),dpi = 1000)

axes = fig_adjusted.add_axes([0.1,0.1,0.8,0.8])

axes.plot(four_atom_cyclic[0],four_atom_cyclic[2],'+-',markersize = 10,label='Cyclic')
axes.plot(four_atom_linear[0],four_atom_linear[2],'o-',markersize = 5,markerfacecolor = 'white',label='Linear')
axes.plot(four_atom_rhombus[0],four_atom_rhombus[2],'^-',markerfacecolor=(1,1,1,0.9),label='Rhombus')
axes.plot(four_atom_tetrahedron[0],four_atom_tetrahedron[2],'s-',markerfacecolor = (1,1,1,0.9),label='Tetrahedron')
axes.plot(four_atom_roadsign[0],four_atom_roadsign[2],'d-',markerfacecolor = (1,1,1,0.9),label='Roadsign')
axes.plot(two_diatomics[0],two_diatomics[2],'2-', markersize = 15,markerfacecolor = (1,1,1,0.9),label='2 Diatomics')

axes.set_xlabel('Number of electrons')
axes.set_ylabel(r'Huckel Energy/ $\beta$')
axes.set_title('Plot of adjusted Huckel energies against number of electrons in system ')

axes.legend(loc=0)

#Saving figure to directory
fig_adjusted.savefig('Adjusted Huckel Energies')






fig_normal = plt.figure(figsize = (8,5),dpi = 1000)

axes = fig_normal.add_axes([0.1,0.1,0.8,0.8])

axes.plot(four_atom_cyclic[0],four_atom_cyclic[1],'+-',markersize = 10,label='Cyclic')
axes.plot(four_atom_linear[0],four_atom_linear[1],'o-',markersize = 5,markerfacecolor = 'white',label='Linear')
axes.plot(four_atom_rhombus[0],four_atom_rhombus[1],'^-',markerfacecolor=(1,1,1,0.9),label='Rhombus')
axes.plot(four_atom_tetrahedron[0],four_atom_tetrahedron[1],'s-',markerfacecolor = (1,1,1,0.9),label='Tetrahedron')
axes.plot(four_atom_roadsign[0],four_atom_roadsign[1],'d-',markerfacecolor = (1,1,1,0.9),label='Roadsign')
axes.plot(two_diatomics[0],two_diatomics[1],'2-', markersize = 15,markerfacecolor = (1,1,1,0.9),label='2 Diatomics')

axes.set_xlabel('Number of electrons')
axes.set_ylabel(r'Huckel Energy/ $\beta$')
axes.set_title('Plot of Huckel energies against number of electrons in system ')

axes.legend(loc=0)

#Saving figure to directory
fig_normal.savefig('Normal Huckel Energies')




fig_normal_three_atom = plt.figure(figsize = (8,5),dpi = 1000)

axes = fig_normal_three_atom.add_axes([0.1,0.1,0.8,0.8])

axes.plot(three_atom_cyclic[0],three_atom_cyclic[1],'+-',markersize = 10,label='Cyclic')
axes.plot(three_atom_linear[0],three_atom_linear[1],'o-',markersize = 5,markerfacecolor = 'white',label='Linear')

axes.set_ylabel(r'Huckel Energy/ $\beta$')
axes.set_title('Plot of Huckel energies against number of electrons in system ')

axes.legend(loc=0)

#Saving figure to directory
fig_normal_three_atom.savefig('Normal Huckel Three atom')



fig_adjusted_three_atom = plt.figure(figsize = (8,5),dpi = 1000)

axes = fig_adjusted_three_atom.add_axes([0.1,0.1,0.8,0.8])

axes.plot(three_atom_cyclic[0],three_atom_cyclic[2],'+-',markersize = 10,label='Cyclic')
axes.plot(three_atom_linear[0],three_atom_linear[2],'o-',markersize = 5,markerfacecolor = 'white',label='Linear')

axes.set_ylabel(r'Huckel Energy/ $\beta$')
axes.set_title('Plot of adjusted Huckel energies against number of electrons in system ')

axes.legend(loc=0)

#Saving figure to directory

fig_adjusted_three_atom.savefig('Adjusted Huckel Three atom')






#Generating the Pandas tables to read into Latex

data = three_atom_linear[0] + three_atom_linear[1]+three_atom_linear[2] + three_atom_cyclic[1]+three_atom_cyclic[2]

data = np.reshape(data,(5,6)).T

outside = ['','Linear Energies','Linear Energies','Cyclic Energies','Cyclic Energies']
inside = ['Electrons','Normal','Adjusted','Normal','Adjusted']
hier_index = list(zip(outside,inside))
hier_index = pd.MultiIndex.from_tuples(hier_index)

pd.options.mode.chained_assignment = None

df = pd.DataFrame(data,index=None,columns=hier_index)

df.loc[:,(slice(None),'Electrons')] = df.loc[:,(slice(None),'Electrons')].astype(int)

with open('three_atom_energies.tex','w') as tf:
    tf.write(df.to_latex())




data2 = four_atom_linear[0] + four_atom_linear[1]+four_atom_linear[2] + four_atom_cyclic[1] + four_atom_cyclic[2]

data2 = np.reshape(data2,(5,8)).T

outside = ['','Linear Energies','Linear Energies','Cyclic Energies','Cyclic Energies']
inside = ['Electrons','Normal','Adjusted','Normal','Adjusted']
hier_index = list(zip(outside,inside))
hier_index = pd.MultiIndex.from_tuples(hier_index)

pd.options.mode.chained_assignment = None

df = pd.DataFrame(data2,index=None,columns=hier_index)

df.loc[:,(slice(None),'Electrons')] = df.loc[:,(slice(None),'Electrons')].astype(int)

with open('four_atom_lin_cyc.tex','w') as tf:
    tf.write(df.to_latex(column_format='cccccc',index=False,multicolumn_format='c'))




data3 = four_atom_tetrahedron[0] + four_atom_tetrahedron[1]+four_atom_tetrahedron[2] + four_atom_rhombus[1] + four_atom_rhombus[2]

data3 = np.reshape(data3,(5,8)).T

outside = ['','Linear Energies','Linear Energies','Cyclic Energies','Cyclic Energies']
inside = ['Electrons','Normal','Adjusted','Normal','Adjusted']
hier_index = list(zip(outside,inside))
hier_index = pd.MultiIndex.from_tuples(hier_index)

pd.options.mode.chained_assignment = None

df = pd.DataFrame(data3,index=None,columns=hier_index)

df.loc[:,(slice(None),'Electrons')] = df.loc[:,(slice(None),'Electrons')].astype(int)

with open('four_atom_tet_rhomb.tex','w') as tf:
    tf.write(df.to_latex(column_format='cccccc',index=False,multicolumn_format='c'))



data4 = four_atom_roadsign[0] + four_atom_roadsign[1]+four_atom_roadsign[2] + two_diatomics[1] + two_diatomics[2]

data4 = np.reshape(data4,(5,8)).T

outside = ['','Roadsign Energies','Roadsign Energies','2 Diatomics Energies','2 Diatomics Energies']
inside = ['Electrons','Normal','Adjusted','Normal','Adjusted']
hier_index = list(zip(outside,inside))
hier_index = pd.MultiIndex.from_tuples(hier_index)

pd.options.mode.chained_assignment = None

df = pd.DataFrame(data4,index=None,columns=hier_index)

df.loc[:,(slice(None),'Electrons')] = df.loc[:,(slice(None),'Electrons')].astype(int)

with open('four_atom_roadsign_diatomics.tex','w') as tf:
    tf.write(df.to_latex(column_format='cccccc',index=False,multicolumn_format='c'))



#Generate the momet table and write it into Latex, finally...
momentdata = four_atom_linear[4]+four_atom_cyclic[4]+four_atom_tetrahedron[4]+four_atom_rhombus[4]+four_atom_roadsign[4]+two_diatomics[4]
momentdata = np.reshape(momentdata,(6,4))
df = pd.DataFrame(momentdata,['Linear', 'Cyclic', 'Tetrahedron', 'Rhombus', 'Roadsign', '2 Diatomics'],['0th moment', '1st moment','2nd moment','3rd moment'])
df.index.name = 'Structure'
df = df.astype(int)

with open('moments.tex','w') as tf:
    tf.write(df.to_latex())



#use usepackage{booktab} to read in the tables into Latex