import numpy as np
from sys import stdout
from openmm.app import *
from openmm import *
from openmm.unit import *
""" program to calculate the atom-center of mass distance and Rg from a pdb file. OA 20230122"""
myPdb = PDBFile('BoltornG4Dry1K.pdb')
#Time of the step in 1k run
time_1K = np.linspace(4*1000/20, 4*1000, 20)

positions_1 = myPdb.getPositions(frame=0)  #get the positions for the frame in nanometer
X_1 = [float(coord[0]/nanometer) for coord in positions_1]
Y_1 = [float(coord[1]/nanometer) for coord in positions_1]
Z_1 = [float(coord[2]/nanometer) for coord in positions_1]
noAtoms = len(X_1)#Assign of the number of atoms
#Assign of atoms
myAtoms = [atom for atom in myPdb.topology.atoms()]
# The elements are returned as <Element carbon>
myElements = [str(atom.element)[9:-1] for atom in myAtoms] #Assign of elemets names
myMasses =[]# atom masses
newX_1 = [] #The coordinates of the atoms with respect to center of mass
newY_1 = []
newZ_1 = []
for i in range(noAtoms):
    if myElements[i] == 'hydrogen':
        myMasses.append(1.00709)
    if myElements[i] == 'carbon':
        myMasses.append(12.011)
    if myElements[i] == 'oxygen':
        myMasses.append(15.999)
    if X_1[i]==X_1[0]:
        newX_1.append(0.0)
    else:
        newX_1.append(X_1[i]-X_1[0])
    if Y_1[i]==Y_1[0]:
        newY_1.append(0.0)
    else:
        newY_1.append(Y_1[i]-Y_1[0])
    if Z_1[i]==Z_1[0]:
        newZ_1.append(0.0)
    else:
        newZ_1.append(Z_1[i]-Z_1[0])
#Total number of mass
totalMass = sum(myMasses)
#(atom-center of mass) distances
myDistances_1 = np.sqrt([(newX_1[i]**2) + (newY_1[i]**2)+ (newZ_1[i]**2) for i in range(noAtoms)])
mySquaredDistances_1 = [] 
for x in myDistances_1:
    mySquaredDistances_1.append(x**2)
#The square of (atom-center of mass) * atom mass
myMassWeightedSqDists_1 = [mySquaredDistances_1[i] * myMasses[i] for i in range(noAtoms)]
mySquare_1Sum = sum(myMassWeightedSqDists_1)
# Radius of gyration of the structure
radiusOfGyration_1 = np.sqrt(mySquare_1Sum/totalMass)
#The average of(atom-center of mass) distance of the structure
print('Average distance_1 = ', sum(myDistances_1)/(noAtoms), 'nm')   
print('Radius of Gyration_1 is:', radiusOfGyration_1, 'nm')
#(atom-center of mass) distance
print('distance_1(atom-center of mass) = ', myDistances_1)

positions_2 = myPdb.getPositions(frame=1)  #get the positions for the frame n nanometer
X_2 = [float(coord[0]/nanometer) for coord in positions_2]
Y_2 = [float(coord[1]/nanometer) for coord in positions_2]
Z_2 = [float(coord[2]/nanometer) for coord in positions_2]
newX_2 = [] #The coordinates of the atoms with respect to center of mass
newY_2 = []
newZ_2 = []
for i in range(noAtoms):
    if X_2[i]==X_2[0]:
        newX_2.append(0.0)
    else:
        newX_2.append(X_2[i]-X_2[0])
    if Y_2[i]==Y_2[0]:
        newY_2.append(0.0)
    else:
        newY_2.append(Y_2[i]-Y_2[0])
    if Z_2[i]==Z_2[0]:
        newZ_2.append(0.0)
    else:
        newZ_2.append(Z_2[i]-Z_2[0])
#(atom-center of mass) distances
myDistances_2 = np.sqrt([(newX_2[i]**2) + (newY_2[i]**2)+ (newZ_2[i]**2) for i in range(noAtoms)])
mySquaredDistances_2 = [] 
for x in myDistances_2:
    mySquaredDistances_2.append(x**2)
#The square of (atom-center of mass) * atom mass
myMassWeightedSqDists_2 = [mySquaredDistances_2[i] * myMasses[i] for i in range(noAtoms)]
mySquare_2Sum = sum(myMassWeightedSqDists_2)
# Radius of gyration of the structure
radiusOfGyration_2 = np.sqrt(mySquare_2Sum/totalMass)
#The average of(atom-center of mass) distance of the structure
print('Average distance_2 = ', sum(myDistances_2)/(noAtoms), 'nm')   
print('Radius of Gyration_2 is:', radiusOfGyration_2, 'nm')
#(atom-center of mass) distances
print('distance_2(atom-center of mass) = ',myDistances_2)

positions_3 = myPdb.getPositions(frame=2) #get the positions for the frame in nanometer
X_3 = [float(coord[0]/nanometer) for coord in positions_3]
Y_3 = [float(coord[1]/nanometer) for coord in positions_3]
Z_3 = [float(coord[2]/nanometer) for coord in positions_3]
newX_3 = [] #The coordinates of the atoms with respect to center of mass
newY_3 = []
newZ_3 = []
for i in range(noAtoms):
    if X_3[i]==X_3[0]:
        newX_3.append(0.0)
    else:
        newX_3.append(X_3[i]-X_3[0])
    if Y_3[i]==Y_3[0]:
        newY_3.append(0.0)
    else:
        newY_3.append(Y_3[i]-Y_3[0])
    if Z_3[i]==Z_3[0]:
        newZ_3.append(0.0)
    else:
        newZ_3.append(Z_3[i]-Z_3[0])
#(atom-center of mass) distances
myDistances_3 = np.sqrt([(newX_3[i]**2) + (newY_3[i]**2)+ (newZ_3[i]**2) for i in range(noAtoms)])
mySquaredDistances_3 = [] 
for x in myDistances_3:
    mySquaredDistances_3.append(x**2)
#The square of (atom-center of mass) * atom mass
myMassWeightedSqDists_3 = [mySquaredDistances_3[i] * myMasses[i] for i in range(noAtoms)]
mySquare_3Sum = sum(myMassWeightedSqDists_3)
# Radius of gyration of the structure
radiusOfGyration_3 = np.sqrt(mySquare_3Sum/totalMass)
#The average of(atom-center of mass) distance of the structure
print('Average distance_3 = ', sum(myDistances_3)/(noAtoms), 'nm')   
print('Radius of Gyration_3 is:', radiusOfGyration_3, 'nm')
#(atom-center of mass) distances
print('distance_3 = ', myDistances_3)


positions_4 = myPdb.getPositions(frame=3) #get the positions for the frame in nanometer
X_4 = [float(coord[0]/nanometer) for coord in positions_4]
Y_4 = [float(coord[1]/nanometer) for coord in positions_4]
Z_4 = [float(coord[2]/nanometer) for coord in positions_4]
newX_4 = [] #The coordinates of the atoms with respect to center of mass
newY_4 = []
newZ_4 = []
for i in range(noAtoms):
    if X_4[i]==X_4[0]:
        newX_4.append(0.0)
    else:
        newX_4.append(X_4[i]-X_4[0])
    if Y_4[i]==Y_4[0]:
        newY_4.append(0.0)
    else:
        newY_4.append(Y_4[i]-Y_4[0])
    if Z_4[i]==Z_4[0]:
        newZ_4.append(0.0)
    else:
        newZ_4.append(Z_4[i]-Z_4[0])
#(atom-center of mass) distances
myDistances_4 = np.sqrt([(newX_4[i]**2) + (newY_4[i]**2)+ (newZ_4[i]**2) for i in range(noAtoms)])
mySquaredDistances_4 = [] 
for x in myDistances_4:
    mySquaredDistances_4.append(x**2)
#The square of (atom-center of mass) * atom mass
myMassWeightedSqDists_4 = [mySquaredDistances_4[i] * myMasses[i] for i in range(noAtoms)]
mySquare_4Sum = sum(myMassWeightedSqDists_4)
# Radius of gyration of the structure
radiusOfGyration_4 = np.sqrt(mySquare_4Sum/totalMass)
#The average of(atom-center of mass) distance of the structure
print('Average distance_4 = ', sum(myDistances_4)/(noAtoms), 'nm')   
print('Radius of Gyration_4 is:', radiusOfGyration_4, 'nm')
#(atom-center of mass) distances
print('distance_4 = ', myDistances_4)


positions_5 = myPdb.getPositions(frame=4) #get the positions for the frame in nanometer
X_5 = [float(coord[0]/nanometer) for coord in positions_5]
Y_5 = [float(coord[1]/nanometer) for coord in positions_5]
Z_5 = [float(coord[2]/nanometer) for coord in positions_5]
newX_5 = [] #The coordinates of the atoms with respect to center of mass
newY_5 = []
newZ_5 = []
for i in range(noAtoms):
    if X_5[i]==X_5[0]:
        newX_5.append(0.0)
    else:
        newX_5.append(X_5[i]-X_5[0])
    if Y_5[i]==Y_5[0]:
        newY_5.append(0.0)
    else:
        newY_5.append(Y_5[i]-Y_5[0])
    if Z_5[i]==Z_5[0]:
        newZ_5.append(0.0)
    else:
        newZ_5.append(Z_5[i]-Z_5[0])
#(atom-center of mass) distances
myDistances_5 = np.sqrt([(newX_5[i]**2) + (newY_5[i]**2)+ (newZ_5[i]**2) for i in range(noAtoms)])
mySquaredDistances_5 = [] 
for x in myDistances_5:
    mySquaredDistances_5.append(x**2)
#The square of (atom-center of mass) * atom mass
myMassWeightedSqDists_5 = [mySquaredDistances_5[i] * myMasses[i] for i in range(noAtoms)]
mySquare_5Sum = sum(myMassWeightedSqDists_5)
# Radius of gyration of the structure
radiusOfGyration_5 = np.sqrt(mySquare_5Sum/totalMass)
#The average of(atom-center of mass) distance of the structure
print('Average distance_5 = ', sum(myDistances_5)/(noAtoms), 'nm')   
print('Radius of Gyration_5 is:', radiusOfGyration_5, 'nm')
#(atom-center of mass) distances
print('distance_5 = ', myDistances_5)


positions_6 = myPdb.getPositions(frame=5) #get the positions for the frame in nanometer
X_6 = [float(coord[0]/nanometer) for coord in positions_6]
Y_6 = [float(coord[1]/nanometer) for coord in positions_6]
Z_6 = [float(coord[2]/nanometer) for coord in positions_6]
newX_6 = [] #The coordinates of the atoms with respect to center of mass
newY_6 = []
newZ_6 = []
for i in range(noAtoms):
    if X_6[i]==X_6[0]:
        newX_6.append(0.0)
    else:
        newX_6.append(X_6[i]-X_6[0])
    if Y_6[i]==Y_6[0]:
        newY_6.append(0.0)
    else:
        newY_6.append(Y_6[i]-Y_6[0])
    if Z_6[i]==Z_6[0]:
        newZ_6.append(0.0)
    else:
        newZ_6.append(Z_6[i]-Z_6[0])
#(atom-center of mass) distances
myDistances_6 = np.sqrt([(newX_6[i]**2) + (newY_6[i]**2)+ (newZ_6[i]**2) for i in range(noAtoms)])
mySquaredDistances_6 = [] 
for x in myDistances_6:
    mySquaredDistances_6.append(x**2)
#The square of (atom-center of mass) * atom mass
myMassWeightedSqDists_6 = [mySquaredDistances_6[i] * myMasses[i] for i in range(noAtoms)]
mySquare_6Sum = sum(myMassWeightedSqDists_6)
# Radius of gyration of the structure
radiusOfGyration_6 = np.sqrt(mySquare_6Sum/totalMass)
#The average of(atom-center of mass) distance of the structure
print('Average distance_6 = ', sum(myDistances_6)/(noAtoms), 'nm')   
print('Radius of Gyration_6 is:', radiusOfGyration_6, 'nm')
#(atom-center of mass) distances
print('distance_6 = ', myDistances_6)


positions_7 = myPdb.getPositions(frame=6) #get the positions for the frame in nanometer
X_7 = [float(coord[0]/nanometer) for coord in positions_7]
Y_7 = [float(coord[1]/nanometer) for coord in positions_7]
Z_7 = [float(coord[2]/nanometer) for coord in positions_7]
newX_7 = [] #The coordinates of the atoms with respect to center of mass
newY_7 = []
newZ_7 = []
for i in range(noAtoms):
    if X_7[i]==X_7[0]:
        newX_7.append(0.0)
    else:
        newX_7.append(X_7[i]-X_7[0])
    if Y_7[i]==Y_7[0]:
        newY_7.append(0.0)
    else:
        newY_7.append(Y_7[i]-Y_7[0])
    if Z_7[i]==Z_7[0]:
        newZ_7.append(0.0)
    else:
        newZ_7.append(Z_7[i]-Z_7[0])
#(atom-center of mass) distances
myDistances_7 = np.sqrt([(newX_7[i]**2) + (newY_7[i]**2)+ (newZ_7[i]**2) for i in range(noAtoms)])
mySquaredDistances_7 = [] 
for x in myDistances_7:
    mySquaredDistances_7.append(x**2)
#The square of (atom-center of mass) * atom mass
myMassWeightedSqDists_7 = [mySquaredDistances_7[i] * myMasses[i] for i in range(noAtoms)]
mySquare_7Sum = sum(myMassWeightedSqDists_7)
# Radius of gyration of the structure
radiusOfGyration_7 = np.sqrt(mySquare_7Sum/totalMass)
#The average of(atom-center of mass) distance of the structure
print('Average distance_7 = ', sum(myDistances_7)/(noAtoms), 'nm')   
print('Radius of Gyration_7 is:', radiusOfGyration_7, 'nm')
#(atom-center of mass) distances
print('distance_7 = ', myDistances_7)


positions_8 = myPdb.getPositions(frame=7) #get the positions for the frame in nanometer
X_8 = [float(coord[0]/nanometer) for coord in positions_8]
Y_8 = [float(coord[1]/nanometer) for coord in positions_8]
Z_8 = [float(coord[2]/nanometer) for coord in positions_8]
newX_8 = [] #The coordinates of the atoms with respect to center of mass
newY_8 = []
newZ_8 = []
for i in range(noAtoms):
    if X_8[i]==X_8[0]:
        newX_8.append(0.0)
    else:
        newX_8.append(X_8[i]-X_8[0])
    if Y_8[i]==Y_8[0]:
        newY_8.append(0.0)
    else:
        newY_8.append(Y_8[i]-Y_8[0])
    if Z_8[i]==Z_8[0]:
        newZ_8.append(0.0)
    else:
        newZ_8.append(Z_8[i]-Z_8[0])
#(atom-center of mass) distances
myDistances_8 = np.sqrt([(newX_8[i]**2) + (newY_8[i]**2)+ (newZ_8[i]**2) for i in range(noAtoms)])
mySquaredDistances_8 = [] 
for x in myDistances_8:
    mySquaredDistances_8.append(x**2)
#The square of (atom-center of mass) * atom mass
myMassWeightedSqDists_8 = [mySquaredDistances_8[i] * myMasses[i] for i in range(noAtoms)]
mySquare_8Sum = sum(myMassWeightedSqDists_8)
# Radius of gyration of the structure
radiusOfGyration_8 = np.sqrt(mySquare_8Sum/totalMass)
#The average of(atom-center of mass) distance of the structure
print('Average distance_8 = ', sum(myDistances_8)/(noAtoms), 'nm')   
print('Radius of Gyration_8 is:', radiusOfGyration_8, 'nm')
#(atom-center of mass) distances
print('distance_8 = ', myDistances_8)


positions_9 = myPdb.getPositions(frame=8) #get the positions for the frame in nanometer
X_9 = [float(coord[0]/nanometer) for coord in positions_9]
Y_9 = [float(coord[1]/nanometer) for coord in positions_9]
Z_9 = [float(coord[2]/nanometer) for coord in positions_9]
newX_9 = [] #The coordinates of the atoms with respect to center of mass
newY_9 = []
newZ_9 = []
for i in range(noAtoms):
    if X_9[i]==X_9[0]:
        newX_9.append(0.0)
    else:
        newX_9.append(X_9[i]-X_9[0])
    if Y_9[i]==Y_9[0]:
        newY_9.append(0.0)
    else:
        newY_9.append(Y_9[i]-Y_9[0])
    if Z_9[i]==Z_9[0]:
        newZ_9.append(0.0)
    else:
        newZ_9.append(Z_9[i]-Z_9[0])
#(atom-center of mass) distances
myDistances_9 = np.sqrt([(newX_9[i]**2) + (newY_9[i]**2)+ (newZ_9[i]**2) for i in range(noAtoms)])
mySquaredDistances_9 = [] 
for x in myDistances_9:
    mySquaredDistances_9.append(x**2)
#The square of (atom-center of mass) * atom mass
myMassWeightedSqDists_9 = [mySquaredDistances_9[i] * myMasses[i] for i in range(noAtoms)]
mySquare_9Sum = sum(myMassWeightedSqDists_9)
# Radius of gyration of the structure
radiusOfGyration_9 = np.sqrt(mySquare_9Sum/totalMass)
#The average of(atom-center of mass) distance of the structure
print('Average distance_9 = ', sum(myDistances_9)/(noAtoms), 'nm')   
print('Radius of Gyration_9 is:', radiusOfGyration_9, 'nm')
#(atom-center of mass) distances
print('distance_9 = ', myDistances_9)


positions_10 = myPdb.getPositions(frame=9) #get the positions for the frame in nanometer
X_10 = [float(coord[0]/nanometer) for coord in positions_10]
Y_10 = [float(coord[1]/nanometer) for coord in positions_10]
Z_10 = [float(coord[2]/nanometer) for coord in positions_10]
newX_10 = [] #The coordinates of the atoms with respect to center of mass
newY_10 = []
newZ_10 = []
for i in range(noAtoms):
    if X_10[i]==X_10[0]:
        newX_10.append(0.0)
    else:
        newX_10.append(X_10[i]-X_10[0])
    if Y_10[i]==Y_10[0]:
        newY_10.append(0.0)
    else:
        newY_10.append(Y_10[i]-Y_10[0])
    if Z_10[i]==Z_10[0]:
        newZ_10.append(0.0)
    else:
        newZ_10.append(Z_10[i]-Z_10[0])
#(atom-center of mass) distances
myDistances_10 = np.sqrt([(newX_10[i]**2) + (newY_10[i]**2)+ (newZ_10[i]**2) for i in range(noAtoms)])
mySquaredDistances_10 = [] 
for x in myDistances_10:
    mySquaredDistances_10.append(x**2)
#The square of (atom-center of mass) * atom mass
myMassWeightedSqDists_10 = [mySquaredDistances_10[i] * myMasses[i] for i in range(noAtoms)]
mySquare_10Sum = sum(myMassWeightedSqDists_10)
# Radius of gyration of the structure
radiusOfGyration_10 = np.sqrt(mySquare_10Sum/totalMass)
#The average of(atom-center of mass) distance of the structure
print('Average distance_10 = ', sum(myDistances_10)/(noAtoms), 'nm')   
print('Radius of Gyration_10 is:', radiusOfGyration_10, 'nm')
#(atom-center of mass) distances
print('distance_10 = ', myDistances_10)


positions_11 = myPdb.getPositions(frame=10) #get the positions for the frame in nanometer
X_11 = [float(coord[0]/nanometer) for coord in positions_11]
Y_11 = [float(coord[1]/nanometer) for coord in positions_11]
Z_11 = [float(coord[2]/nanometer) for coord in positions_11]
newX_11 = [] #The coordinates of the atoms with respect to center of mass
newY_11 = []
newZ_11 = []
for i in range(noAtoms):
    if X_11[i]==X_11[0]:
        newX_11.append(0.0)
    else:
        newX_11.append(X_11[i]-X_11[0])
    if Y_11[i]==Y_11[0]:
        newY_11.append(0.0)
    else:
        newY_11.append(Y_11[i]-Y_11[0])
    if Z_11[i]==Z_11[0]:
        newZ_11.append(0.0)
    else:
        newZ_11.append(Z_11[i]-Z_11[0])
#(atom-center of mass) distances
myDistances_11 = np.sqrt([(newX_11[i]**2) + (newY_11[i]**2)+ (newZ_11[i]**2) for i in range(noAtoms)])
mySquaredDistances_11 = [] 
for x in myDistances_11:
    mySquaredDistances_11.append(x**2)
#The square of (atom-center of mass) * atom mass
myMassWeightedSqDists_11 = [mySquaredDistances_11[i] * myMasses[i] for i in range(noAtoms)]
mySquare_11Sum = sum(myMassWeightedSqDists_11)
# Radius of gyration of the structure
radiusOfGyration_11 = np.sqrt(mySquare_11Sum/totalMass)
#The average of(atom-center of mass) distance of the structure
print('Average distance_11 = ', sum(myDistances_11)/(noAtoms), 'nm')   
print('Radius of Gyration_11 is:', radiusOfGyration_11, 'nm')
#(atom-center of mass) distances
print('distance_11 = ', myDistances_11)

positions_12 = myPdb.getPositions(frame=11) #get the positions for the frame in nanometer
X_12 = [float(coord[0]/nanometer) for coord in positions_12]
Y_12 = [float(coord[1]/nanometer) for coord in positions_12]
Z_12 = [float(coord[2]/nanometer) for coord in positions_12]
newX_12 = [] #The coordinates of the atoms with respect to center of mass
newY_12 = []
newZ_12 = []
for i in range(noAtoms):
    if X_12[i]==X_12[0]:
        newX_12.append(0.0)
    else:
        newX_12.append(X_12[i]-X_12[0])
    if Y_12[i]==Y_12[0]:
        newY_12.append(0.0)
    else:
        newY_12.append(Y_12[i]-Y_12[0])
    if Z_12[i]==Z_12[0]:
        newZ_12.append(0.0)
    else:
        newZ_12.append(Z_12[i]-Z_12[0])
#(atom-center of mass) distances
myDistances_12 = np.sqrt([(newX_12[i]**2) + (newY_12[i]**2)+ (newZ_12[i]**2) for i in range(noAtoms)])
mySquaredDistances_12 = [] 
for x in myDistances_12:
    mySquaredDistances_12.append(x**2)
#The square of (atom-center of mass) * atom mass
myMassWeightedSqDists_12 = [mySquaredDistances_12[i] * myMasses[i] for i in range(noAtoms)]
mySquare_12Sum = sum(myMassWeightedSqDists_12)
# Radius of gyration of the structure
radiusOfGyration_12 = np.sqrt(mySquare_12Sum/totalMass)
#The average of(atom-center of mass) distance of the structure
print('Average distance_12 = ', sum(myDistances_12)/(noAtoms), 'nm')   
print('Radius of Gyration_12 is:', radiusOfGyration_12, 'nm')
#(atom-center of mass) distances
print('distance_12 = ', myDistances_12)


positions_13 = myPdb.getPositions(frame=12) #get the positions for the frame in nanometer
X_13 = [float(coord[0]/nanometer) for coord in positions_13]
Y_13 = [float(coord[1]/nanometer) for coord in positions_13]
Z_13 = [float(coord[2]/nanometer) for coord in positions_13]
newX_13 = [] #The coordinates of the atoms with respect to center of mass
newY_13 = []
newZ_13 = []
for i in range(noAtoms):
    if X_13[i]==X_13[0]:
        newX_13.append(0.0)
    else:
        newX_13.append(X_13[i]-X_13[0])
    if Y_13[i]==Y_13[0]:
        newY_13.append(0.0)
    else:
        newY_13.append(Y_13[i]-Y_13[0])
    if Z_13[i]==Z_13[0]:
        newZ_13.append(0.0)
    else:
        newZ_13.append(Z_13[i]-Z_13[0])
#(atom-center of mass) distances
myDistances_13 = np.sqrt([(newX_13[i]**2) + (newY_13[i]**2)+ (newZ_13[i]**2) for i in range(noAtoms)])
mySquaredDistances_13 = [] 
for x in myDistances_13:
    mySquaredDistances_13.append(x**2)
#The square of (atom-center of mass) * atom mass
myMassWeightedSqDists_13 = [mySquaredDistances_13[i] * myMasses[i] for i in range(noAtoms)]
mySquare_13Sum = sum(myMassWeightedSqDists_13)
# Radius of gyration of the structure
radiusOfGyration_13 = np.sqrt(mySquare_13Sum/totalMass)
#The average of(atom-center of mass) distance of the structure
print('Average distance_13 = ', sum(myDistances_13)/(noAtoms), 'nm')   
print('Radius of Gyration_13 is:', radiusOfGyration_13, 'nm')
#(atom-center of mass) distances
print('distance_13 = ', myDistances_13)


positions_14 = myPdb.getPositions(frame=13) #get the positions for the frame in nanometer
X_14 = [float(coord[0]/nanometer) for coord in positions_14]
Y_14 = [float(coord[1]/nanometer) for coord in positions_14]
Z_14 = [float(coord[2]/nanometer) for coord in positions_14]
newX_14 = [] #The coordinates of the atoms with respect to center of mass
newY_14 = []
newZ_14 = []
for i in range(noAtoms):
    if X_14[i]==X_14[0]:
        newX_14.append(0.0)
    else:
        newX_14.append(X_14[i]-X_14[0])
    if Y_14[i]==Y_14[0]:
        newY_14.append(0.0)
    else:
        newY_14.append(Y_14[i]-Y_14[0])
    if Z_14[i]==Z_14[0]:
        newZ_14.append(0.0)
    else:
        newZ_14.append(Z_14[i]-Z_14[0])
#(atom-center of mass) distances
myDistances_14 = np.sqrt([(newX_14[i]**2) + (newY_14[i]**2)+ (newZ_14[i]**2) for i in range(noAtoms)])
mySquaredDistances_14 = [] 
for x in myDistances_14:
    mySquaredDistances_14.append(x**2)
#The square of (atom-center of mass) * atom mass
myMassWeightedSqDists_14 = [mySquaredDistances_14[i] * myMasses[i] for i in range(noAtoms)]
mySquare_14Sum = sum(myMassWeightedSqDists_14)
# Radius of gyration of the structure
radiusOfGyration_14 = np.sqrt(mySquare_14Sum/totalMass)
#The average of(atom-center of mass) distance of the structure
print('Average distance_14 = ', sum(myDistances_14)/(noAtoms), 'nm')   
print('Radius of Gyration_14 is:', radiusOfGyration_14, 'nm')
#(atom-center of mass) distances
print('distance_14 = ', myDistances_14)


positions_15 = myPdb.getPositions(frame=14) #get the positions for the frame in nanometer
X_15 = [float(coord[0]/nanometer) for coord in positions_15]
Y_15 = [float(coord[1]/nanometer) for coord in positions_15]
Z_15 = [float(coord[2]/nanometer) for coord in positions_15]
newX_15 = [] #The coordinates of the atoms with respect to center of mass
newY_15 = []
newZ_15 = []
for i in range(noAtoms):
    if X_15[i]==X_15[0]:
        newX_15.append(0.0)
    else:
        newX_15.append(X_15[i]-X_15[0])
    if Y_15[i]==Y_15[0]:
        newY_15.append(0.0)
    else:
        newY_15.append(Y_15[i]-Y_15[0])
    if Z_15[i]==Z_15[0]:
        newZ_15.append(0.0)
    else:
        newZ_15.append(Z_15[i]-Z_15[0])
#(atom-center of mass) distances
myDistances_15 = np.sqrt([(newX_15[i]**2) + (newY_15[i]**2)+ (newZ_15[i]**2) for i in range(noAtoms)])
mySquaredDistances_15 = [] 
for x in myDistances_15:
    mySquaredDistances_15.append(x**2)
#The square of (atom-center of mass) * atom mass
myMassWeightedSqDists_15 = [mySquaredDistances_15[i] * myMasses[i] for i in range(noAtoms)]
mySquare_15Sum = sum(myMassWeightedSqDists_15)
# Radius of gyration of the structure
radiusOfGyration_15 = np.sqrt(mySquare_15Sum/totalMass)
#The average of(atom-center of mass) distance of the structure
print('Average distance_15 = ', sum(myDistances_15)/(noAtoms), 'nm')   
print('Radius of Gyration_15 is:', radiusOfGyration_15, 'nm')
#(atom-center of mass) distances
print('distance_15 = ', myDistances_15)


positions_16 = myPdb.getPositions(frame=15) #get the positions for the frame in nanometer
X_16 = [float(coord[0]/nanometer) for coord in positions_16]
Y_16 = [float(coord[1]/nanometer) for coord in positions_16]
Z_16 = [float(coord[2]/nanometer) for coord in positions_16]
newX_16 = [] #The coordinates of the atoms with respect to center of mass
newY_16 = []
newZ_16 = []
for i in range(noAtoms):
    if X_16[i]==X_16[0]:
        newX_16.append(0.0)
    else:
        newX_16.append(X_16[i]-X_16[0])
    if Y_16[i]==Y_16[0]:
        newY_16.append(0.0)
    else:
        newY_16.append(Y_16[i]-Y_16[0])
    if Z_16[i]==Z_16[0]:
        newZ_16.append(0.0)
    else:
        newZ_16.append(Z_16[i]-Z_16[0])
#(atom-center of mass) distances
myDistances_16 = np.sqrt([(newX_16[i]**2) + (newY_16[i]**2)+ (newZ_16[i]**2) for i in range(noAtoms)])
mySquaredDistances_16 = [] 
for x in myDistances_16:
    mySquaredDistances_16.append(x**2)
#The square of (atom-center of mass) * atom mass
myMassWeightedSqDists_16 = [mySquaredDistances_16[i] * myMasses[i] for i in range(noAtoms)]
mySquare_16Sum = sum(myMassWeightedSqDists_16)
# Radius of gyration of the structure
radiusOfGyration_16 = np.sqrt(mySquare_16Sum/totalMass)
#The average of(atom-center of mass) distance of the structure
print('Average distance_16 = ', sum(myDistances_16)/(noAtoms), 'nm')   
print('Radius of Gyration_16 is:', radiusOfGyration_16, 'nm')
#(atom-center of mass) distances
print('distance_16 = ', myDistances_16)


positions_17 = myPdb.getPositions(frame=16) #get the positions for the frame in nanometer
X_17 = [float(coord[0]/nanometer) for coord in positions_17]
Y_17 = [float(coord[1]/nanometer) for coord in positions_17]
Z_17 = [float(coord[2]/nanometer) for coord in positions_17]
newX_17 = [] #The coordinates of the atoms with respect to center of mass
newY_17 = []
newZ_17 = []
for i in range(noAtoms):
    if X_17[i]==X_17[0]:
        newX_17.append(0.0)
    else:
        newX_17.append(X_17[i]-X_17[0])
    if Y_17[i]==Y_17[0]:
        newY_17.append(0.0)
    else:
        newY_17.append(Y_17[i]-Y_17[0])
    if Z_17[i]==Z_17[0]:
        newZ_17.append(0.0)
    else:
        newZ_17.append(Z_17[i]-Z_17[0])
#(atom-center of mass) distances
myDistances_17 = np.sqrt([(newX_17[i]**2) + (newY_17[i]**2)+ (newZ_17[i]**2) for i in range(noAtoms)])
mySquaredDistances_17 = [] 
for x in myDistances_17:
    mySquaredDistances_17.append(x**2)
#The square of (atom-center of mass) * atom mass
myMassWeightedSqDists_17 = [mySquaredDistances_17[i] * myMasses[i] for i in range(noAtoms)]
mySquare_17Sum = sum(myMassWeightedSqDists_17)
# Radius of gyration of the structure
radiusOfGyration_17 = np.sqrt(mySquare_17Sum/totalMass)
#The average of(atom-center of mass) distance of the structure
print('Average distance_17 = ', sum(myDistances_17)/(noAtoms), 'nm')   
print('Radius of Gyration_17 is:', radiusOfGyration_17, 'nm')
#(atom-center of mass) distances
print('distance_17 = ', myDistances_17)


positions_18 = myPdb.getPositions(frame=17) #get the positions for the frame in nanometer
X_18 = [float(coord[0]/nanometer) for coord in positions_18]
Y_18 = [float(coord[1]/nanometer) for coord in positions_18]
Z_18 = [float(coord[2]/nanometer) for coord in positions_18]
newX_18 = [] #The coordinates of the atoms with respect to center of mass
newY_18 = []
newZ_18 = []
for i in range(noAtoms):
    if X_18[i]==X_18[0]:
        newX_18.append(0.0)
    else:
        newX_18.append(X_18[i]-X_18[0])
    if Y_18[i]==Y_18[0]:
        newY_18.append(0.0)
    else:
        newY_18.append(Y_18[i]-Y_18[0])
    if Z_18[i]==Z_18[0]:
        newZ_18.append(0.0)
    else:
        newZ_18.append(Z_18[i]-Z_18[0])
#(atom-center of mass) distances
myDistances_18 = np.sqrt([(newX_18[i]**2) + (newY_18[i]**2)+ (newZ_18[i]**2) for i in range(noAtoms)])
mySquaredDistances_18 = [] 
for x in myDistances_18:
    mySquaredDistances_18.append(x**2)
#The square of (atom-center of mass) * atom mass
myMassWeightedSqDists_18 = [mySquaredDistances_18[i] * myMasses[i] for i in range(noAtoms)]
mySquare_18Sum = sum(myMassWeightedSqDists_18)
# Radius of gyration of the structure
radiusOfGyration_18 = np.sqrt(mySquare_18Sum/totalMass)
#The average of(atom-center of mass) distance of the structure
print('Average distance_18 = ', sum(myDistances_18)/(noAtoms), 'nm')   
print('Radius of Gyration_18 is:', radiusOfGyration_18, 'nm')
#(atom-center of mass) distances
print('distance_18 = ', myDistances_18)


positions_19 = myPdb.getPositions(frame=18) #get the positions for the frame in nanometer
X_19 = [float(coord[0]/nanometer) for coord in positions_19]
Y_19 = [float(coord[1]/nanometer) for coord in positions_19]
Z_19 = [float(coord[2]/nanometer) for coord in positions_19]
newX_19 = [] #The coordinates of the atoms with respect to center of mass
newY_19 = []
newZ_19 = []
for i in range(noAtoms):
    if X_19[i]==X_19[0]:
        newX_19.append(0.0)
    else:
        newX_19.append(X_19[i]-X_19[0])
    if Y_19[i]==Y_19[0]:
        newY_19.append(0.0)
    else:
        newY_19.append(Y_19[i]-Y_19[0])
    if Z_19[i]==Z_19[0]:
        newZ_19.append(0.0)
    else:
        newZ_19.append(Z_19[i]-Z_19[0])
#(atom-center of mass) distances
myDistances_19 = np.sqrt([(newX_19[i]**2) + (newY_19[i]**2)+ (newZ_19[i]**2) for i in range(noAtoms)])
mySquaredDistances_19 = [] 
for x in myDistances_19:
    mySquaredDistances_19.append(x**2)
#The square of (atom-center of mass) * atom mass
myMassWeightedSqDists_19 = [mySquaredDistances_19[i] * myMasses[i] for i in range(noAtoms)]
mySquare_19Sum = sum(myMassWeightedSqDists_19)
# Radius of gyration of the structure
radiusOfGyration_19 = np.sqrt(mySquare_19Sum/totalMass)
#The average of(atom-center of mass) distance of the structure
print('Average distance_19 = ', sum(myDistances_19)/(noAtoms), 'nm')   
print('Radius of Gyration_19 is:', radiusOfGyration_19, 'nm')
#(atom-center of mass) distances
print('distance_19 = ', myDistances_19)


positions_20 = myPdb.getPositions(frame=19) #get the positions for the frame in nanometer
X_20 = [float(coord[0]/nanometer) for coord in positions_20]
Y_20 = [float(coord[1]/nanometer) for coord in positions_20]
Z_20 = [float(coord[2]/nanometer) for coord in positions_20]
newX_20 = [] #The coordinates of the atoms with respect to center of mass
newY_20 = []
newZ_20 = []
for i in range(noAtoms):
    if X_20[i]==X_20[0]:
        newX_20.append(0.0)
    else:
        newX_20.append(X_20[i]-X_20[0])
    if Y_20[i]==Y_20[0]:
        newY_20.append(0.0)
    else:
        newY_20.append(Y_20[i]-Y_20[0])
    if Z_20[i]==Z_20[0]:
        newZ_20.append(0.0)
    else:
        newZ_20.append(Z_20[i]-Z_20[0])
#(atom-center of mass) distances
myDistances_20 = np.sqrt([(newX_20[i]**2) + (newY_20[i]**2)+ (newZ_20[i]**2) for i in range(noAtoms)])
mySquaredDistances_20 = [] 
for x in myDistances_20:
    mySquaredDistances_20.append(x**2)
#The square of (atom-center of mass) * atom mass
myMassWeightedSqDists_20 = [mySquaredDistances_20[i] * myMasses[i] for i in range(noAtoms)]
mySquare_20Sum = sum(myMassWeightedSqDists_20)
# Radius of gyration of the structure
radiusOfGyration_20 = np.sqrt(mySquare_20Sum/totalMass)
#The average of(atom-center of mass) distance of the structure
print('Average distance_20 = ', sum(myDistances_20)/(noAtoms), 'nm')   
print('Radius of Gyration_20 is:', radiusOfGyration_20, 'nm')
#(atom-center of mass) distances
print('distance_20 = ', myDistances_20)

averageDistances_1K = [sum(myDistances_1)/noAtoms, sum(myDistances_2)/noAtoms, sum(myDistances_3)/noAtoms, sum(myDistances_4)/noAtoms, sum(myDistances_5)/noAtoms, sum(myDistances_6)/noAtoms, sum(myDistances_7)/noAtoms, sum(myDistances_8)/noAtoms, sum(myDistances_9)/noAtoms, sum(myDistances_10)/noAtoms, sum(myDistances_11)/noAtoms, sum(myDistances_12)/noAtoms, sum(myDistances_13)/noAtoms, sum(myDistances_14)/noAtoms, sum(myDistances_15)/noAtoms, sum(myDistances_16)/noAtoms, sum(myDistances_17)/noAtoms, sum(myDistances_18)/noAtoms, sum(myDistances_19)/noAtoms, sum(myDistances_20)/noAtoms]      

radiiOfGyration_1K = [radiusOfGyration_1, radiusOfGyration_2, radiusOfGyration_3, radiusOfGyration_4, radiusOfGyration_5, radiusOfGyration_6, radiusOfGyration_7, radiusOfGyration_8, radiusOfGyration_9, radiusOfGyration_10, radiusOfGyration_11, radiusOfGyration_12, radiusOfGyration_13, radiusOfGyration_14, radiusOfGyration_15, radiusOfGyration_16, radiusOfGyration_17, radiusOfGyration_18, radiusOfGyration_19, radiusOfGyration_20]

