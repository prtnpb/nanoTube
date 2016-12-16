#!/usr/bin/python

'''
Program: nanoTube.py
Author: Patrick Ryan Thomas
        prtnpb@mail.umkc.edu

Purpose: The purpose of this program is to take the smallest
         translationally invariant portion of a structure and
         through a series of translations and rotations, make a
         nanotube structure that meets periodic boundary conditions
         in the z-direction

Inputs: The program takes in a list of element names corresponding
        to a list of a list of atomic position. The program also
        requires the initial cell dimensions.

Calculations: The program will then calculate all of the plausible 
              radii and rotation angles by rules of polygon geometry
              and trigonometry.

              DETERMINATION OF TUBE RADII:

              Given the y-width of the cell as the length of one
              side of a polygon.  The total sum of the angles of the
              polygon is pi*(n-2) where n is the number of sides.

              Dividing this angle by the number of angles which is
              equal to n, we get each internal angle around the
              perimeter of the polygon. Further, dividing by two
              creates a right triangle. Therefore,

                 ( pi     )   r 
              tan| --(n-2)| = --
                 ( 2n     }   l/2

              DETERMINATION OF ROTATION ANGLES:

              For the angle of rotation, since we have a right
              triangle:

              angle_of_rotation=2(pi/2 - (pi/2n)*(n-2))

              Traditional rotation matrices are then used for
              rotations about the principle axis both in setting up
              the initial structures and in the duplication and 
              rotation of cells.

              Arbitrarily number of sides can be calculated for all
              the orientations from 15 to 80 sides. It then becomes
              the job of the investigator to choose the correct
              structure.

              CHIRALITY:

              For a given chirality, we want to ensure that z-
              periodicity is maintained. We accomplish this by
              limiting the height to one unit cell in the z-
              direction. We then define a total shift to be
              accomplished 

              FUNCTIONALIZATION:
              
              Have to consider width of atom or particle.
              Calculate random angle between 0 and 2*pi and
              random height.


Outputs: Skeleton file structure for use in OLCAO.
'''

import os
import sys
import numpy as np
import time
import copy

#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D


start=time.time()

def printLine():
  print "============================================================================"


# Rotation about x-axis clockwise by angle ang
def Rotx(ang):
  cosang=np.cos(ang)
  sinang=np.sin(ang)

  Rx=np.zeros((3,3),dtype='d')
  Rx[0,0]=1.0
  Rx[1,1]=cosang
  Rx[2,1]=sinang
  Rx[1,2]=-1.0*sinang
  Rx[2,2]=cosang

  return Rx

# Rotation about y-axis clockwise by angle ang
def Roty(ang):
  cosang=np.cos(ang)
  sinang=np.sin(ang)

  Ry=np.zeros((3,3),dtype='d')
  Ry[0,0]=cosang
  Ry[2,0]=-1.0*sinang
  Ry[1,1]=1.0
  Ry[0,2]=sinang
  Ry[2,2]=cosang

  return Ry

# Rotation about z-axis clockwise by angle ang
def Rotz(ang):
  cosang=np.cos(ang)
  sinang=np.sin(ang)

  Rz=np.zeros((3,3),dtype='d')
  Rz[0,0]=cosang
  Rz[1,0]=sinang
  Rz[0,1]=-1.0*sinang
  Rz[1,1]=cosang
  Rz[2,2]=1.0

  return Rz

''' 
This next class will open an olcao.mi and get the resulting
cell vectors, angles, atom names, positions.


''' 
class olcaomiFile:
  def __init__(self,filePath):
    self.atomicList=[]
    self.positionList=[]
    self.cellVec=np.zeros(3,dtype='d')
    self.cellAng=np.zeros(3,dtype='d')
    self.fracorcart=''
    self.numAtoms=0
    self.fP=filePath

  def getFile(self):
    olcaomi=open(self.fP,'r')
    allLines=olcaomi.readlines()
    olcaomi.close()


    for i in range(len(allLines)):

      if len(self.atomicList)==self.numAtoms and self.numAtoms!=0:
        break
    
      if allLines[i].strip()=='cell':
        tempLine=(allLines[i+1].strip()).split()
        self.cellVec[0]=tempLine[0]
        self.cellVec[1]=tempLine[1]
        self.cellVec[2]=tempLine[2]
        self.cellAng[0]=tempLine[3]
        self.cellAng[1]=tempLine[4]
        self.cellAng[2]=tempLine[5]

        tempLine=(allLines[i+2].strip()).split()
        self.fracorcart=tempLine[0]
        self.numAtoms=int(tempLine[1])        

        if self.fracorcart=='cart':
          for j in range(i+3,i+self.numAtoms+3):
            tempLine=(allLines[j].strip()).split()
            self.atomicList.append(tempLine[0])     
            self.positionList.append([tempLine[1],tempLine[2],tempLine[3]])        
        else:
          for j in range(i+3,i+self.numAtoms+3):
            tempLine=(allLines[j].strip()).split()
            self.atomicList.append(tempLine[0])     
            self.positionList.append([float(tempLine[1])*self.cellVec[0],float(tempLine[2])*self.cellVec[1],float(tempLine[3])*self.cellVec[2]])        

'''
  Rotation Descriptions:

  Each face of a unit cell has four different orientations,
  therefore we must consider twenty four posssible configurations.

  Similarly, two chiralities must also be taken into account.

  



 

  Rotation 1:
    Since all of the positions are in the first octant, i.e.
    (0<x,0<y,0<z). We simply take all positions in the provided
    translationally invariant structure and subtract y-width/2
    from all of the y values.

  Rotation 2:
   For this rotation, we do the following:

  Rotation 3:
    For this rotation, we do the following:

  Rotation 4:
    For this rotation, we do the following:

  Rotation 5:
    For this rotation, we do the following:

  Rotation 6:
    For this rotation, we do the following:
'''

class initStructure:
  '''
  CLASS VARIABLES
    initPositions: list of positions [[x_1,y_1,z_1],[x_2,y_2,z_2]...]
    corr_atoms: list of element names associated with initPositions
    width: vector containing the initial cell vectors.
    orientpos1: similar to initPositions, but with Rotation 1 
                (see above) applied
    width1: In this case, the x,y,z widths do not change because we aren't
            rotating the cell so this is the same as the initial.
 
    orientpos2: similar to initPositions, but with Rotation 2 applied.
    width2: x,y, and z widths after Rotation 2 is applied
    orientpos3: similar to initPositions, but with Rotation 3 applied.
    width3: x,y, and z widths after Rotation 3 is applied

  CLASS FUNCTIONS
    Initalization is accomplished by providing the inital cell widths.

    addAtom: This allows for adding in atoms and positions individually.
             Takes string of element type and list of x,y,z position.

    addEntireStructure: Allows for passing of list of element names and
                        list of list of atomic positions. 
                        See initPositions.
                        Takes list of element names and list of list
                        of atomic positions.

    setinitPos1: Performs the shifting described Rotation 1 on atoms in the
                 initial structure
      
    setinitPos2: Peforms the rotations and shifting described in
                 Rotation 2.

    setinitPos3: Performs the rotations and shifting described in
                 Rotation 3.

    setAllThreePos: Peforms all three functions of setinitPos1, 2, and 3.

    doRot1/2/3: 1) The loop determines the radii and the angle of rotation.
                2) Shift all x positions by the radii amount.
                3) Peform number sides -1 rotations about the z-axis.
    
    doAllRots: Will make all three functions

  '''
  
  def __init__(self):
    self.initPositions=[]           # List of initial positions
    self.corr_atoms=[]              # List of the atom symbols corresponding
                                    #  to initPositions.

    self.width=[]                   # Initial cell dimensions
    self.chir=0                     # Initial chirality set to zero.
    self.scrnum=1                   # Number of loops in scrolling
    self.scrdist=0                  # Distance to extend in scrolling
    self.func=0                     # Functionalization 0=No, 1=Yes
    self.funcLoc=0                  # Determines whether in, out, or both
    self.funcNumPart=0              # Number nanoparticles.
    self.funcPartWidth=[]           # Cell dimensions for particle
    self.funcinitPos=[]             # Atomic positions for particle
    self.funccorr_atoms=[]          # Atomic names for particle.
    self.funcElemNums=[]            # Number of times to replicate atoms.

    self.height=1                   # Number of times to replicate in
                                    #   z-direction.
    self.name=''                    # Name of Tube
    self.filePath=''                # File path of the olcao.mi file
    self.tuberadbounds=[]           # Number of replicated cells
    self.outfile=''                 # Name of the output file

    # List of positions of each atom when each face of the parallelpiped
    #  is facing the origin.
    self.orientpos1=np.zeros(1,dtype='d')
    self.orientpos1a=np.zeros(1,dtype='d')
    self.orientpos2=np.zeros(1,dtype='d')
    self.orientpos2a=np.zeros(1,dtype='d')
    self.orientpos3=np.zeros(1,dtype='d')
    self.orientpos3a=np.zeros(1,dtype='d')
    self.orientpos4=np.zeros(1,dtype='d')
    self.orientpos4a=np.zeros(1,dtype='d')
    self.orientpos5=np.zeros(1,dtype='d')
    self.orientpos5a=np.zeros(1,dtype='d')
    self.orientpos6=np.zeros(1,dtype='d')
    self.orientpos6a=np.zeros(1,dtype='d')

    # List of the cell widths in cartesian coordinates.
    self.width1=[]
    self.width1a=[]
    self.width2=[]
    self.width2a=[]
    self.width3=[]
    self.width3a=[]
    self.width4=[]
    self.width4a=[]
    self.width5=[]
    self.width5a=[]
    self.width6=[]
    self.width6a=[]

  def addAtom(self,atom,pos):
    self.initPositions.append(pos[:])
    self.corr_atoms.append(atom)

  def addEntireStructure(self,atomlist,poslist):
    self.initPositions=list(copy.deepcopy(poslist))
    self.corr_atoms=list(copy.deepcopy(atomlist))


  '''
  This functions parses the input control file 'nanoInCtrl' for
  the various needed controlling factors for creation of the
  nanotubes.
  '''
  def parseInputFile(self):
    # Open Input Control File and read lines to list.
    f=open('nanoInCtrl','r')
    ctrl=f.read().splitlines()
    f.close()

    # List of all current flags.
    flags=['NAME','IN','TUBERAD','HEIGHT','CHIRALITY','SCROLLING','FUNCTIONALIZATION','OUT','END']


    #Check for END
    if flags[8] in ctrl:

      # Check for NAME  
      if flags[0] in ctrl:
        self.name=ctrl[ctrl.index(flags[0])+1]
      else:
        self.name="MODEL"        
 
      # Check for IN
      if flags[1] in ctrl:
        # If it is a file, get the relevant information.
        if ctrl[ctrl.index(flags[1])+1]=='file':
          nf=olcaomiFile(ctrl[ctrl.index(flags[1])+2])
          nf.getFile()
          self.width=[nf.cellVec[0],nf.cellVec[1],nf.cellVec[2]]
          self.addEntireStructure(nf.atomicList,nf.positionList)

        # If using the adhoc method of list, take in the
        #  relevant parameters and coordinates.
        elif ctrl[ctrl.index(flags[1])+1]=='list':
          i=ctrl.index(flags[1])+2
          line=ctrl[i].split()
          self.width=[float(line[0]),float(line[1]),float(line[2])]
          i+=1
          while(ctrl[i] not in flags):
            line=ctrl[i].split()        
            self.addAtom(line[0],[float(line[1]),float(line[2]),float(line[3])])
            i+=1 
      else:
        sys.exit('IN must be included in the input control file.\n')
 
      # Check for TUBERAD
      if flags[2] in ctrl:
        self.tuberadbounds=map(int, ctrl[ctrl.index(flags[2])+1].split())
      else:
        self.tuberadbounds=[15,80]

      # Check for HEIGHT
      if flags[3] in ctrl:
        self.height=int(ctrl[ctrl.index(flags[3])+1])

      # Check for CHIRALITY
      if flags[4] in ctrl:
        self.chir=float(ctrl[ctrl.index(flags[4])+1])

      # Check for SCROLLING
      if flags[5] in ctrl:
        line=ctrl[ctrl.index(flags[5])+1].split()

        if len(line)==2:
          self.scrdist=float(line[0])
          self.scrnum=int(line[2])
        else:
          self.scrdist=float(line[0])
 
      # Check for FUNCTIONALIZATION
      if flags[6] in ctrl:
        self.func=1
        self.funcLoc=int(ctrl[ctrl.index(flags[6])+1]) 
        line=ctrl[ctrl.index(flags[6])+2]

        # Check for 'list' specification

        if line=='list':
          line2=ctrl[ctrl.index(flags[6])+3]
          # Nanoparticle Designation
          if len(line2)==1:
            self.funcNumPart=int(line2[0]) 
            i=ctrl.index(flags[6])+4
            line3=ctrl[i].split()
            self.funcPartWidth=[float(line3[0]),float(line3[1]),float(line3[2])]
            i+=1
            line3=ctrl[i].split()
            while(line3[0] not in flags):
              self.funccorr_atoms.append(line3[0])
              self.funcinitPos.append([line3[1],line3[2],line3[3]])
              i+=1
              line3=ctrl[i].split()

          # Single Element Designation
          else:
            i=ctrl.index(flags[6])+3
            line2=ctrl[i].split()

            while(line2[0] not in flags):
              self.funccorr_atoms.append(line2[1])
              self.funcElemNums.append(int(line2[0]))
              i+=1
              line2=ctrl[i].split()
           
        # Do file specification when list isn't present.
        elif line=='file':
          line2=ctrl[ctrl.index(flags[6])+3].split()
          self.funcNumPart=int(line2[0])
          nap=olcaomiFile(line2[1])
          nap.getFile()
          self.funccorr_atoms.extend(nap.atomicList)
          self.funcinitPos.extend(nap.positionList)
          self.funcPartWidth=[nap.cellVec[0],nap.cellVec[1],nap.cellVec[2]]

      # Check for OUT
      if flags[7] in ctrl:
        self.outfile=ctrl[ctrl.index(flags[7])+1]
      else:
        self.outfile="olcao.skl"

    else:
      sys.exit("END must be included in the input control file.\n") 



    # Do a quick print out of the input parameters.

    printLine()
    print "INPUT SUMMARY:"
    printLine()

    print "Initial Structure:"
    print "  Atoms in model:\t"+str(len(self.corr_atoms))
    print "  Atoms Names:\t\t"+' '.join(set(self.corr_atoms))
    print "  Cell Dimensions:\t"+' '.join(str(e) for e in self.width)

    printLine()

    print "Nanotube Properties:"
    print "  Tube Radii:\t"+' '.join(str(e) for e in self.tuberadbounds)
    print "  Tube Height:\t"+str(self.height)
    print "  Chirality:\t"+str(self.chir)
    print "  Scrolling:"
    print "    Number of turns:\t"+str(self.scrnum)
    print "    Scroll distance:\t"+str(self.scrdist)
    print "  Functionalization:"
    if self.func==1:
     if len(self.funcElemNums)==0:
       print "    Nanoparticle Decoration"
       if(self.funcPartWidth[2] > self.width[0]*self.height or self.funcPartWidth[2] > self.width[1]*self.height or self.funcPartWidth[2] > self.width[2]*self.height):
         print 
         sys.exit("       ERROR: Nanoparticle height is bigger than tube height. Please change the tube height by an integer value.")
       print "      Number of Nanoparticles:\t"+str(self.funcNumPart)
       print "      Nanoparticle Dimensions:\t"+' '.join(str(e) for e in self.funcPartWidth)
       print "      Atoms in Nanoparticle:\t"+str(len(self.funcinitPos))
       print "      Atom Names:\t\t"+' '.join(set(self.funccorr_atoms))


     else:
       print "    Single Atom Decoration"
       if self.funcLoc==1:
         print "    Location:\t\tInside"
       elif self.funcLoc==2:
         print "    Location:\t\tOutside"
       else:
         print "    Location:\t\tBoth"
       print "    Number of atoms:\t"+str(sum(self.funcElemNums))
       print "    Atom names:\t\t"+' '.join(str(e) for e in self.funccorr_atoms)
    else:
      print "    None"

    printLine()

    print "Output file: "+self.outfile

    printLine()


  # DEFINE THE FOLLOWING ORIENTATIONS 
  def setinitPos1(self):
    self.orientpos1=np.transpose(np.matrix(self.initPositions,dtype='d'))
    self.orientpos1[1,:]=self.orientpos1[1,:]-self.width[1]/2.0
    self.width1=copy.deepcopy(self.width)


    self.width1a=[self.width1[0],self.width1[2],self.width1[1]]
    self.orientpos1a=np.transpose(np.matrix(self.initPositions,dtype='d'))
    self.orientpos1a=np.dot(Rotx(np.pi/2.0),self.orientpos1a)
    self.orientpos1a[1,:]=self.orientpos1a[1,:]+self.width1a[1]/2.0
        

  def setinitPos2(self):
    self.orientpos2=np.transpose(np.matrix(self.initPositions,dtype='d'))
    self.orientpos2=np.dot(Rotz(np.pi/2.0),self.orientpos2)
    self.width2=[self.width[1],self.width[0],self.width[2]]
    self.orientpos2[1,:]=self.orientpos2[1,:]+self.width2[1]/2.0

   
    self.width2a=[self.width2[0],self.width2[2],self.width2[1]]
    self.orientpos2a=copy.deepcopy(self.orientpos2)
    self.orientpos2a=np.dot(Rotx(np.pi/2.0),self.orientpos2a)
    self.orientpos2a[2,:]=self.orientpos2a[2,:]+self.width2a[2]/2.0
    self.orientpos2a[1,:]=self.orientpos2a[1,:]+self.width2a[1]/2.0



  def setinitPos3(self):
    self.orientpos3=np.transpose(np.matrix(self.initPositions,dtype='d'))
    self.orientpos3=np.dot(Rotz(np.pi),self.orientpos3)
    self.width3=copy.deepcopy(self.width)
    self.orientpos3[0,:]=self.orientpos3[0,:]+self.width[0]
    self.orientpos3[1,:]=self.orientpos3[1,:]+(self.width[1]/2.0)


    self.width3a=[self.width3[0],self.width3[2],self.width3[1]]
    self.orientpos3a=copy.deepcopy(self.orientpos3)
    self.orientpos3a=np.dot(Rotx(np.pi/2.0),self.orientpos3a)
    self.orientpos3a[2,:]=self.orientpos3a[2,:]+self.width3a[2]/2.0
    self.orientpos3a[1,:]=self.orientpos3a[1,:]+self.width3a[1]/2.0

  def setinitPos4(self):
    self.orientpos4=np.transpose(np.matrix(self.initPositions,dtype='d'))
    self.orientpos4=np.dot(Rotz(0.75*np.pi),self.orientpos4)
    self.width4=[self.width[1],self.width[0],self.width[2]]
    self.orientpos4[0,:]=self.orientpos4[0,:]+self.width4[0]   
    self.orientpos4[1,:]=self.orientpos4[1,:]-self.width4[1]/2.0

    self.width4a=[self.width4[0],self.width4[2],self.width4[1]]
    self.orientpos4a=copy.deepcopy(self.orientpos4)
    self.orientpos4a=np.dot(Rotx(np.pi/2.0),self.orientpos4a)
    self.orientpos4a[2,:]=self.orientpos4a[2,:]+self.width4a[2]/2.0
    self.orientpos4a[1,:]=self.orientpos4a[1,:]+self.width4a[1]/2.0


  def setinitPos5(self):
    self.orientpos5=np.transpose(np.matrix(self.initPositions,dtype='d'))
    self.orientpos5=np.dot(Roty(np.pi/2),self.orientpos5)
    self.width5=[self.width[2],self.width[1],self.width[0]]
    self.orientpos5[0,:]=self.orientpos5[0,:]+self.width5[0]
    self.orientpos5[1,:]=self.orientpos5[1,:]-self.width5[1]/2.0

    self.width5a=[self.width5[0],self.width5[2],self.width5[1]]
    self.orientpos5a=copy.deepcopy(self.orientpos5)
    self.orientpos5a=np.dot(Rotx(np.pi/2.0),self.orientpos5a)
    self.orientpos5a[2,:]=self.orientpos5a[2,:]+self.width5a[2]/2.0
    self.orientpos5a[1,:]=self.orientpos5a[1,:]+self.width5a[1]/2.0



  def setinitPos6(self):
    self.orientpos6=np.transpose(np.matrix(self.initPositions,dtype='d'))
    self.orientpos6=np.dot(Roty(0.75*np.pi),self.orientpos6)
    self.width6=[self.width[2],self.width[1],self.width[0]]
    self.orientpos6[2,:]=self.orientpos6[2,:]+self.width6[2]
    self.orientpos6[1,:]=self.orientpos6[1,:]-self.width6[1]/2.0

    self.width6a=[self.width6[0],self.width6[2],self.width6[1]]
    self.orientpos6a=copy.deepcopy(self.orientpos6)
    self.orientpos6a=np.dot(Rotx(np.pi/2.0),self.orientpos6a)
    self.orientpos6a[2,:]=self.orientpos6a[2,:]+self.width6a[2]/2.0
    self.orientpos6a[1,:]=self.orientpos6a[1,:]+self.width6a[1]/2.0



  def setAllPos(self):
    self.setinitPos1()
    self.setinitPos2()
    self.setinitPos3()
    self.setinitPos4()
    self.setinitPos5()
    self.setinitPos6()


  def doRot(self,atompositions,celldims,n):
    tuberad=0.0
    ang=0.0
    rotang=0.0

    for i in range(self.tuberadbounds[0],self.tuberadbounds[1]):
      # Get rotation angle and tube radius
      ang=(np.pi*(i-2.0))/(2.0*(i))
      rotang=2*((np.pi/2)-ang)

      tuberad=(celldims[1]/2.0)*np.tan(ang)
      
      # Make a temporary copy of the atomic positions
      temp=copy.deepcopy(atompositions)
      
      # Check for chirality
      if self.chir>0:
        zshift=(self.chir*celldims[2])/float(i)
        for j in range(len(self.corr_atoms)):
          temp[2,j]=temp[2,j]+((temp[1,j]*zshift)/celldims[1])
          if temp[2,j] > celldims[2]:
            temp[2,j]=temp[2,j]-celldims[2]
          elif temp[2,j]<0.0:
            temp[2,j]=temp[2,j]+celldims[2]


      # Shift the cell x-values out such that the inner cell face is
      #  at the specified radius
      temp[0,:]=temp[0,:]+tuberad


      # Make the tube coordinates the new coordinates
      tubecoords=copy.deepcopy(temp)
      elemlist=list(copy.deepcopy(self.corr_atoms)) 


      # Do n-1 duplication / rotations of the structure until the ring is formed.
      for j in range(1,i):
        temp=np.dot(Rotz(rotang),temp)

        # If chirality is present, keep z-periodicity by moving
        #  anything above or below the cell height.
        if self.chir>0:
          temp[2,:]=temp[2,:]-zshift

          for k in range(len(self.corr_atoms)):
            if temp[2,k]>celldims[2]:
              temp[2,k]=temp[2,k]-celldims[2]
            elif temp[2,k]<0.0:
              temp[2,k]=temp[2,k]+celldims[2]


        # Add the new piece to the existing structure.
        tubecoords=np.append(tubecoords,temp,1)
        elemlist.extend(list(copy.deepcopy(self.corr_atoms)))


 
      # Extend the nanotube in the z-direction based on the height.
      if self.height>1:
        temp2=copy.deepcopy(tubecoords)
        temp3=copy.deepcopy(elemlist)

        for k in range(self.height+1):
          temp2[2,:]=temp2[2,:]+celldims[2]
        
          tubecoords=np.append(tubecoords,temp2,1)
          elemlist.extend(list(copy.deepcopy(temp3)))


      newtubeheight=celldims[2]*self.height

      # Functionalize the nanotube
      if self.func==1:
        # Define functionalization radii
        delx=0.0025
        funcRadiusIn=tuberad-delx
        funcRadiusOut=tuberad+celldims[0]+delx

        # Check for functionalization location
        # INSIDE
        if self.funcLoc==1:

          # Single Elems or Nanoparticle.
          #  NANOPARTICLE
          if len(self.funcElemNums)==0:
            funcAngs=np.random.uniform(low=0.0,high=2.0*np.pi,size=(self.funcNumPart,))
            xpos=funcRadiusIn*np.cos(funcAngs)
            ypos=funcRadiusIn*np.sin(funcAngs)
            zpos=np.random.uniform(low=0.0,high=newtubeheight-self.funcPartWidth[2],size=(self.funcNumPart,))
           
            for j in range(len(self.funcNumPart)):
              tempPartCoords=copy.deepcopy(np.transpose(np.matrix(self.funcinitPos,dtype='d')))
              tempPartCoords[1,:]=tempPartCoords[1,:]-self.funcPartWidth[1]/2.0
              tempPartCoords[0,:]=tempPartCoords[0,:]+(funcRadiusIn-self.funcPartWidth[0])
              tempPartCoords=np.dot(Rotz(funcAngs[j]),tempPartCoords)
              tempPartCoords[2,:]=tempPartCoords[2,:]+zpos[j]
              tubecoords=np.append(tubecoords,tempPartCoords,1)
              elemlist.extend(self.funccorr_atoms)
 
          # SINGLE ELEMENTS
          else:  
            for j in range(len(self.funccorr_atoms)):
              for k in range(self.funcElemNums[j]):
                funcAngs=np.random.uniform(low=0.0,high=2.0*np.pi,size=(1,))
                xpos=funcRadiusIn*np.cos(funcAngs)
                ypos=funcRadiusIn*np.sin(funcAngs)
                zpos=np.random.uniform(low=0.0,high=newtubeheight,size=(1,))

                quickvec=np.zeros((3,1),dtype='d')
                quickvec[0,0]=xpos[0]
                quickvec[1,0]=ypos[0]
                quickvec[2,0]=zpos[0]

                tubecoords=np.append(tubecoords,quickvec,1)
                elemlist.append(self.funccorr_atoms[j])
             
        # OUTSIDE
        elif self.funcLoc==2:
          # Single Elems or Nanoparticle.
          # NANOPARTICLE
          if len(self.funcElemNums)==0:
            funcAngs=np.random.uniform(low=0.0,high=2.0*np.pi,size=(self.funcNumPart,))
            xpos=funcRadiusOut*np.cos(funcAngs)
            ypos=funcRadiusOut*np.sin(funcAngs)
            zpos=np.random.uniform(low=0.0,high=newtubeheight-self.funcPartWidth[2],size=(self.funcNumPart,))

            for j in range(self.funcNumPart):
              tempPartCoords=np.transpose(np.matrix(self.funcinitPos,dtype='d'))
              tempPartCoords[1,:]=tempPartCoords[1,:]-self.funcPartWidth[1]/2.0
              tempPartCoords[0,:]=tempPartCoords[0,:]+funcRadiusOut
              tempPartCoords=np.dot(Rotz(funcAngs[j]),tempPartCoords)
              tempPartCoords[2,:]=tempPartCoords[2,:]+zpos[j]
              tubecoords=np.append(tubecoords,tempPartCoords,1)
              elemlist.extend(self.funccorr_atoms)

          # SINGLE ELEMENTS
          else:  
            for j in range(len(self.funccorr_atoms)):
              for k in range(self.funcElemNums[j]):
                funcAngs=np.random.uniform(low=0.0,high=2.0*np.pi,size=(1,))
                xpos=funcRadiusOut*np.cos(funcAngs)
                ypos=funcRadiusOut*np.sin(funcAngs)
                zpos=np.random.uniform(low=0.0,high=newtubeheight,size=(1,))

                quickvec=np.zeros((3,1),dtype='d')
                quickvec[0,0]=xpos[0]
                quickvec[1,0]=ypos[0]
                quickvec[2,0]=zpos[0]
                tubecoords=np.append(tubecoords,quickvec,1)
                elemlist.append(self.funccorr_atoms[j])


        # BOTH
        elif self.funcLoc==3:
          # Single Elems or Nanoparticle.
          if len(self.funcElemNums)==0:
            for j in range(self.funcNumPart+1):
              if np.random.randint(2)==0:
                funcAngs=np.random.uniform(low=0.0,high=2.0*np.pi,size=(1,))
                xmove=funcRadiusIn-self.funcPartWidth[0]
                zpos=np.random.uniform(low=0.0,high=newtubeheight-self.funcPartWidth[2],size=(1,))
              
              else: 
                funcAngs=np.random.uniform(low=0.0,high=2.0*np.pi,size=(1,))
                xmove=funcRadiusOut
                zpos=np.random.uniform(low=0.0,high=newtubeheight-self.funcPartWidth[2],size=(1,))
         
              tempPartCoords=copy.deepcopy(np.transpose(np.matrix(self.funcinitPos,dtype='d')))
              tempPartCoords[1,:]=tempPartCoords[1,:]-self.funcPartWidth[1]/2.0
              tempPartCoords[0,:]=tempPartCoords[0,:]+xmove
              tempPartCoords=np.dot(Rotz(funcAngs[0]),tempPartCoords)
              tempPartCoords[2,:]=tempPartCoords[2,:]+zpos[0]

              tubecoords=np.append(tubecoords,tempPartCoords,1)
              elemlist.extend(self.funccorr_atoms)
              
          else:
            for j in range(len(self.funccorr_atoms)):
              for k in range(self.funcElemNums[j]):
                funcAngs=np.random.uniform(low=0.0,high=2.0*np.pi,size=(1,))
                zpos=np.random.uniform(low=0.0,high=newtubeheight,size=(1,))

                if np.random.randint(2)==0:
                  xpos=funcRadiusIn*np.cos(funcAngs)
                  ypos=funcRadiusIn*np.sin(funcAngs)
                else:
                  xpos=funcRadiusOut*np.cos(funcAngs)
                  ypos=funcRadiusOut*np.sin(funcAngs)

                quickvec=np.zeros((3,1),dtype='d')
                quickvec[0,0]=xpos[0]
                quickvec[1,0]=ypos[0]
                quickvec[2,0]=zpos[0]


                tubecoords=np.append(tubecoords,quickvec,1)
                elemlist.append(self.funccorr_atoms[j])
                  


      # STARTING OUTPUT HERE

      xmax=0.0
      ymax=0.0

      newcellwidth=2.0*tuberad+2.0*celldims[0]+30.0
      fileout=open('NT_'+self.name+'_Face'+str(n)+'_'+str(i)+'C'+str(self.chir).replace('.','p')+'_Func'+str(self.func)+'.dat','w')

      tubecoords[0,:]=tubecoords[0,:]+newcellwidth/2.0
      tubecoords[1,:]=tubecoords[1,:]+newcellwidth/2.0


      fileout.write('title\n')
      fileout.write(self.name+' '+str(i)+'\n')
      fileout.write('end\n')
      fileout.write('cell\n')

      fileout.write(str(newcellwidth)+'  '+str(newcellwidth)+'  '+str(newtubeheight)+' 90.000 90.000 90.000\n')
      fileout.write('cartesian '+str(len(elemlist))+'\n')
      for j in range(len(elemlist)):
        fileout.write(elemlist[j]+'\t'+str(tubecoords[0,j])+'\t'+str(tubecoords[1,j])+'\t'+str(tubecoords[2,j])+'\n')
      fileout.write('space 1_a\n')
      fileout.write('supercell 1 1 1\n')
      fileout.write('full\n')
      fileout.close()




  def doAllRots(self):
    self.doRot(self.orientpos1,self.width1,1)
    self.doRot(self.orientpos1a,self.width1a,'1a')
    self.doRot(self.orientpos2,self.width2,2)
    self.doRot(self.orientpos2a,self.width2a,'2a')
    self.doRot(self.orientpos3,self.width3,3)
    self.doRot(self.orientpos3a,self.width3a,'3a')
    self.doRot(self.orientpos4,self.width4,4)
    self.doRot(self.orientpos4a,self.width4a,'4a')
    self.doRot(self.orientpos5,self.width5,5)
    self.doRot(self.orientpos5a,self.width5a,'5a')
    self.doRot(self.orientpos6,self.width6,6)
    self.doRot(self.orientpos6a,self.width6a,'6a')






print "\n" 
printLine()
printLine()
print "Beginning Nanotube Creation"
printLine()

'''     
Basic Idea:
 1) Acquire the structure
 2) Acquire desired parameters
 3) Place initial piece with correct orientation and at correct distance
 4) Make first ring
 5) Replicate ring upwards
'''

tube=initStructure()
tube.parseInputFile()
tube.setAllPos()
tube.doAllRots()


print "Time to complete "+str(time.time()-start)+"s."
printLine()
printLine()



#    fig=plt.figure()
#    ax=fig.add_subplot(111,projection='3d')
#
#    for k in range(len(self.corr_atoms)):
#      ax.scatter(self.orientpos3[0,k],self.orientpos3[1,k],self.orientpos3[2,k])
#
#    plt.show()
