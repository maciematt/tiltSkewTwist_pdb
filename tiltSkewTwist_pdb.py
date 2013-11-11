#!/usr/bin/env python


"""
This script takes in pdb's and spits out the skequential tilt, skew, and twist angles between component domains, as defined by the user. It requires numpy to run.
Usage: ./tiltSkewTwist.py [infile] [dom1 residues] [dom1 atoms] [dom1 ref] [dom2 residues] [dom2 atoms] [dom2_ref]
[dom1/dom2 residues] is the residue range of a given domain given in PyMOL format, i.e. low1-high1+low2-high2 etc., same goes for [dom1/dom2 atoms], eg. n+ca+c (write 'all' to include all atom types), and [dom1/dom2 ref] is the seq. number of the reference residue.

Note: review the code to choose the coordinate set type (right- or left-handed; right-handed by default), the angle convention (0 - 360 degrees or -180 - 180 degrees; the latter by default) etc., as annotated within the code.

"""


__author__ = "Mateusz Maciejewski"
__email__ = "matt@mattmaciejewski.com"
__date__ = "December 2010"
__copyright__ = "Copyright (c) 2010, 2011, 2012 Mateusz Maciejewski, Barlow Lab, Biomolecular NMR Unit, University of Edinburgh"
__license__ = "MIT"


import sys
import re
import numpy
from math import cos, sin, sqrt, acos, degrees, pi


if len(sys.argv) != 8:
    print >> sys.stderr, __doc__
    sys.exit(1)



class AtomCollector:

    """
    This class collects info on each atom in target pdb.
    It creates a molecule object that stores atoms as attributes identified
    by atom numbers from pdb:
    .model - molecular model identifier
    .resid - type of residue the atom belongs to
    .chain - chain ID
    .residnum - sequence number of the residue
    .xcoor - x coordinate
    .ycoor - y coordinate
    .zcoor - z coordinate
    .occup - atom occupancy
    .tempfact - temperature factor
    .mass - atomic mass
    """

    def __init__(self):

        self.atom = {}
        self.atomNumber = {}
        self.resid = {}
        self.chain = {}
        self.residnum = {}
        self.xcoor = {}
        self.ycoor = {}
        self.zcoor = {}
        self.occup = {}
        self.tempfact = {}
        self.mass = {}
        self.model = {}
        self.atomCounter = 0
        whiteSpace = re.compile(r'\s+')

    def collect(self,fileToWorkOn):
                
        pattern1 = re.compile('^ATOM')
        pattern2 = re.compile('^MODEL')

        for line in fileToWorkOn.readlines():

            if pattern2.search(line):
                current_model = int(line.split()[1])
            
            if pattern1.search(line):
                
                    
                self.atomNumber[self.atomCounter] = line[6:11]
                self.atom[self.atomCounter] = line[12:16].strip()
                self.resid[self.atomCounter] = line[17:20].strip()
                self.chain[self.atomCounter] = line[21].strip()
                self.residnum[self.atomCounter] = int(line[22:26])
                self.xcoor[self.atomCounter] = float(line[30:38])
                self.ycoor[self.atomCounter] = float(line[38:46])
                self.zcoor[self.atomCounter] = float(line[46:54])
                self.occup[self.atomCounter] = line[54:60]
                self.tempfact[self.atomCounter] = line[60:66]

                try:
                    self.model[self.atomCounter] = current_model
                except NameError:
                    self.model[self.atomCounter] = 1
        

                if self.atom[self.atomCounter][0].lower() == 'c': self.mass[self.atomCounter] = 12.0
                elif self.atom[self.atomCounter][0].lower() == 'o': self.mass[self.atomCounter] = 16.0
                elif self.atom[self.atomCounter][0].lower() == 'n': self.mass[self.atomCounter] = 15.0
                elif self.atom[self.atomCounter][0].lower() == 's': self.mass[self.atomCounter] = 32.0
                elif self.atom[self.atomCounter][0].lower() == 'h': self.mass[self.atomCounter] = 1.0
                elif len(self.atom[self.atomCounter]) > 1 and self.atom[self.atomCounter][1].lower() == 'h': self.mass[self.atomCounter] = 1.0
                elif self.atom[self.atomCounter][0].lower() == 'd': self.mass[self.atomCounter] = 2.0
                elif self.atom[self.atomCounter][0].lower() == '2': self.mass[self.atomCounter] = 1.0
                elif self.atom[self.atomCounter][0].lower() == '4': self.mass[self.atomCounter] = 1.0
                elif self.atom[self.atomCounter][0].lower() == 'p': self.mass[self.atomCounter] = 31.0
                elif self.atom[self.atomCounter][0].lower() == 'z': self.mass[self.atomCounter] = 65.0
                elif self.atom[self.atomCounter][0].lower() == 'f': self.mass[self.atomCounter] = 55.9

            
                #self.mass[self.atomCounter] = 1.0 # uniform mass switch

                self.atomCounter += 1


fileToWorkOn = open(sys.argv[1],'r')

id_mat = numpy.array([(1,0,0),(0,1,0),(0,0,1)])

molecule = AtomCollector()
molecule.collect(fileToWorkOn)


def norm(vekt):
    normal = sqrt(vekt[0]**2+vekt[1]**2+vekt[2]**2)
    return vekt[0]/normal, vekt[1]/normal, vekt[2]/normal

def angle(a,b):
    return acos(a[0]*b[0]+a[1]*b[1]+a[2]*b[2])    

def rotate(rot1,rot2):
    rotation = numpy.matrix([(cos(rot1),(-1)*sin(rot1),0),
                           (sin(rot1),cos(rot1),0),
                           (0,0,1)])
    rotation = rotation * numpy.matrix([(cos(rot2),0,(-1)*sin(rot2)),
                                        (0,1,0),
                                        (sin(rot2),0,cos(rot2))])
    return rotation.T



def parseInput(resids, atoms):
    
    if resids != "all":
        residues = []
        if resids.find("+") != -1:
            resids = resids.split("+")

            for res in resids:
                if res.find("-") != -1:
                    residues += range(int(res.split("-")[0]),int(res.split("-")[1])+1)
                else:
                    residues += [int(res)]
        else:
            if resids.find("-") != -1:
                residues += range(int(resids.split("-")[0]),int(resids.split("-")[1])+1)
            else:
                residues += [int(resids)]
    else:
        residues = "all"

    if atoms != "all":
        atoms = atoms.split("+")

    return residues, atoms

def calcTensor(res, ats, trp, mod, molecule):
    
    COM = [0,0,0]    #center of mass
    spec = [0,0,0]
    comdiv = 0

    for entry in molecule.atom.keys():
        if (int(molecule.residnum[entry]) in res) and \
           (molecule.atom[entry].lower() in ats or ats == "all") and \
           (molecule.resid[entry] != "ANI") and molecule.model[entry] == mod:

            COM[0] += molecule.xcoor[entry]*molecule.mass[entry]
            COM[1] += molecule.ycoor[entry]*molecule.mass[entry]
            COM[2] += molecule.zcoor[entry]*molecule.mass[entry]
            comdiv += molecule.mass[entry]

        if int(molecule.residnum[entry]) == trp and molecule.atom[entry].lower() == "ca" and \
                molecule.model[entry] == mod:
            spec[0] = molecule.xcoor[entry]
            spec[1] = molecule.ycoor[entry]
            spec[2] = molecule.zcoor[entry]

        if int(molecule.residnum[entry]) == min(res) and \
                molecule.atom[entry].lower() == "ca" and molecule.model[entry] == mod:
            N_term = [molecule.xcoor[entry], molecule.ycoor[entry], molecule.zcoor[entry]]

        if int(molecule.residnum[entry]) == max(res) and \
                molecule.atom[entry].lower() == "ca" and molecule.model[entry] == mod:
            C_term = [molecule.xcoor[entry], molecule.ycoor[entry], molecule.zcoor[entry]]

    orient_vect = norm([C_term[0]-N_term[0], C_term[1]-N_term[1], C_term[2]-N_term[2]])


    COM[0] /= comdiv
    COM[1] /= comdiv
    COM[2] /= comdiv


    spec[0] -= COM[0]
    spec[1] -= COM[1]
    spec[2] -= COM[2]

    norm(spec)

    I=[]

    for index in range(9):
        I.append(0)

    for entry in molecule.atom.keys():

        if (int(molecule.residnum[entry]) in res) and \
           (molecule.atom[entry].lower() in ats or ats == "all") and \
           (molecule.resid[entry] != "ANI") and molecule.model[entry] == mod:
            
            temp_x, temp_y, temp_z = molecule.xcoor[entry], molecule.ycoor[entry], molecule.zcoor[entry]
            temp_x-=COM[0]; temp_y-=COM[1]; temp_z-=COM[2]

            I[0] += molecule.mass[entry] * (temp_y**2 + temp_z**2)
            I[4] += molecule.mass[entry] * (temp_x**2+ temp_z**2)
            I[8] += molecule.mass[entry] * (temp_x**2+ temp_y**2)
            I[1] -= molecule.mass[entry] * temp_x * temp_y
            I[3] -= molecule.mass[entry] * temp_x * temp_y
            I[2] -= molecule.mass[entry] * temp_x * temp_z
            I[6] -= molecule.mass[entry] * temp_x * temp_z
            I[5] -= molecule.mass[entry] * temp_y * temp_z
            I[7] -= molecule.mass[entry] * temp_y * temp_z

    tensor = numpy.array([(I[0:3]),(I[3:6]),(I[6:9])])
    vals,vects = numpy.linalg.eig(tensor)
    eig_ord = numpy.argsort(vals) # COLUMN i corrensponds to eigenvalue i.
    ord_vects = vects[:,eig_ord].T
    ref = ord_vects[0]


    if degrees(angle(ref, orient_vect)) > 90:
        for index in range(3): ref[index] *= (-1)


    return ref, spec
            

domainsToDo = 2

unique_models = list(set(molecule.model.values())) # this will remove the repetitions from values upon creation of the set, and then turn that into a list


for mod in unique_models:

    ref = []
    spec = []
    refT = []
    specT = []

    for i in range(domainsToDo):

        b = i*3+2
        res, ats = parseInput(sys.argv[b], sys.argv[b+1])
        trp = int(sys.argv[b+2])

        ref_ten, spec_ten = calcTensor(res, ats, trp, mod, molecule)
        ref.append(ref_ten)
        spec.append(spec_ten)

    temp = []

    refNormd = [norm(ref[ii]) for ii in range(domainsToDo)]
    specNormd = [norm(spec[ii]) for ii in range(domainsToDo)]

    print sys.argv[1], "model", mod,


    y_ax = norm(numpy.cross(refNormd[0],specNormd[0])) # y axis for reference (based on z and Trp)
    x_ax = norm(numpy.cross(y_ax,refNormd[0]))  # the true x axis for reference

    for index in range(3): temp.append(x_ax[index]*(-1))   # uncomment these two lines to use a left-handed coordinate set, otherwise it's right-handed.
    x_ax = temp[:]

    temp = []

    test_y = norm(numpy.cross(refNormd[1],specNormd[1]))
    test_x = norm(numpy.cross(test_y,refNormd[1]))

    for index in range(3): temp.append(test_x[index]*(-1))   # uncomment these two lines to use a left-handed coordinate set, otherwise it's right-handed.
    text_x = temp[:]

    coord_ref = numpy.array([x_ax,y_ax,refNormd[0]]).T
    trans = numpy.matrix(coord_ref).T

    # In the next move trans is used for transforming angles into
    # Cartesian set (M*M.T = I)

    refT.append( numpy.matrix(norm((numpy.matrix(trans) * numpy.matrix(refNormd[0]).T).A[:,0])) )
    refT.append( numpy.matrix(norm((numpy.matrix(trans) * numpy.matrix(refNormd[1]).T).A[:,0])) )

    specT.append( numpy.matrix(norm((numpy.matrix(trans) * numpy.matrix(specNormd[0]).T).A[:,0])) )
    specT.append( numpy.matrix(norm((numpy.matrix(trans) * numpy.matrix(specNormd[1]).T).A[:,0])) )

    plane_xy = norm(numpy.cross((1,0,0),(0,1,0)))
    plane_zx = norm(numpy.cross((0,0,1),(1,0,0)))
    plane_yz = norm(numpy.cross((0,1,0),(0,0,1)))

    ref_on_plane_xy = norm(numpy.cross(plane_xy,numpy.cross(plane_xy,refT[1])[0]))

    twist = angle(id_mat[0],ref_on_plane_xy)
    aux = angle(id_mat[1],ref_on_plane_xy)

    if aux > pi/2:
    #    twist = pi*2 - twist # one of the possible conventions
        twist = (-1) * twist  # one of the possible conventions
        pass

    tilt = angle(refT[0].T,refT[1].T)

    rotation = rotate(twist,tilt)

    spec_r = rotation*specT[1].T
    spec_on_plane_zx = norm(numpy.cross(plane_zx,numpy.cross(plane_zx,spec_r.A[:,0]))) # this filters out the xy-plane projection of spec1

    skew = angle(plane_yz,spec_on_plane_zx)
    aux = angle(plane_zx,spec_on_plane_zx)

    if aux > pi/2:
    #    skew = pi*2 - skew # one of the possible conventions
        skew = (-1) * skew # one of the possible conventions
        pass


    print "Tilt:", degrees(tilt), "Twist:", degrees(twist),"Skew:", degrees(skew),

    print

