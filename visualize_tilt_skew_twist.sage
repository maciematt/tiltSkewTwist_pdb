
###################################################################
#                                                                 #
#  (c) Mateusz Maciejewski, 2010                                  #
#  License: MIT                                                   #
#  This script requires SAGE (http://www.sagemath.org/) to        #
#  visualize the output of tilt_skew_twist.py:                    #
#                                                                 #
#  sage visualize_tilt_skew_twist.sage                            #
#                                                                 #
#  substitute infilename and outfilename below with desired       #
#  names.                                                         #
#                                                                 #
###################################################################


from math import radians
from os import getcwd

working_dir = getcwd()+"/"
infilename = "tilt_skew_twist.txt"
outfilename = "tilt_skew_twist.eps"

zeros=[]

for line in open(working_dir+infilename).readlines():
    line = [a.strip() for a in line.split()]
    zeros.append([float(line[4]),float(line[6]),float(line[8])])
            

zero_tilts12 = point([(1.3*cos(radians(ang[0])),1.3*sin(radians(ang[0]))) for ang in zeros],rgbcolor='blue', size=25,zorder=1)
zero_skews12 = point([(1.1*cos(radians(ang[1])),1.1*sin(radians(ang[1]))) for ang in zeros],rgbcolor='blue', size=25,zorder=1)
zero_twists12 = point([(0.9*cos(radians(ang[2])),0.9*sin(radians(ang[2]))) for ang in zeros],rgbcolor='blue', size=25,zorder=1)

cir = circle((0,0),1.5,thickness=1.5,rgbcolor=Color('black'))
cir_tilt = circle((0,0),1.3,rgbcolor=Color('black').lighter(0.7),zorder=0)
cir_skew = circle((0,0),1.1,rgbcolor=Color('black').lighter(0.7),zorder=0)
cir_twist = circle((0,0),0.9,rgbcolor=Color('black').lighter(0.7),zorder=0)
lin_x = line2d([(-1.5,0), (1.5,0)],rgbcolor='black')
lin_y = line2d([(0,-1.5), (0,1.5)],rgbcolor='black')

label1 = text('0',(1.65,0),fontsize=20,rgbcolor='black')
label2 = text('90',(0,1.7),fontsize=20,rgbcolor='black')
label3 = text('180',(-1.8,0),fontsize=20,rgbcolor='black')
label4 = text('270',(0,-1.7),fontsize=20,rgbcolor='black')

lab12_tilt = text('Tilt',(-0.22,1.23),fontsize=15,rgbcolor='black')
lab12_skew = text('Skew',(-0.14,1),fontsize=15,rgbcolor='black')
lab12_twist = text('Twist',(-0.26,0.7),fontsize=15,rgbcolor='black')

between_12 = (cir+zero_tilts12+zero_skews12+zero_twists12+\
              label1+label2+label3+label4+lin_x+lin_y+cir_tilt+cir_skew+cir_twist)
between_12.show(figsize=[6,6],xmax=1.9,ymax=1.9,xmin=-1.9,ymin=-1.9,ticks=[[],[]],axes=False)


between_12.save(working_dir+outfilename,figsize=[6,6],xmax=1.9,ymax=1.9,xmin=-1.9,ymin=-1.9,ticks=[[],[]],axes=False,dpi=300)
