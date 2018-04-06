#!/usr/bin/python

"""
    ConfBuster-Single-Molecule-Minimization.py
    Copyright (C) 2017  Xavier Barbeau, Antony T. Vincent and Patrick Lague

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/.
"""

import os, sys, math, random, getopt
import filecmp

def main(argv):
    ### arguments
    
    PYMOL_PATH = '' # path to the pymol executable. Leave empty if in $PATH 
    BABEL_PATH = '' # path to the babel executables. Leave empty if in $PATH

    steps           = 5000
    FF              = 'MMFF94s'
    in_name = ''
    threads = 2
    mandatory = ''
    out_name_prefix = ''
    try:
        opts, args = getopt.getopt(argv,"hi:T:o:F:s:")
    except getopt.GetoptError:
        print 'ConfBuster-Single-Molecule-Minimization.py -i inputfile [mandatory] -o output name prefix [default: replace input file]'
        sys.exit(2)

    for opt, arg in opts:
        print "%s %s" %(opt, arg)
        if opt == '-h':
            print 
            print 'ConfBuster-Single-Molecule-Minimization.py -i inputfile [mandatory] -o output name prefix [default: replace input file]'
            sys.exit()
        elif opt in ("-i"):
            in_name = arg
            mandatory = '1'
        elif opt in ("-o"):
            out_name_prefix = arg
        elif opt in ("-T"):
            threads = str(arg)
        elif opt in ("-s"):
            steps = int(arg)
    if mandatory != '1':
        print 'ConfBuster-Single-Molecule-Minimization.py -i inputfile [mandatory] -o output name prefix [default: replace input file]'
        sys.exit()

    BABELKEYWORD    = " -c 1e-8 -sd -n %s "%(steps)


    # remove WARNINGS from xyz file
    def remove_warning_xyz(tmpfilename):
        tmp_f = open(tmpfilename,'r')
        tmp_r = tmp_f.read()
        tmp_f.close()
        tmp_s = tmp_r.splitlines()
    
        tmp_out = open(tmpfilename,'w')
        for x in tmp_s:
            if "WARNING" not in x:
                tmp_out.write(x+'\n')
        tmp_out.close()

    command = '%sobminimize -ff %s %s %s > tmp.pdb'%(BABEL_PATH,FF, BABELKEYWORD,in_name)
    print command
    os.system(command)
    
    #remove_warning_xyz('tmp.xyz')   
 
    command = '%sobenergy -ff %s tmp.pdb > tmp.energy'%(BABEL_PATH,FF)
    os.system(command)
    command = '%sobabel tmp.pdb -O tmp.mol2'%(BABEL_PATH)
    os.system(command)
    
    f1 = open('tmp.mol2','r')
    r1 = f1.read()
    f1.close()
    s1 = r1.splitlines()
    
    tmp_f = open('tmp.energy','r')
    tmp_r = tmp_f.read()
    tmp_f.close()
    tmp_s = tmp_r.splitlines()
    for x in tmp_s:
        if 'TOTAL ENERGY =' in x:
            energy = x.split()[3]
    
    #print 'energy is : ',energy
    
    out = open('tmp.mol2','w')
    for i in range(len(s1)):
        if i != 1:
            out.write(s1[i]+'\n')
        else:
            out.write(energy+'\n')
    out.close()

    if out_name_prefix == '':
        try:
            outname = in_name.rsplit('/',1)[1].rsplit('.',1)[0]
        except:
            outname = in_name.rsplit('.',1)[0]
    else:
        outname = out_name_prefix
    
    command = 'cp tmp.mol2 %s.mol2'%(outname)
    os.system(command)
    
    

if __name__ == "__main__":
    main(sys.argv[1:])

