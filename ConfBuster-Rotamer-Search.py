#!/usr/bin/python

import os, sys, math, random, getopt
import filecmp

"""
    ConfBuster-Rotamer-Search.py
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


def main(argv):
    ### arguments

    PYMOL_PATH = '' # path to the pymol executable. Leave empty if in path 
    BABEL_PATH = '' # path to the babel executables. Leave empty if in path

    molecule    = '' 
    threads = 2
    number = 100
    E_CUTOFF = 50 
    user_outdir = ''
    mandatory = ''
    outformat = 'mol2'
    FF = 'MMFF94s'
    gentype = 'energy'

    try:
        opts, args = getopt.getopt(argv,"hi:T:g:e:o:f:F:G:")
    except getopt.GetoptError:
        print 'ConfBuster-Rotamer-Search.py -i input filename in mol2 [mandatory] -g # of generations [100] -e energy cutoff to discriminate conformations in kcal/mol [default: 50] -d output directory [default: use the prefix of the input filename] -f format of outputted molecules [xyz,mol2(default)]'
        sys.exit(2)

    for opt, arg in opts:
        print "%s %s" %(opt, arg)
        if opt == '-h':
            print 'ConfBuster-Rotamer-Search.py -i input filename in mol2 [mandatory] -g # of generations [100] -e energy cutoff to discriminate conformations in kcal/mol [default: 50] -d output directory [default: use the prefix of the input filename] -f format of outputted molecules [xyz,mol2(default)]'
            sys.exit()
        elif opt in ("-i"):
            molecule = arg
            mandatory = '1'
        elif opt in ("-g"):
            number = int(arg)
        elif opt in ("-G"):
            if arg in ['energy','rmsd']:
                gentype = str(arg)
            else:
                sys.exit('-G arguments are : energy or rmsd')
        elif opt in ("-e"):
            E_CUTOFF = int(arg)
        elif opt in ("-o"):
            user_outdir = arg
        elif opt in ("-f"):
            outformat = str(arg)
            if outformat not in ['mol2','xyz']:
                mandatory = 0
                print 'Accepted format of output molecules is mol2 or xyz'
    if mandatory != '1':
        print 'ConfBuster-Rotamer-Search.py -i input filename in mol2 [mandatory] -g # of generations [100] -e energy cutoff to discriminate conformations in kcal/mol [default: 50] -d output directory [default: use the prefix of the input filename] -f format of outputted molecules [xyz,mol2(default)]'
        sys.exit()

    
    # setting the number of threads
    os.environ['OMP_NUM_THREADS']="%d" %(threads)
    ENERGY    = " -ff %s "%(FF)   # ff used for energy calculation which is stored in the molecule title
    
    #
    #   1   Generate conformations
    #
   
    outmolecule = molecule.rsplit('.',1)[0]

    if user_outdir != '': 
        outdir = user_outdir+'/'
    else:
        outdir      = 'Rotamers_for_'+outmolecule+'/'
    print outdir
    os.system('mkdir -p %s'%(outdir))

    command = '%sobabel %s -O %s-conformers.mol2 --conformer --nconf %s --score %s --writeconformers --original --log '%(BABEL_PATH,molecule,outmolecule,number,gentype)
    print command
    os.system(command)
    
    command2 = '%sobabel %s-conformers.mol2 -O %s/%s-c.%s -m > /dev/null'%(BABEL_PATH,outmolecule,outdir,outmolecule,outformat)
    #print command2
    os.system(command2)
    
    #
    #   2   get energy of conformations and store in xyz file
    #
    
    print 'Storing molecule energy in individual conformation files'
    
    molecules_energies_list = []
    count = 0
    for x in os.listdir(outdir):
        if '.%s'%(outformat) in x:
            count += 1
            if count % 100 == 0:
                print count
            command = '%sobenergy %s %s 1> tmp.mmff94 2> /dev/null'%(BABEL_PATH,ENERGY,outdir+'/'+x) 
            #print command
            os.system(command)
    
            # read energy tmp.mmff94
    
            energy = 9999999
    
            tmp_f = open('tmp.mmff94','r')
            tmp_r = tmp_f.read()
            tmp_f.close()
            tmp_s = tmp_r.splitlines()
            for y in tmp_s:
                if 'TOTAL ENERGY =' in y:
                    energy = float(y.split()[3])
    
            # write energy to xyz
    
            f = open('%s'%(outdir+'/'+x),'r')
            r = f.read()
            s = r.splitlines()
            f.close()
    
            head = s.pop(0)
            title  = s.pop(0)
    
            out = open('%s'%(outdir+'/'+x),'w')
    
            out.write(head+'\n')
            out.write(str(energy)+'\n')
    
            for y in s:
                out.write(y+'\n')
    
            out.close()
    
    #
    #   3. remove conformation with E > E_CUTOFF
    #
    
    energy_list = []
    
    ls = os.listdir(outdir)
    
    for x in ls:
        if '.%s'%(outformat) in x:
            f = open('%s/%s'%(outdir,x),'r') 
            r = f.read()
            f.close()
            s = r.splitlines()
            e = s[1]
            energy_list.append([float(e),x])
    
    energy_list.sort() # sorting. best energy is the lowest, thus the first of list
    
    #
    #   4. get lowest energy. remove molecules with energy E_CUTOFF kcal higher
    #
    
    min_e = 1000000
    
    
    for x in energy_list:
        e = x[0]
        if e < min_e:
            min_e = e
    
    over = 0
    
    for x in energy_list:
        e = x[0]
        if e > min_e + E_CUTOFF:
            over += 1
            command = 'rm %s/%s'%(outdir,x[1])
            os.system(command)
    
    print 'Lowest energy is : %s'%(str(min_e))
    print 'There are %s conformations with energy %s kcal higher than lowest energy conformation. \nThese conformations were removed.'%(str(over), str(E_CUTOFF))
    
    
    ######  write pymol follow-up code
    
    Pymol_follow = """                              
cmd.reinitialize()
ls = os.listdir('.')                   

load_order = []

for x in ls:                     
    if '.%s' in x:
        f = open(x,'r')
        r = f.read()
        f.close()
        s = r.splitlines()
        E = float(s[1])
        load_order.append([E,x])

load_order.sort()
                
refname = load_order[0][1].rsplit('.',1)[0]

cmd.load(load_order[0][1])
                             
for x in range(len(load_order)-1):                     
    cmd.load(load_order[x+1][1])                
    cmd.pair_fit(load_order[x+1][1].rsplit('.',1)[0],load_order[0][1].rsplit('.',1)[0])
    loaded = load_order[x][1].rsplit('.',1)[0]     

cmd.center()
cmd.color('gray')
cmd.do('util.cnc')
cmd.color('green','elem c and '+refname)
cmd.show('sticks',refname)
cmd.set('stick_radius','0.1')
    """%(outformat)
                                                 
    write_pymol_script = open('%s/Follow-%s.py'%(outdir,outmolecule),'w')
    write_pymol_script.write(Pymol_follow)          
    write_pymol_script.close()                   


    try:
        os.remove('%s-conformers.%s'%(outmolecule,outformat))
        os.remove('tmp.mmff94')
    except:
        null=0
    print '== ConfBuster Version 1.0'
    print "== Conformations are in the directory %s" %(molecule.rsplit('.',1)[0])
    print "== Enter pymol Follow-%s.py to visualize conformations with pymol" %(molecule.rsplit('.',1)[0])
    print '== Follow updates at github.com/patricklague/ConfBuster'

if __name__ == "__main__":
    main(sys.argv[1:])
