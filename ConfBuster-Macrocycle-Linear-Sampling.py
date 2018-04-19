#!/usr/bin/python
"""
    Macrocycle-linear-sampling.py
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

import os, sys, math, getopt
import networkx as nx

#
#   ADVANCED OPTION especially when smile from mol2 leads to NH+.
#
#   -x bolean (0,1) 1 uses xyz to calculate smile reference instead of mol2
#


def main(argv):

    PYMOL_PATH = '' # path to the pymol executable. Leave empty if in path 
    BABEL_PATH = '' # path to the babel executables. Leave empty if in path

    REPEAT      =   5    ### change with -n ### number of rotamer conformational searches to perform
    TMP2_NCONF  =   5    ### change with -N ### number of conformations generated per conformational search
                                            ### total number of generated conformations per cliving point = REPEAT * TMP2_NCONF  
    REMOVE_DOUBLE_BOND_CUT = True   # if false, single bonds between conjugated double bonds will be used a cliving points.
                                    # in some cases, this greatly increases the number of cliving points 

    RMSCUTOFF = 0.5      # rms deviation cutoff in A.


    # diheral angle values used for sampling 
    DH_VALUES = ['40','80','160','200','280','320']

    mandatory = ''
    in_name = ''
    threads = 2 #babel verion 2.4.1 does not seem to work with threads...
    out_dir_name = ''
    USE_XYZ_FOR_COMP = False 
    FF        = "MMFF94s" # The forceField used in Babel may be changed here
    rotate_conjugated = False 
    try:
        opts, args = getopt.getopt(argv,"hi:T:r:g:b:o:s:n:N:")
    except getopt.GetoptError:
        print 'ConfBuster-Macrocycle-Linear-Sampling.py -i inputfile [mandatory] -r rms deviation cutoff [0.5] -n for each cliving point, number of rotamer searches performed [5] -N for each cliving point, number of conformations extracted from each rotamer search [5] -o output directory name' 
        sys.exit(2)

    for opt, arg in opts:
        print "%s %s" %(opt, arg)
        if opt == '-h':
            print 'ConfBuster-Macrocycle-Linear-Sampling.py -i inputfile [mandatory] -r rms deviation cutoff [0.5] -n for each cliving point, number of rotamer searches performed [5] -N for each cliving point, number of conformations extracted from each rotamer search [5] -o output directory name' 
            sys.exit()
        elif opt in ("-i"):
            in_name = arg
            mandatory = '1'
        elif opt in ("-T"):
            threads = int(arg)
        elif opt in ("-r"):
            RMSCUTOFF = float(arg)
        elif opt in ("-n"):
            REPEAT = int(arg)
        elif opt in ("-N"):
            TMP2_NCONF = int(arg)
        elif opt in ("-o"):   
            out_dir_name = arg
        elif opt in ("-R"):
            if int(arg) in [0,1]:    
                if int(arg) == 1:
                    rotate_conjugated = True
                else:
                    rotate_conjugated = False
            else:
                sys.exit('-R values are 0 or 1')
        elif opt in ("-x"):
            if int(arg) in [0,1]:    
                if int(arg) == 1:
                    USE_XYZ_FOR_COMP = True
                else:
                    USE_XYZ_FOR_COMP = False
            else:
                sys.exit('-x values are 0 or 1')
                
    if mandatory != '1':
        print 'ConfBuster-Macrocycle-Linear-Sampling.py -i inputfile [mandatory] -r rms deviation cutoff [0.5] -n for each cliving point, number of rotamer searches performed [5] -N for each cliving point, number of conformations extracted from each rotamer search [5] -o output directory name' 
        sys.exit()

    # setting the number of threads
    os.environ['OMP_NUM_THREADS']="%d" %(threads)

    
    NUMBER = REPEAT*TMP2_NCONF


    print '%s conformations will be built for each cliving point'%(str(NUMBER))

    # out directory                           
    if out_dir_name != '':                    
        out_dir = out_dir_name+'/'            
        os.system('mkdir -p '+out_dir)        
    else:                                     
        out_dir = in_name.rsplit('.',1)[0]+'/'
        os.system('mkdir -p '+out_dir)        


    #
    #  write pymol RMS code
    #

    Pymol_RMS = """
dirtolist = '%s'

ls = os.listdir(dirtolist)

cmd.load('tmp.mol2')
refname = 'tmp'

out = open('tmp_rms.dat','w')
for x in range(len(ls)):
    if 'mol2' in ls[x]:
        cmd.load(dirtolist+ls[x])
        loaded = ls[x].rsplit('.',1)[0]
        rmsd = cmd.pair_fit(loaded,refname)
        out.write( loaded+'\t'+refname+'\t'+str(rmsd)+'\\n')
out.close()"""%(out_dir)

    write_pymol_script = open('Pymol-RMS.py','w')
    write_pymol_script.write(Pymol_RMS)
    write_pymol_script.close()

    #
    #   write pymol follow code
    #

    pymol_follow = """
cmd.reinitialize()
ls = os.listdir('.')

load_order = []

for x in ls:
    if 'mol2' in x:
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
cmd.color('black')
cmd.bg_color('white')
cmd.set('line_width','0.5')
cmd.do('util.cnc')
cmd.color('green','elem c and '+refname)
cmd.show('sticks',refname)
cmd.set('stick_radius','0.1')
"""

    out = open(out_dir+'/Follow-%s.py'%(in_name.split('.mol2')[0]),'w')
    out.write(pymol_follow)
    out.close()


    #
    #   load molecule
    #



    f_in = open(in_name,'r')
    r_in = f_in.read()
    f_in.close()
    s_in = r_in.splitlines()

    #
    # make smile ref
    #

    if USE_XYZ_FOR_COMP == True:
        make_ref = '%sobabel %s -O ref.xyz 2> /dev/null'%(BABEL_PATH, in_name)
        os.system(make_ref)
        make_ref = '%sobprop ref.xyz > ref.prop 2> /dev/null'%(BABEL_PATH)
        os.system(make_ref)
    else:
        make_ref = '%sobprop %s > ref.prop 2> /dev/null'%( BABEL_PATH, in_name)
        os.system(make_ref)
    f_ref = open('ref.prop','r')
    r_ref = f_ref.read()
    f_ref.close()
    s_ref = r_ref.splitlines()
    for x in s_ref:
        if 'canonical_SMILES' in x:
            ref_smi = x.split()[1]

    #print ref_smi

    #
    #   get header
    #

    head_list = []

    start_head = 0
    end_head = 0

    for x in range(len(s_in)):
        if start_head != 0:
            if len(s_in[x].split()) == 0:
                if end_head == 0:
                    end_head = x

        if '@<TRIPOS>MOLECULE' in s_in[x]:
            start_head = x+1

    print 'heads'
    print start_head, end_head

    for x in range(start_head, end_head):
        head_list.append(s_in[x])


    #
    #   get atom list
    #


    atom_list = []

    start_atom = 0
    end_atom = 0

    for x in range(len(s_in)):
        if start_atom != 0:
            if len(s_in[x].split()) != 9:
                if end_atom == 0:
                    end_atom = x

        if '@<TRIPOS>ATOM' in s_in[x]:
            start_atom = x+1

    print 'atoms'
    print start_atom, end_atom

    for x in range(start_atom, end_atom):
        atom_list.append(s_in[x])

    #
    #  get bonds
    #


    bond_list = []

    start_bond = 0
    end_bond = len(s_in)

    for x in range(len(s_in)):

        if start_bond != 0:
            if len(s_in[x].split()) != 4:
                end_bond = x

        if '@<TRIPOS>BOND' in s_in[x]:
            start_bond = x+1

    print 'bonds'
    print start_bond, end_bond

    for x in range(start_bond, end_bond):
        bond_list.append(s_in[x])



    #
    #   find main cycle
    #

    G=nx.Graph()

    for x in range( start_bond,end_bond ):
        sp = s_in[x].split()
        a1 = sp[1]
        a2 = sp[2]
        bond_type = sp[3]
        G.add_edge(a1, a2, bonding=bond_type)
    cycle_list_a = (nx.cycle_basis(G) )

    cycle_list = sorted(cycle_list_a, key = len, reverse=True)

    if len(cycle_list) == 0:
        sys.exit('This is not a macrocyclic molecule... Not even a cyclic one ! Have Fun.')
    else:
        if len(cycle_list[0]) <= 6:
            sys.exit('This does not seem to be macrocyclic molecule... ! Have Fun.')
            

    small_cycle_list = []
        

    for x in cycle_list:
        if len(x) <= 6:
            for y in x:
                small_cycle_list.append(y)
    


    #
    #   exclude list
    #

    #   exclude atoms involved in double bonds



    ######
    #
    #    exclusion START
    #
    #####


    exclude_list = []

    # commentout for extended search ... does generate much more amount of not accepted results but you might get lucky

    if REMOVE_DOUBLE_BOND_CUT == True:
        for x in bond_list:
            sp = x.split()
            bond_type = sp[3]
            if bond_type == '2':
                is_a_c = 0
                for y in atom_list:
                    spp = y.split()
                    if sp[1] == spp[0] or sp[2] == spp[0]:
                        if spp[1][0] == 'C':
                            is_a_c += 1
                if is_a_c == 2:
                    exclude_list.append(sp[1])
                    exclude_list.append(sp[2])


    ######
    #
    #    exclusion  END
    #
    #####



    #
    #   obchiral
    #

    command = '%sobchiral %s > tmp.chiral'%(BABEL_PATH, in_name)
    os.system(command)

    f_chiral = open('tmp.chiral','r')
    r_chiral = f_chiral.read()
    f_chiral.close()
    s_chiral = r_chiral.splitlines()

    s_chiral.pop(0) # remove title

    chiral_list = []

    for x in s_chiral:
        atom = x.split()[1]
        chiral_list.append(atom)

    #
    #   DEFINITIONS
    #


    def write_cut_bond_pymol(name,a1,a2):
        pymol_cut_bond = "cmd.load('%s')\n"%(name)
        pymol_cut_bond += "cmd.unbond('id %s','id %s')\n"%(a1,a2)
        pymol_cut_bond += "cmd.h_add('id %s')\n"%(a1)
        pymol_cut_bond += "cmd.extract('h1','id %s')\n"%(str(len(atom_list)+1))
        pymol_cut_bond += "cmd.h_add('id %s')\n"%(a2)
        pymol_cut_bond += "cmd.extract('h2','id %s')\n"%(str(len(atom_list)+2))
        pymol_cut_bond += "cmd.save('tmp_h1.pdb','h1')\n"
        pymol_cut_bond += "cmd.save('tmp_h2.pdb','h2')\n"
        pymol_cut_bond += "cmd.quit()\n"
        out = open('tmp_pymol_cut_bond.py','w')
        out.write(pymol_cut_bond)
        out.close()


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

    def get_E(name):
        tmp_f = open(name,'r')
        tmp_r = tmp_f.read()
        tmp_f.close()
        tmp_s = tmp_r.splitlines()
        for x in tmp_s:
            if 'TOTAL ENERGY =' in x:
                energy = x.split()[3]
                return energy

    #
    #   identify cliving position
    #
    print '== Cliving positions:'
    for x in range( start_bond,end_bond ):
        sp = s_in[x].split()
        a1 = sp[1]
        a2 = sp[2]
        bond_type = sp[3]
        if bond_type in ['1','am'] \
            and a1 not in chiral_list \
            and a2 not in chiral_list \
            and atom_list[int(a1)-1].split()[1] != 'H' \
            and atom_list[int(a2)-1].split()[1] != 'H' \
            and atom_list[int(a1)-1].split()[1][0] != 'S' \
            and atom_list[int(a2)-1].split()[1][0] != 'S' \
            and a1 not in exclude_list \
            and a2 not in exclude_list \
            and a1 in cycle_list[0] \
            and a2 in cycle_list[0] \
            and a1 not in small_cycle_list\
            and a2 not in small_cycle_list:

            print a1,a2


    print '== Starting with position:'
    for x in range( start_bond,end_bond ):
        sp = s_in[x].split()
        a1 = sp[1]
        a2 = sp[2]
        bond_type = sp[3]
        if bond_type in ['1','am'] \
            and a1 not in chiral_list \
            and a2 not in chiral_list \
            and atom_list[int(a1)-1].split()[1] != 'H' \
            and atom_list[int(a2)-1].split()[1] != 'H' \
            and atom_list[int(a1)-1].split()[1][0] != 'S' \
            and atom_list[int(a2)-1].split()[1][0] != 'S' \
            and a1 not in exclude_list \
            and a2 not in exclude_list \
            and a1 in cycle_list[0] \
            and a2 in cycle_list[0] \
            and a1 not in small_cycle_list\
            and a2 not in small_cycle_list:

            print 'using positions',a1,a2


            command = '%sobabel %s -O tmp.mol2 2> /dev/null'%(BABEL_PATH,in_name)
            os.system(command)

            #
            #   this uses pymol to cut the desired bond and add h. is saves the two added h individually.
            #


            write_cut_bond_pymol('tmp.mol2',a1,a2)
    
            os.system('%spymol -cQ tmp_pymol_cut_bond.py'%(PYMOL_PATH))

            

            #
            #   convert hs in mol2
            #

            command = '%sobabel tmp_h1.pdb -O tmp_h1.mol2'%(BABEL_PATH)
            os.system(command)
            command = '%sobabel tmp_h2.pdb -O tmp_h2.mol2'%(BABEL_PATH)
            os.system(command)

            out = open('tmp.mol2','w')
            out.write('@<TRIPOS>MOLECULE\n')
            for a in head_list:
                sp = a.split()
                if len(sp) == 5:
                    new_line = ' '+str(int(sp[0])+2)+' '+str(int(sp[1])+1)+' '+sp[2]+' '+sp[3]+' '+sp[4]+'\n'
                    out.write(new_line)
                else:
                    out.write(a+'\n')

            out.write('\n@<TRIPOS>ATOM\n')
            for a in atom_list:
                out.write(a+'\n')

            f_h1 = open('tmp_h1.mol2','r')
            r_h1 = f_h1.read()
            f_h1.close()
            s_h1 = r_h1.splitlines()
            for a in s_h1:
                sp = a.split()
                if len(sp) == 9:
                    new_line = ' '*(7 - len(str(len(atom_list)+1))) + str(len(atom_list)+1) + a[8:] + '\n'
                    out.write(new_line)

            f_h2 = open('tmp_h2.mol2','r')
            r_h2 = f_h2.read()
            f_h2.close()
            s_h2 = r_h2.splitlines()
            for a in s_h2:
                sp = a.split()
                if len(sp) == 9:
                    new_line = ' '*(7 - len(str(len(atom_list)+2))) + str(len(atom_list)+2) + a[8:] + '\n'
                    out.write(new_line)

            out.write('@<TRIPOS>BOND\n')
            new_bond_list = []
            for a in bond_list:
                sp = a.split()
                if sp[1] in [a1,a2] and sp[2] in [a1,a2]:
                    null = 0
                else:
                    new_bond_list.append([sp[1],sp[2],sp[3]])

            new_bond_list.append([a1, str(len(atom_list)+1), '1'])
            new_bond_list.append([a2, str(len(atom_list)+2), '1'])

            for a in range(len(new_bond_list)):
                line  = ' '*(6 - len(str(a+1))) + str(a+1)
                line += ' '*(6 - len(str(new_bond_list[a][0]))) + str(new_bond_list[a][0])
                line += ' '*(6 - len(str(new_bond_list[a][1]))) + str(new_bond_list[a][1])
                line += ' '*(6 - len(str(new_bond_list[a][2]))) + str(new_bond_list[a][2])
                line += '\n'
                out.write(line)
            out.close()
            


            print '== Minimization of the molecule'
            command = 'ConfBuster-Single-Molecule-Minimization.py -i tmp.mol2 -F %s -T %s 2> /dev/null'%(FF,threads)
            os.system(command)
    

            #
            #   Make rotamers for current clived bond
            #
            try:
                os.system('rm -r tmp2_out_conf/')
            except:
                null = 0
            os.system('mkdir tmp2_out_conf/')

            for x in range(REPEAT):
                if os.path.exists("tmp_out_conf"):
                    os.system('rm -r tmp_out_conf/')

                print '== Looking for rotamers'
                rotamers = '50' #the number of rotamers used for the rotamer search may be changed here
                command = 'ConfBuster-Rotamer-Search.py -i tmp.mol2 -g %s -G energy -e 1000 -o tmp_out_conf -f xyz -T %s '%(rotamers, threads)
                print command
                os.system(command)
    
    
                dir1      = 'tmp_out_conf/'
    
                atom1 = int(a1)
                atom2 = int(a2)
                remove1 = len(atom_list) + 1
                remove2 = len(atom_list) + 2
    
                #
                #   Remake cycle and minimize
                #
    
    
                #
                # find rotatable bonds in the linearized cycle 
                # including bonds between double bonds. (are excluded of conservec at the end of this segment.) 
                # these bonds are rotated systematiiclay with pymol to minimize de distance between atoms belonging to the
                # bond which was cut in order to linearize.
    
                r_bond_l = []            
    
                for a in new_bond_list:
                    aa1 = a[0]
                    aa2 = a[1]
                    part_of_big_macro = 0
                    count_macro = 0
                    for xxx in range(len(cycle_list)):
                        if len(cycle_list[xxx]) > 6:
                            count_macro += 1
                            if aa1 in cycle_list[xxx] and aa2 in cycle_list[xxx]:
                                part_of_big_macro += 1
                    if part_of_big_macro == count_macro:
                        isH = 0
                        for y in atom_list:
                            sp = y.split()
                            if aa1 == sp[0] and sp[1][0] == 'H':
                                isH += 1 
                        if isH == 0:
                            cycle_check = 0
                            for x in cycle_list:
                                if len(x) <= 7:
                                    if aa1 in x and aa2 in x:
                                        cycle_check += 1
                            if cycle_check <= 1:
                                r_bond_l.append([a[0],a[1]])
    
    
                print
                print
                print

                angles_l = []
    
                for x in r_bond_l:
                    done = 0
                    atom = x[0]
                    for y in new_bond_list:
                        if done == 0:
                            if atom == y[0] and y[1] not in x:
                                for z in atom_list:
                                    sp = z.split()
                                    if y[1] == sp[0] and sp[1][0] != 'H':
                                        angles_l.append([y[1]]+x)
                                        done = 1
                                        break
                            elif atom == y[1] and y[0] not in x:
                                for z in atom_list:
                                    sp = z.split()
                                    if y[0] == sp[0] and sp[1][0] != 'H':
                                        angles_l.append([y[0]]+x)
                                        done = 1
                                        break

    
                dh_l = []
                for x in angles_l:
                    done = 0
                    atom = x[2]
                    for y in new_bond_list:
                        if done == 0:
                            if atom == y[0] and y[1] not in x:
                                for z in atom_list:
                                    sp = z.split()
                                    if y[1] == sp[0] and sp[1][0] != 'H':
                                        dh_l.append(x+[y[1]])
                                        done = 1
                                        break
                            elif atom == y[1] and y[0] not in x:
                                for z in atom_list:
                                    sp = z.split()
                                    if y[0] == sp[0] and sp[1][0] != 'H':
                                        dh_l.append(x+[y[0]])
                                        done = 1
                                        break
    
    
                dh_l_2 = []
                for x in dh_l:
                    btype0 = ''
                    btype1 = ''
                    btype2 = ''
                    for y in new_bond_list:
                        if x[1] in y and x[2] in y:
                            btype0 = y[2]
                    for y in new_bond_list:
                        if x[0] in y and x[1] in y:
                            btype1 = y[2]
                    for y in new_bond_list:
                        if x[2] in y and x[3] in y:
                            btype2 = y[2]
                    if btype0 == '1':
                        if btype1 == '1' or btype2 == '1':
                            dh_l_2.append(x)
    
                # conserve or exclude bond between double bonds as rotatable. This should be excluded as usually, conjugated
                # double bonds will remain in a minimal energy conformation. FF are bad at describing the conjugation energy
                # thus it is better to keep fixed. (hopefully they are well placed after the obabe --rotamer.)

                if rotate_conjugated == True:
                    dh_l_2 = dh_l 


                ter_bond_w_h = []

        
                for x in new_bond_list:
                    if x[1] == str(len(atom_list)+1) or x[1] == str(len(atom_list)+2):
                        ter_bond_w_h.append([x[0],x[1]])

                ter_angle_w_h = []

                for y in ter_bond_w_h:
                    for x in new_bond_list:
                        if x[1] == y[0]:
                            if x[0] not in y:
                                for z in atom_list:
                                    sp = z.split()
                                    if x[0] == sp[0] and sp[1][0] != 'H':
                                        ter_angle_w_h.append([sp[0],y[0],y[1]])
                        if x[0] == y[0]:
                            if x[1] not in y:
                                for z in atom_list:
                                    sp = z.split()
                                    if x[1] == sp[0] and sp[1][0] != 'H':
                                        ter_angle_w_h.append([sp[0],y[0],y[1]])

                ter_dh_w_h = []
                
                for y in ter_angle_w_h:
                    done = 0
                    for x in new_bond_list:
                        if done == 0:
                            if x[1] == y[0]:
                                if x[0] not in y:
                                    for z in atom_list:
                                        sp = z.split()
                                        if x[0] == sp[0] and sp[1][0] != 'H':
                                            ter_dh_w_h.append([sp[0],y[0],y[1],y[2]])
                                            done += 1
                                            break
                            if x[0] == y[0]:
                                if x[1] not in y:
                                    for z in atom_list:
                                        sp = z.split()
                                        if x[1] == sp[0] and sp[1][0] != 'H':
                                            ter_dh_w_h.append([sp[0],y[0],y[1],y[2]])
                                            done += 1
                                            break

                for x in ter_dh_w_h:
                    dh_l_2.append(x) 
                #print ter_bond_w_h
                #print
                #print ter_angle_w_h
                #print
                #print ter_dh_w_h


                ##for debuging 
    
                #print 'r bonds'
                #for x in r_bond_l:
                #    print x
        
                #print
    
                #print 'r angles'
                #for x in angles_l:
                #    print x
        
                #print
    
                #print 'r dh'
                #for x in dh_l:
                #    print x
    
                #print 'r dh 2'
                #for x in dh_l_2:
                #    print x
    

                ls = os.listdir(dir1)
                #dist_d = {}
                #d_l = []
                to_open_l = []
    
                to_sort = []
                for x in ls:
                    if 'mol2' in x or 'xyz' in x:
                        rank = int(x.split('.')[0].split('-c')[1])  
                        to_sort.append([rank,x])
                to_sort.sort()
                
    
                # change the range to use more than only the best conformation obtained from the obabel --rotamer   
                for x in range(1):
                    to_open_l.append(dir1+to_sort[x][1])

                for x in to_open_l:
                                
                    pymol_dihedral_opt = """
from random import shuffle
    
dihedral_list_one = ["""
    
                    for n in dh_l_2:
                        pymol_dihedral_opt += "("
                        for nn in n:
                            pymol_dihedral_opt += "'"+nn+"'," 
                        pymol_dihedral_opt = pymol_dihedral_opt[:-1]
                        pymol_dihedral_opt += "),"
    
                    pymol_dihedral_opt = pymol_dihedral_opt[:-1]
                    pymol_dihedral_opt += """]
#dihedral_values = ['30','60','90','150','180','210','270','300','330']
dihedral_values = %s 
a1 = '%s'
a2 = '%s'
a3 = '%s'
a4 = '%s'
count = 0

final_dh = []
final_dh_w_H = []
dihedral_list = []

for x in dihedral_list_one:
    if a1 in x or a2 in x:
        if a3 not in x and a4 not in x:
            final_dh.append(x)
        else:
            final_dh_w_H.append(x)
    else:
        dihedral_list.append(x)

#print final_dh
print
#print final_dh_w_H


cmd.reinitialize()
cmd.load('%s','tmp')
                    """%(str(DH_VALUES),str(a1),str(a2),str(len(atom_list)+1),str(len(atom_list)+2),x)
    
                    pymol_dihedral_opt += """
counting = 0
pass_list = []
for dangle in dihedral_list:
    for value in dihedral_values:
        counting += 1
        if counting %  1000 == 0:
            print counting
        ini_dh = cmd.get_dihedral('id %s'%(dangle[0]) ,'id %s'%(dangle[1]),'id %s'%(dangle[2]),'id %s'%(dangle[3]))
        cmd.set_dihedral('id %s'%(dangle[0]) ,'id %s'%(dangle[1]),'id %s'%(dangle[2]),'id %s'%(dangle[3]), '%s'%(value))
        d1 = float(cmd.distance('id %s'%(a1),'id %s'%(a2)))
        pass_list.append([d1,value,(dangle[0],dangle[1],dangle[2],dangle[3])])
        cmd.set_dihedral('id %s'%(dangle[0]) ,'id %s'%(dangle[1]),'id %s'%(dangle[2]),'id %s'%(dangle[3]), '%s'%(ini_dh))

for dangle in dihedral_list:
    for value in dihedral_values:
        for dangle2 in dihedral_list:
            for value2 in dihedral_values:"""
                    pymol_dihedral_opt += """
                cmd.reinitialize()
                cmd.load('%s','tmp')"""%(x)
                    pymol_dihedral_opt += """
                counting += 1
                if counting %  1000 == 0:
                    print counting
                cmd.set_dihedral('id %s'%(dangle[0]) ,'id %s'%(dangle[1]),'id %s'%(dangle[2]),'id %s'%(dangle[3]), '%s'%(value))
                cmd.set_dihedral('id %s'%(dangle2[0]) ,'id %s'%(dangle2[1]),'id %s'%(dangle2[2]),'id %s'%(dangle2[3]), '%s'%(value2))
                d1 = float(cmd.distance('id %s'%(a1),'id %s'%(a2)))
                pass_list.append([d1, value, (dangle[0],dangle[1],dangle[2],dangle[3]), value2, (dangle2[0],dangle2[1],dangle2[2],dangle2[3]) ])



pass_list.sort()"""

                    pymol_dihedral_opt += """
get_mol = 0
count_mol = 0
while get_mol < %s:
    cmd.reinitialize()
    cmd.load('%s','tmp')"""%(str(TMP2_NCONF),x)

                    pymol_dihedral_opt += """
    print count_mol, get_mol
    x = count_mol
    try :
        if len(pass_list[x]) == 3:
            cmd.set_dihedral('id %s'%(pass_list[x][2][0]) ,'id %s'%(pass_list[x][2][1]),'id %s'%(pass_list[x][2][2]),'id %s'%(pass_list[x][2][3]), '%s'%(pass_list[x][1]))
        if len(pass_list[x]) == 5:
            cmd.set_dihedral('id %s'%(pass_list[x][2][0]) ,'id %s'%(pass_list[x][2][1]),'id %s'%(pass_list[x][2][2]),'id %s'%(pass_list[x][2][3]), '%s'%(pass_list[x][1]))
            cmd.set_dihedral('id %s'%(pass_list[x][4][0]) ,'id %s'%(pass_list[x][4][1]),'id %s'%(pass_list[x][4][2]),'id %s'%(pass_list[x][4][3]), '%s'%(pass_list[x][3]))
        if len(pass_list[x]) == 7:
            cmd.set_dihedral('id %s'%(pass_list[x][2][0]) ,'id %s'%(pass_list[x][2][1]),'id %s'%(pass_list[x][2][2]),'id %s'%(pass_list[x][2][3]), '%s'%(pass_list[x][1]))
            cmd.set_dihedral('id %s'%(pass_list[x][4][0]) ,'id %s'%(pass_list[x][4][1]),'id %s'%(pass_list[x][4][2]),'id %s'%(pass_list[x][4][3]), '%s'%(pass_list[x][3]))
            cmd.set_dihedral('id %s'%(pass_list[x][6][0]) ,'id %s'%(pass_list[x][6][1]),'id %s'%(pass_list[x][6][2]),'id %s'%(pass_list[x][6][3]), '%s'%(pass_list[x][5]))

        agin_dh_list = []
        dh_values = ['15','30','45','60','75','90','105','135','150','165','180','195','210','225','255','270','285','300','315','330','345']
        #dh_values = ['40','60','80','160','180','200','280','300','320']
        for fdh in dihedral_list:
            ini_value = cmd.get_dihedral('id %s'%(fdh[0]) ,'id %s'%(fdh[1]),'id %s'%(fdh[2]),'id %s'%(fdh[3]))
            for value in dh_values:
                cmd.set_dihedral('id %s'%(fdh[0]) ,'id %s'%(fdh[1]),'id %s'%(fdh[2]),'id %s'%(fdh[3]), '%s'%(value))
                d1 = float(cmd.distance('id %s'%(a1),'id %s'%(a2)))
                if 1.5 < d1 :
                    agin_dh_list.append([d1,['%s'%(fdh[0]) ,'%s'%(fdh[1]),'%s'%(fdh[2]),'%s'%(fdh[3]), '%s'%(value)]])
            cmd.set_dihedral('id %s'%(fdh[0]) ,'id %s'%(fdh[1]),'id %s'%(fdh[2]),'id %s'%(fdh[3]), '%s'%(ini_value))

        agin_dh_list.sort()
        cmd.set_dihedral('id %s'%(agin_dh_list[0][1][0]) ,'id %s'%(agin_dh_list[0][1][1]),'id %s'%(agin_dh_list[0][1][2]),'id %s'%(agin_dh_list[0][1][3]), '%s'%(agin_dh_list[0][1][4]))


        clash = 0"""
                    pymol_dihedral_opt += """
        sel_one   = cmd.select('selone','elem H')
        sel_check = cmd.select('check','selone around 0.9 ')
        if sel_check > 0:
            clash += 1
        if clash == 0:

            
            dh_values = ['0','30','60','90','120','150','180','210','240','270','300','330','360']
            for fdh in final_dh:
                best_d = float(cmd.distance('id %s'%(a1),'id %s'%(a2)))
                best_value = cmd.get_dihedral('id %s'%(fdh[0]) ,'id %s'%(fdh[1]),'id %s'%(fdh[2]),'id %s'%(fdh[3]))
                for value in dh_values:
                    counting += 1
                    cmd.set_dihedral('id %s'%(fdh[0]) ,'id %s'%(fdh[1]),'id %s'%(fdh[2]),'id %s'%(fdh[3]), '%s'%(value))
                    d1 = float(cmd.distance('id %s'%(a1),'id %s'%(a2)))
                    if d1 < best_d :
                        best_d = d1
                        best_value = value
                cmd.set_dihedral('id %s'%(fdh[0]) ,'id %s'%(fdh[1]),'id %s'%(fdh[2]),'id %s'%(fdh[3]), '%s'%(best_value))

            dh_values = ['0','15','30','45','60','75','90','105','120','135','150','165','180','195','210','225','240','255','270','285','300','315','330','345','360']
            for fdh in final_dh_w_H:
                best_d = float(cmd.distance('id %s'%(a3),'id %s'%(a4)))
                best_value = cmd.get_dihedral('id %s'%(fdh[0]) ,'id %s'%(fdh[1]),'id %s'%(fdh[2]),'id %s'%(fdh[3]))
                for value in dh_values:
                    counting += 1
                    cmd.set_dihedral('id %s'%(fdh[0]) ,'id %s'%(fdh[1]),'id %s'%(fdh[2]),'id %s'%(fdh[3]), '%s'%(value))
                    d1 = float(cmd.distance('id %s'%(a3),'id %s'%(a4)))
                    if d1 < best_d :
                        best_d = d1
                        best_value = value
                cmd.set_dihedral('id %s'%(fdh[0]) ,'id %s'%(fdh[1]),'id %s'%(fdh[2]),'id %s'%(fdh[3]), '%s'%(best_value))


            lsx = os.listdir('tmp2_out_conf')

            is_already_there = 0
            for conf in lsx:
                if 'mol2' in conf or 'xyz' in conf or 'pdb' in conf:
                    cmd.load('tmp2_out_conf/%s'%(conf),'there')
                    rms = cmd.pair_fit('there','tmp')"""
                    pymol_dihedral_opt += """
                    if rms < %s:"""%(RMSCUTOFF)
                    pymol_dihedral_opt += """
                        is_already_there += 1
                    cmd.delete('there')
            if is_already_there == 0:
                cmd.save('tmp2_out_conf/tmp2-%s.pdb'%(str(len(lsx))))
                get_mol += 1
        count_mol += 1
    except:
        break
                """
                    write_pymol_script = open('tmp_pymol_dihedral_opt.py','w')
                    write_pymol_script.write(pymol_dihedral_opt)
                    write_pymol_script.close()
                    print 'building pool of %s linear conformations for cyclisation. %s of %s done'%(str(NUMBER), str(len(os.listdir('tmp2_out_conf/'))), str(NUMBER) )
                    pymol_dh_opt = '%spymol -cQ tmp_pymol_dihedral_opt.py'%(PYMOL_PATH)
                    os.system(pymol_dh_opt)
    


            count = 0
         
            to_open_new_l = []
            for x in os.listdir('tmp2_out_conf/'):
                if '.mol2' in x or '.xyz' in x or '.pdb' in x:
                    to_open_new_l.append(x)
        
            for x in to_open_new_l:

                count += 1
                print 'Conformation %s of %s'%(str(count),str(len(to_open_new_l)))

                name = x.rsplit('.',1)[0]
                print name
                command = '%sobabel %s -O tmp.mol2'%(BABEL_PATH, 'tmp2_out_conf/'+x)
                os.system(command)


                # remove last two hydrogenes to make bond

                f_in_tmp = open('tmp.mol2','r')
                r_in_tmp = f_in_tmp.read()
                f_in_tmp.close()
                s_in_tmp = r_in_tmp.splitlines()

                #
                #   get header
                #

                head_list_tmp = []

                start_head_tmp = 0
                end_head_tmp = 0

                for x in range(len(s_in_tmp)):
                    if start_head_tmp != 0:
                        if len(s_in_tmp[x].split()) == 0:
                            if end_head_tmp == 0:
                                end_head_tmp = x

                    if '@<TRIPOS>MOLECULE' in s_in_tmp[x]:
                        start_head_tmp = x+1

                #print 'heads'
                #print start_head_tmp, end_head_tmp

                for x in range(start_head_tmp, end_head_tmp):
                    head_list_tmp.append(s_in_tmp[x])


                #
                #   get atom list
                #


                atom_list_tmp = []

                start_atom_tmp = 0
                end_atom_tmp = 0

                for x in range(len(s_in_tmp)):
                    if start_atom_tmp != 0:
                        if len(s_in_tmp[x].split()) != 9:
                            if end_atom_tmp == 0:
                                end_atom_tmp = x

                    if '@<TRIPOS>ATOM' in s_in_tmp[x]:
                        start_atom_tmp = x+1

                #print 'atoms'
                #print start_atom_tmp, end_atom_tmp

                for x in range(start_atom_tmp, end_atom_tmp):
                    atom_list_tmp.append(s_in_tmp[x])

                #
                #  get bonds
                #


                bond_list_tmp = []

                start_bond_tmp = 0
                end_bond_tmp = len(s_in_tmp)

                for x in range(len(s_in_tmp)):

                    if start_bond_tmp != 0:
                        if len(s_in_tmp[x].split()) != 4:
                            end_bond_tmp = x

                    if '@<TRIPOS>BOND' in s_in_tmp[x]:
                        start_bond_tmp = x+1

                #print 'bonds'
                #print start_bond_tmp, end_bond_tmp

                for x in range(start_bond_tmp, end_bond_tmp):
                    bond_list_tmp.append(s_in_tmp[x])

                out = open('tmp.mol2','w')
                out.write('@<TRIPOS>MOLECULE\n')
                for a in head_list_tmp:
                    sp = a.split()
                    if len(sp) == 5:
                        new_line = ' '+str(len(atom_list))+' '+str(len(bond_list))+' '+sp[2]+' '+sp[3]+' '+sp[4]+'\n'
                        out.write(new_line)
                    else:
                        out.write(a+'\n')

                out.write('\n@<TRIPOS>ATOM\n')
                for a in range(len(atom_list)):
                    out.write(atom_list_tmp[a]+'\n')

                out.write('@<TRIPOS>BOND\n')
                new_bond_list = []
                for a in bond_list:
                    sp = a.split()
                    new_bond_list.append([sp[1],sp[2],sp[3]])

                new_bond_list.append([a1, a2, 1])

                for a in range(len(new_bond_list)):
                    line  = ' '*(6 - len(str(a+1))) + str(a+1)
                    line += ' '*(6 - len(str(new_bond_list[a][0]))) + str(new_bond_list[a][0])
                    line += ' '*(6 - len(str(new_bond_list[a][1]))) + str(new_bond_list[a][1])
                    line += ' '*(6 - len(str(new_bond_list[a][2]))) + str(new_bond_list[a][2])
                    line += '\n'
                    out.write(line)
                out.close()


                command = '%sobminimize -ff %s -c 1e-8 -sd -n 5000 tmp.mol2 1> tmp1.pdb 2> /dev/null'%(BABEL_PATH,FF)
                print command
                os.system(command)

                # remove WARNING from xyz
                #remove_warning_xyz('tmp1.xyz')


                # do the rotamer here bellow twice if need be, once with energy and second with rmsd 
                for rotamer_keyword in ['energy']:

                    command = '%sobabel tmp1.pdb -O tmp.mol2 --conformer --nconf 50 --score %s 1> tmp 2> /dev/null'%(BABEL_PATH, rotamer_keyword)
                    print command
                    os.system(command)

                    command = '%sobminimize -ff %s -c 1e-8 -sd -n 5000 tmp.mol2 1> tmp1.pdb 2> /dev/null'%(BABEL_PATH,FF)
                    print command
                    os.system(command)

                    command = '%sobabel tmp1.pdb -O tmp1.xyz'%(BABEL_PATH)
                    print command
                    os.system(command)
                    command = '%sobabel tmp1.xyz -O tmp.mol2'%(BABEL_PATH)
                    print command
                    os.system(command)
    
                    command = '%sobenergy -ff %s tmp.mol2 > tmp.ff'%(BABEL_PATH,FF)
                    print command
                    os.system(command)
    
                    # read energy from file and put as xyz title
    
                    energy = get_E('tmp.ff')
                    if energy == None:
                        sys.exit('energy is none')
                    print 'energy :',energy
    
                    f = open('tmp.mol2','r')
                    r = f.read()
                    f.close()
    
                    if len(r) > 0:
    
                        s = r.splitlines()
    
                        natoms = s.pop(0)
                        title  = s.pop(0)
    
                        out = open('tmp.mol2','w')
                        out.write(natoms+'\n')
                        out.write(str(energy)+'\n')
    
                        for y in s:
                            out.write(y+'\n')
    
                        out.close()
    
    
                        # make smile to compare
                        if USE_XYZ_FOR_COMP == True:
                            os.system('%sobabel tmp.mol2 -O tmp.xyz'%(BABEL_PATH))
                            make_check = '%sobprop tmp.xyz > tmp.prop 2> /dev/null'%(BABEL_PATH)
                            os.system(make_check)
                        else:
                            make_check = '%sobprop tmp.mol2 > tmp.prop 2> /dev/null'%(BABEL_PATH)
                            os.system(make_check)
    
                        fsmi = open('tmp.prop','r')
                        rsmi = fsmi.read()
                        fsmi.close()
                        ssmi = rsmi.splitlines()
                        for a in ssmi:
                            if 'canonical_SMILES' in a:
                                tmp_smi = a.split()[1]
                        print
                        print 'ref'
                        print ref_smi
                        print tmp_smi
                        print 'tmp'
                        print
                        if tmp_smi == ref_smi:
                            print 'smile comparison passed'
                            smile_cmp = 1
    
                            #command = '%sobabel tmp.xyz -O tmp.mol2'%(BABEL_PATH)
                            #os.system(command)
    
                            #
                            #   rmsd sorting
                            #
    
                            rmsd_pymol = '%spymol -cQ Pymol-RMS.py'%(PYMOL_PATH)
                            os.system(rmsd_pymol)
    
                            #   Read tmp_rms.dat
                            #   If rms is below rms cutoff, keep the conformatmion
                            f_rms = open('tmp_rms.dat','r')
                            r_rms = f_rms.read()
                            f_rms.close()
                            s_rms = r_rms.splitlines()
    
                            rms_check = []
                            for line in s_rms:
                                sp = line.split('\t')
                                line_rms = sp[2]
                                if  float(sp[2]) < RMSCUTOFF:
                                    rms_check.append([float(sp[2]),sp[0],sp[1]])
                            rms_check.sort()
    
    
                            if len(rms_check) == 0: # conformation is new and kept
                                print 'rms_check passed'
                                print 'keeping conformation'
                                ls1 = os.listdir(out_dir) # conformation is numbered as the number of files in out_dir + 1
                                command = 'cp tmp.mol2 %sconf-%s.mol2'%(out_dir,str(len(ls1)))
                                print command
                                os.system(command)
                            if len(rms_check) > 0: # compare energies with similar (rms < RMSCUTOFF), keep the best
                                print 'found similar conf'
                                print 'rms, conf, tmp'
                                for line in rms_check:
                                    if line[0] < RMSCUTOFF:
                                        print line
                                print
                                for conformations in rms_check:
                                    e_f = open(out_dir+conformations[1]+'.mol2')
                                    e_r = e_f.read()
                                    e_f.close()
                                    e_s = e_r.splitlines()
                                    kept_conf_E = float(e_s[1])
                                    if float(energy) < float(kept_conf_E):
                                        os.system('mv tmp.mol2 %s.mol2'%(out_dir+conformations[1]))
                                        print '%s E = %s'%(conformations[1],kept_conf_E)
                                        print 'current conf E = %s'%(str(energy))
                                        print 'replacing conformation %s with current one'%(conformations[1])
                                        break
                                    else:
                                        print 'kept conf %s is better. Checking with next similar conf if the case.'%(conformations[1])
    
    
    
    
                        else:
                            print 'smile comparison failed'
    
                        print '##############################'
                        print
                        print



    #
    #   remove tmp files
    #

    for rm in ['tmp*','Pymol-RMS.py','ref.prop','used-macro-start.txt','ref.xyz']:
        try:
            os.system("rm -rf %s"%(rm))
        except:
            null = 0
    print '== ConfBuster Version 1.0'
    print "== Conformations are in the directory %s" %(in_name.split('.mol2')[0])
    print "== Enter pymol Follow-%s.py to visualize conformations with pymol" %(in_name.split('.mol2')[0])
    print '== Follow updates at github.com/patricklague/ConfBuster'
if __name__ == "__main__":
    main(sys.argv[1:])
