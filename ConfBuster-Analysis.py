#!/usr/bin/python


import os, sys, math, random, getopt
import filecmp

"""
    ConfBuster-Analysis.py
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
    #
    #   START Matrix Colors
    #

    PYMOL_PATH = '' # path to the pymol executable. Leave empty if in path 
    BABEL_PATH = '' # path to the babel executables. Leave empty if in path
    
    
    COLOR_MAT_MIN   =   '#0270ff'
    COLOR_MAT_HALF  =   'white'
    COLOR_MAT_MAX   =   '#ff9102'
    
    COLOR_E_MIN     =   '#5dde00'
    COLOR_E_HALF    =   'gray95'
    COLOR_E_MAX     =   '#9102ff'
    
    # rmsd midpoint = 1 + rmsd midpoint
    RMSD_MIN        =   1 
    RMSD_half       =   2 
    RMSD_MAX        =   5
    
    # if E_MIDPOINT > 0
    # e_half = e_min + E_MIDPOINT

    # if E_MIDPOINT == 0 : e_half =  (e_min + e_max) / 2
    E_MIDPOINT = 0
    
    #
    #   END Matrix Colors
    #
    ### arguments

    Pair_Fit_selection = ''




    DIRECTORY = ''
    mandatory = ''
    RANGE = 10000    
    FF = 'MMFF94s'
    REF = ''

    try:
        opts, args = getopt.getopt(argv,"hi:e:r:n:Z:R:")
    except getopt.GetoptError:
        print ' -i input filename in mol2 [mandatory] -r rmsd cut off -n number of conformations to include in analysis [default: all] -e id-point value for the scale [default: 0]' 
        sys.exit(2)

    for opt, arg in opts:
        print "%s %s" %(opt, arg)
        if opt == '-h':
            print ' -i input filename in mol2 [mandatory] -r rmsd cut off -n number of conformations to include in analysis [default: all] -e id-point value for the scale [default: 0]' 
            sys.exit()
        elif opt in ("-i"):
            DIRECTORY = arg
            mandatory = '1'
        elif opt in ("-e"):
            E_MIDPOINT = float(arg)
        elif opt in ("-R"):
            REF = str(arg)
        elif opt in ("-n"):
            RANGE = int(arg)
    if mandatory != '1':
        print ' -i input filename in mol2 [mandatory] -r rmsd cut off -n number of conformations to include in analysis [default: all] -e mid-point value for the scale [default: 0]' 
        sys.exit()


    
    #### deal with ref

    # read first conf from DIRECTORY to check atom types (either they are conserved or they are changed from convertion mol2 to xyz resulting in a single letter atom type)
    ls = os.listdir(DIRECTORY)

    test = []
    for x in ls:
        if '.mol2' in x or '.pdb' in x or '.xyz' in x:
            test.append(x)

    f = open('%s/%s'%(DIRECTORY,test[0]),'r')
    r = f.read()
    f.close()
    check_lines = r.splitlines()
    #
    #   get atom list
    #
    
    atom_list = []
    
    start_atom = 0
    end_atom = 0
    
    for x in range(len(check_lines)):
        if start_atom != 0:
            if len(check_lines[x].split()) != 9:
                if end_atom == 0:
                    end_atom = x
    
        if '@<TRIPOS>ATOM' in check_lines[x]:
            start_atom = x+1
    
    #print 'atoms'
    #print start_atom, end_atom
    
    for x in range(start_atom, end_atom):
        atom_list.append(check_lines[x])

    tot_len = 0    

    for x in atom_list:
        sp = x.split()
        tot_len += len(sp[1])
    
    reformat = 0  

    if int(tot_len/float(len(atom_list) )) == 1:
        reformat = 1



    # reformat ref accordingly to above


    if REF != '':
        if reformat == 1:
            if '.mol2' in REF:
                os.system('%sobabel %s -O tmp_ref.xyz'%(BABEL_PATH,REF))
                os.system('%sobabel tmp_ref.xyz -O tmp_ref.mol2'%(BABEL_PATH))
                os.system('rm tmp_ref.xyz')
                os.system('mv tmp_ref.mol2 %s'%(DIRECTORY))

            else:
                sys.exit('use a .mol2 ref')
        else:
                os.system('cp %s %s/ref.mol2'%(REF,DIRECTORY))
        
    # important to redefine the ls variable after having added the ref
    ls = os.listdir(DIRECTORY)




    # 1 get conformations and energy

    mol_list1 = []

    ref_entry = []
    
    for x in ls:
        if 'mol2' in x or 'xyz' in x:
            f = open(DIRECTORY+'/'+x,'r')
            r = f.read()
            f.close()
            s = r.splitlines()
            energy = float(s[1])
            mol_list1.append([x,energy])
            if ref_entry == [] and 'tmp_ref' in x:
                ref_entry.append([x,energy])


    # get best e conf

    mol_list1.sort(key=lambda x: x[1])
    mol_list = mol_list1[:RANGE]

    #print mol_list[0]

    if REF != '':
        is_there = 0
        for x in mol_list:
            if 'tmp_ref' in x[0]:
                is_there = 1
        if is_there == 0:
            popped = mol_list.pop()
            mol_list.append(ref_entry[0])

    #print mol_list

    #print len(mol_list)
    if RANGE > len(mol_list):
        RANGE = len(mol_list)
        print 'Adjusting number of molecules to %s'%(str(len(mol_list)))
    # 2 get rmsd matrix

    Pymol_RMS = """
dirtolist = '%s'
RANGE = '%s'
REF = '%s'

ls = os.listdir(dirtolist) 
out = open('analyse_rms.dat','w')

mol_list1 = []

ref_entry = []

for x in ls:
    if 'mol2' in x or 'xyz' in x:
        f = open(dirtolist+'/'+x,'r')
        r = f.read()
        f.close()
        s = r.splitlines()
        energy = float(s[1])
        mol_list1.append([x,energy])
        if ref_entry == [] and 'tmp_ref' in x:
            ref_entry.append([x,energy])


# get best e conf

mol_list1.sort(key=lambda x: x[1])
mol_list = mol_list1[:int(RANGE)]

if REF != '':
    is_there = 0
    for x in mol_list:
        if 'tmp_ref' in x[0]:
            is_there = 1
    if is_there == 0:
        popped = mol_list.pop()
        mol_list.append(ref_entry[0])

for x in range(int(RANGE)):
    if 'mol2' in mol_list[x][0] or 'xyz' in mol_list[x][0]:
        for y in range(int(RANGE)): 
            if 'mol2' in mol_list[y][0] or 'xyz' in mol_list[y][0]:
                cmd.load(dirtolist+'/'+mol_list[x][0])
                loaded1 = mol_list[x][0].rsplit('.',1)[0]
                cmd.load(dirtolist+'/'+mol_list[y][0])
                loaded2 = mol_list[y][0].rsplit('.',1)[0]
                rmsd = cmd.pair_fit(loaded1 + ' %s',loaded2 + ' %s')
                out.write(mol_list[x][0] + '\\t' + mol_list[y][0] + '\\t' + str(rmsd)+'\\n')
                cmd.reinitialize()
out.close()

"""%(DIRECTORY,str(RANGE),str(REF),Pair_Fit_selection,Pair_Fit_selection)
    tmp_f = open('analyse_rms.py','w')
    tmp_f.write(Pymol_RMS)
    tmp_f.close()
    os.system('%spymol -cQ analyse_rms.py'%(PYMOL_PATH))

    rms_f = open('analyse_rms.dat','r')
    rms_r = rms_f.read()
    rms_f.close()
    rms_s = rms_r.splitlines()

    rms_d = {}

    for x in rms_s:
        sp = x.split()
        rms_d[(sp[0],sp[1])] = float(sp[2])

    print
    print


    mat_list = []
    e_list = []

    mat_list.append(['corner'])
    for x in range(RANGE):
        e_list.append([mol_list[x][0],str(mol_list[x][1])])
        mat_list[len(mat_list)-1].append(mol_list[x][0])
    mat_list.append([])

    for x in range(RANGE):
        mat_list[len(mat_list)-1].append(mol_list[x][0])
        for y in range(RANGE):
            mat_list[len(mat_list)-1].append( rms_d[(mol_list[x][0],mol_list[y][0])] )
        mat_list.append([])

    mat_list.pop()

    out_mat = open('rmsd_matrix_%s.dat'%(str(RANGE)),'w')
    out_energy = open('energy_list_%s.dat'%(str(RANGE)),'w')

    check_e_min = 100000000
    check_e_max = 0

    for x in range(len(mat_list)):
        for y in range(len(mat_list[x])):
            if y == 0:
                for a in mol_list:
                    if mat_list[x][y] == a[0]:
                        out_energy.write(str(a[0]) + '\t' + str(a[1]) + '\n')
                        check_e = float(a[1])
                        if check_e < check_e_min :
                            check_e_min = check_e
                        if check_e > check_e_max :
                            check_e_max = check_e
            out_mat.write(str(mat_list[x][y])+'\t')
        out_mat.write('\n')
 
    out_mat.close()
    out_energy.close()


    
    
    E_MIN           =   check_e_min 
    E_MAX           =   check_e_max
    if E_MIDPOINT == 0:
        E_half = (check_e_min + check_e_max) /2.0  
    else:
        E_half          =   check_e_min + E_MIDPOINT 

    print E_MIN, E_MAX, E_half
    
    # 3 make rmsd matrix 

    
    rcode = """


library(ComplexHeatmap)
library(circlize)

dataset <- read.table("rmsd_matrix_%s.dat", row.names = 1, header = TRUE, check.names = FALSE)
dend = hclust(dist(dataset))
conf_order <- unlist(dend["order"])
dataset2<- dataset[conf_order,conf_order]

ncol  = ncol(dataset)
ncol  = ncol/100

NCOL_MAT  = ncol*6
NCOL_FONT = 3.5
NCOL_D    = NCOL_FONT*2
NCOL_GRID = NCOL_FONT/1.8
NCOL_E  = 0.2
NCOL_W  = NCOL_MAT*1.2 + 1.5
NCOL_H  = NCOL_MAT + 0.6

Heatmap_data<- Heatmap(dataset2, column_title_gp = gpar(fontsize = NCOL_FONT, fontface = "bold"), cluster_rows = TRUE, cluster_columns = FALSE, row_dend_reorder = FALSE, column_dend_reorder = FALSE, col = colorRamp2(c(%s,%s,%s), c("%s", "%s", "%s")), column_title = "RMSD", show_column_dend = FALSE, show_column_names = FALSE, show_row_names = FALSE, name = "RMSD", row_dend_width = unit(NCOL_D, "mm"), width = unit(NCOL_MAT, "in"), heatmap_legend_param = list(labels_gp = gpar(fontsize = NCOL_FONT), title_gp = gpar(fontsize = NCOL_FONT, fontface = "bold"), color_bar = "continuous", legend_direction = "vertical", grid_height = unit(NCOL_GRID, "mm"), grid_width = unit(NCOL_GRID, "mm")))

energy <- read.table("energy_list_%s.dat", header = FALSE, row.names = 1, check.names = FALSE)
energy$V3 = energy$V2
energy_ordered <- energy[conf_order,]
energy_ordered <- energy_ordered[1]
Heatmap_energy<- Heatmap(row_names_gp = gpar(fontsize = NCOL_FONT), energy_ordered, column_title_gp = gpar(fontsize = NCOL_FONT, fontface = "bold"), col = colorRamp2(c(%s,%s,%s), c("%s","%s", "%s")), column_title = "Energy", width = unit(NCOL_E, "in"), cluster_rows = FALSE, show_column_names = FALSE, show_row_names = TRUE, name = "Energy", heatmap_legend_param = list(labels_gp = gpar(fontsize = NCOL_FONT), title_gp = gpar(fontsize = NCOL_FONT, fontface = "bold"), color_bar = "continuous", legend_direction = "vertical", grid_height = unit(NCOL_GRID, "mm"), grid_width = unit(NCOL_GRID, "mm")))

pdf("Heatmap_%s.pdf", height=NCOL_H, width=NCOL_W)
draw(Heatmap_data+Heatmap_energy)
dev.off()

"""%(str(RANGE), str(RMSD_MIN), str(RMSD_half), str(RMSD_MAX), COLOR_MAT_MIN, COLOR_MAT_HALF, COLOR_MAT_MAX, str(RANGE), str(E_MIN), str(E_half), str(E_MAX), COLOR_E_MIN, COLOR_E_HALF, COLOR_E_MAX, str(RANGE))

    out = open('make-matrix.R','w')
    out.write(rcode)
    out.close()

    # run R script
    os.system('R < make-matrix.R --no-save 1> /dev/null 2> /dev/null')


    for rm in ['analyse_rms.*','make-matrix.R']:
        try:
            os.system('rm %s'%(rm))
        except:
            null=0

if __name__ == "__main__":
    main(sys.argv[1:])
