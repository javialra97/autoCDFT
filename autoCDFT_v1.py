#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
===============================================================================
Created on    : Fri Nov  8  2019
Author        : Javier Emilio Alfonso Ramos
Mail          : [javier.alfonso@estudiantes.fq.uh.cu, javialra97@gmail.com]
Affiliation   : Laboratory of Computational and Theoretical Chemistry,
Affiliation   : Faculty of Chemistry, University of Havana, Cuba.
===============================================================================

"""



import os, fnmatch, openpyxl


def finder(pattern, root=os.curdir):
    '''Generator expression that locates all files matching **pattern**
    argument in and inside dir **root** argument.
    Args:
        pattern (str): specified pattern.
        root (str)   : path where to begin search.
    '''
    for path, dirs, files in os.walk(os.path.abspath(root)):
        for filename in fnmatch.filter(files, pattern):
            yield os.path.join(path, filename)



"s == softness, w == electrophilicity index, N = nucleophilicity index"

def local_descriptor(charges_files, s, w, N):
    wb = openpyxl.Workbook()
    wb.get_active_sheet()
    sheet = wb.get_sheet_by_name('Sheet')
    sheet['C4'] = 'ATOM'
    sheet['D4'] = 'anion'
    sheet['E4'] = 'neutro'
    sheet['F4'] = 'cation'
    sheet['G4'] = 'f+'
    sheet['H4'] = 'f-'
    sheet['I4'] = 'DD'
    sheet['J4'] = 's+'
    sheet['K4'] = 's-'
    sheet['L4'] = 's+/s-'
    sheet['M4'] = 's-/s+'
    sheet['N4'] = 'wk'
    sheet['O4'] = 'Nk'
    sheet['J2'] = 'Softness (eV-1):'
    sheet['K2'] = s
    H = 5
    fukui_p = []
    fukui_n = []
    dd = []
    s_p = []
    s_n = []
    wk = []
    Nk = []
    descriptors = []
    for i in charges_files:
        with open(i, 'r') as file:
            if('anion') in i:
                anion = file.readlines()
            if('cation') in i:
                cation = file.readlines()
            if('neutro') in i:
                neutro = file.readlines()
    for (x, y, z) in zip(anion, neutro, cation):
        fukui_p.append(str(round(float(y.split()[4]) - float(x.split()[4]), 5)))
        fukui_n.append(str(round(float(z.split()[4]) - float(y.split()[4]), 5)))
        dd.append(str(round(2*float(y.split()[4])-float(z.split()[4]) - float(x.split()[4]),5 )))
    for (x, y) in zip(fukui_p, fukui_n):
        s_p.append(float(x)*s)
        s_n.append(float(y)*s)
        wk.append(str(float(x)*w))
        Nk.append(str(float(y)*N))
    for (p ,q ,s ,t ,u, v, w, x, y, z) in zip(s_p, s_n, wk, Nk, anion, neutro, cation, fukui_p, fukui_n, dd):
        sheet['C' + str(H)] = u.split()[0]
        sheet['D' + str(H)] = str(round(float(u.split()[4]), 5))
        sheet['E' + str(H)] = str(round(float(v.split()[4]), 5))
        sheet['F' + str(H)] = str(round(float(w.split()[4]), 5))
        sheet['G' + str(H)] = x
        sheet['H' + str(H)] = y
        sheet['I' + str(H)] = z
        sheet['J' + str(H)] = p
        sheet['K' + str(H)] = q
        sheet['L' + str(H)] = p/q
        sheet['M' + str(H)] = q/p
        sheet['N' + str(H)] = s
        sheet['O' + str(H)] = t
        descriptors.append(x + '  ' + y + '  ' + z + '  ' + str(p) + '  ' + str(q) + '  ' + str(p/q) + '  ' + str(q/p) + '  ' + s + '  ' + t)
        H = H + 1
    information = os.path.split(i)
    name = information[1]
    address = information[0]
    name = name[0: name.find('_')]
    wb.save(address + '/' + name + '_summary.xlsx')
    return descriptors


"dipolarophile == other"

"YOU SHOULD CALCULATE THIS VALUE BY YOURSELF (eV) "
"B3LYP/6-31G* "
E_HOMO_TCE = -9.121274

# ------ Extracting information about reagents ------ #
with open('information.txt', 'r') as file:
    all_information = file.readlines()
del all_information[0]
for p in all_information:
    dipole = list(finder(p.split()[0]))
    other = list(finder(p.split()[1]))
    atom_dipole = p.split()[2]
    atom_other = p.split()[3]
    with open(dipole[0], 'r') as f_dipole:
        inf_dipole = f_dipole.readlines()
    with open(other[0], 'r') as f_other:
        inf_other = f_other.readlines()
    name_dipole = os.path.split(dipole[0])[1]
    name_other = os.path.split(other[0])[1]
# -------------------------------------------------------------------------- #

# ------ Extracting energy of OF ------ #
    for j in inf_dipole:
        if j.startswith(' Alpha virt. eigenvalues'):
            break
    a = inf_dipole.index(j)
    b = a - 1
    E_HOMO_dipole = (float(inf_dipole[b].split()[len(inf_dipole[b].split())-1])*27.2116)
    E_HOMOless1_dipole = (float(inf_dipole[b].split()[len(inf_dipole[b].split())-2])*27.2116)
    E_LUMO_dipole = (float(inf_dipole[a].split()[4])*27.2116)
    E_LUMOplus1_dipole = (float(inf_dipole[a].split()[5])*27.2116)
    for k in inf_other:
        if k.startswith(' Alpha virt. eigenvalues'):
            break
    a = inf_other.index(k)
    b = a - 1
    E_HOMO_other = (float(inf_other[b].split()[len(inf_other[b].split())-1])*27.2116)
    E_HOMOless1_other = (float(inf_other[b].split()[len(inf_other[b].split())-2])*27.2116)
    E_LUMO_other = (float(inf_other[a].split()[4])*27.2116)
    E_LUMOplus1_other = (float(inf_other[a].split()[5])*27.2116)
    gap1_typeI = E_LUMO_other - E_HOMO_dipole
    gap2_typeI = E_LUMO_other - E_HOMOless1_dipole
    gap3_typeI = E_LUMOplus1_other - E_HOMO_dipole
    gap4_typeI = E_LUMOplus1_other - E_HOMOless1_dipole
    gap1_typeIII = E_LUMO_dipole - E_HOMO_other
    gap2_typeIII = E_LUMO_dipole - E_HOMOless1_other
    gap3_typeIII = E_LUMOplus1_dipole - E_HOMO_other
    gap4_typeIII = E_LUMOplus1_dipole - E_HOMOless1_other
    gaps = [gap1_typeI, gap2_typeI, gap3_typeI, gap4_typeI, gap1_typeIII,
            gap2_typeIII, gap3_typeIII, gap4_typeIII]
    "m == minus      looking for the lowest gap"
    m = gap1_typeI
    for i in gaps:
        if m > i:
            m = i
# -------------------------------------------------------------------------- #

# ------ Extracting Molecular Orbital Coefficients ------ #
    atom_dipole1 = []
    for i in range(int(atom_dipole[0])):
        aaa = atom_dipole.find('-', 2)
        c = atom_dipole[2:aaa]
        atom_dipole = atom_dipole.replace(atom_dipole[2:aaa + 1], '')
        atom_dipole1.append(c)
    atom_other1 = []
    for i in range(int(atom_other[0])):
        aaa = atom_other.find('-', 2)
        c = atom_other[2:aaa]
        atom_other = atom_other.replace(atom_other[2:aaa + 1], '')
        atom_other1.append(c)
# -------------------------------------------------------------------------- #

# ------ Calculating c-DFT ------ #
    chemical_potential_dipole = (E_HOMO_dipole + E_LUMO_dipole)/2
    chemical_potential_other = (E_HOMO_other + E_LUMO_other)/2
    hardness_dipole = (E_LUMO_dipole - E_HOMO_dipole)
    hardness_other = (E_LUMO_other - E_HOMO_other)
    softness_dipole = 1/hardness_dipole
    softness_other = 1/hardness_other
    electrophilicity_index_dipole = (chemical_potential_dipole*chemical_potential_dipole)/(2 * hardness_dipole)
    electrophilicity_index_other = (chemical_potential_other*chemical_potential_other)/(2 * hardness_other)
    nucleophilicity_index_dipole = E_HOMO_dipole - E_HOMO_TCE
    nucleophilicity_index_other = E_HOMO_other - E_HOMO_TCE
    charges_files_d = list(finder(name_dipole[0: name_dipole.find('.')] + '*.chg'))
    charges_files_o = list(finder(name_other[0: name_other.find('.')] + '*.chg'))
    cdft_dipole = local_descriptor(charges_files_d, softness_dipole, electrophilicity_index_dipole, nucleophilicity_index_dipole)
    cdft_other = local_descriptor(charges_files_o, softness_other, electrophilicity_index_other, nucleophilicity_index_other)
# -------------------------------------------------------------------------- #

# ------ Writing the output ------ #
    name = name_dipole[0: name_dipole.find('.')] + '_' + name_other[0: name_other.find('.')] + '.txt'
    with open(name, 'w') as at_last:
        at_last.write('\n\n 1,3 Dipole:     ' + name_dipole[0: name_dipole.find('.')])
        at_last.write('\n Dipolarophile:  ' + name_other[0: name_other.find('.')])
        at_last.write('\n\n\n')
        at_last.write(' **********************************************************************\n')
        at_last.write('\n                        Orbital Energies\n\n')
        at_last.write(' **********************************************************************\n\n')
        at_last.write('                  1,3 dipole         ||          dipolarophile\n')
        if (E_LUMOplus1_dipole < 0):
            d = 0
        else:
            d = 1
        if (E_LUMOplus1_other < 0):
            e = 0
        else:
            e = 1
        at_last.write(' LUMO + 1          ' + d*' ' + str(round(E_LUMOplus1_dipole, 2)) + '             ||             ' + e*' ' + str(round(E_LUMOplus1_other, 2)) + '\n')
        if (E_LUMO_dipole < 0):
            d = 0
        else:
            d = 1
        if (E_LUMO_other < 0):
            e = 0
        else:
            e = 1
        at_last.write(' LUMO              ' + d*' ' + str(round(E_LUMO_dipole, 2)) + '             ||             ' + e*' ' + str(round(E_LUMO_other, 2)) + '\n')
        at_last.write(' HOMO              ' + str(round(E_HOMO_dipole, 2)) + '             ||             ' + str(round(E_HOMO_other, 2)) + '\n')
        at_last.write(' HOMO - 1          ' + str(round(E_HOMOless1_dipole, 2)) + '             ||             ' + str(round(E_HOMO_other, 2)) + '\n')
        at_last.write('\n\n **********************************************************************\n')
        at_last.write('\n                               Gaps\n\n')
        at_last.write(' **********************************************************************\n\n')
        at_last.write(' LUMO    dipolarophile -  HOMO    dipole:         ')
        at_last.write(str(round(gap1_typeI, 2)) + '\n')
        at_last.write(' LUMO    dipolarophile - (HOMO-1) dipole:         ')
        at_last.write(str(round(gap2_typeI, 2)) + '\n')
        at_last.write('(LUMO+1) dipolarophile -  HOMO    dipole:         ')
        at_last.write(str(round(gap3_typeI, 2)) + '\n')
        at_last.write('(LUMO+1) dipolarophile - (HOMO-1) dipole:         ')
        at_last.write(str(round(gap4_typeI, 2)) + '\n')
        at_last.write(' LUMO    dipole        -  HOMO    dipolarophile:  ')
        at_last.write(str(round(gap1_typeIII, 2)) + '\n')
        at_last.write(' LUMO    dipole        - (HOMO-1) dipolarophile:  ')
        at_last.write(str(round(gap2_typeIII, 2)) + '\n')
        at_last.write('(LUMO+1) dipole        -  HOMO    dipolarophile:  ')
        at_last.write(str(round(gap3_typeIII, 2)) + '\n')
        at_last.write('(LUMO+1) dipole        - (HOMO-1) dipolarophile:  ')
        at_last.write(str(round(gap4_typeIII, 2)) + '\n')
        at_last.write('\n\n **********************************************************************\n')
        at_last.write('\n                               c-DFT\n\n')
        at_last.write(' **********************************************************************\n\n')
        at_last.write('\n    Global Indices\n')
        at_last.write(' Chemical Potential dipole:         ' + str(round(chemical_potential_dipole, 3)) + '\n')
        at_last.write(' Chemical Potential dipolarophile:  ' + str(round(chemical_potential_other, 3)) + '\n')
        at_last.write(' Hardness dipole:         ' + str(round(hardness_dipole, 3)) + '\n')
        at_last.write(' Hardness dipolarophile:  ' + str(round(hardness_other, 3)) + '\n')
        at_last.write(' Softness dipole:         ' + str(round(softness_dipole, 3)) + '\n')
        at_last.write(' Softness dipolarophile:  ' + str(round(softness_other, 3)) + '\n')
        at_last.write(' Electrophilicity Index dipole:         ' + str(round(electrophilicity_index_dipole, 3)) + '\n')
        at_last.write(' Electrophilicity Index dipolarophile:  ' + str(round(electrophilicity_index_other, 3)) + '\n')
        at_last.write('\n    Atom Indices Dipole\n')
        at_last.write('\n    Atom         f+         f-           dd          s+          s-        s+/s-       s-/s+         w            N')
        for i in atom_dipole1:
            jj = 5 + len (i) - 1
            if (float(cdft_dipole[int(i)-1].split()[0]) > 0):
                aa = 1
            else: aa = 0
            if (float(cdft_dipole[int(i)-1].split()[1]) > 0):
                bb = 1
            else: bb = 0
            if (float(cdft_dipole[int(i)-1].split()[2]) > 0):
                cc = 1
            else: cc = 0
            if (float(cdft_dipole[int(i)-1].split()[3]) > 0):
                dd = 1
            else: dd = 0
            if (float(cdft_dipole[int(i)-1].split()[4]) > 0):
                ee = 1
            else: ee = 0
            if (float(cdft_dipole[int(i)-1].split()[5]) > 0):
                ff = 1
            else: ff = 0
            if (float(cdft_dipole[int(i)-1].split()[6]) > 0):
                gg = 1
            else: gg = 0
            if (float(cdft_dipole[int(i)-1].split()[7]) > 0):
                hh = 1
            else: hh = 0
            if (float(cdft_dipole[int(i)-1].split()[8]) > 0):
                ii = 1
            else: ii = 0
            at_last.write('\n' + jj*' ' + i + '       ' + aa*' ' + str(round(float(cdft_dipole[int(i)-1].split()[0]), 5)))
            at_last.write(bb*' ' + '    ' +str(round(float(cdft_dipole[int(i)-1].split()[1]), 5)))
            at_last.write(cc*' ' + '    ' +str(round(float(cdft_dipole[int(i)-1].split()[2]), 5)))
            at_last.write(dd*' ' + '    ' +str(round(float(cdft_dipole[int(i)-1].split()[3]), 5)))
            at_last.write(ee*' ' + '    ' +str(round(float(cdft_dipole[int(i)-1].split()[4]), 5)))
            at_last.write(ff*' ' + '    ' +str(round(float(cdft_dipole[int(i)-1].split()[5]), 5)))
            at_last.write(gg*' ' + '    ' +str(round(float(cdft_dipole[int(i)-1].split()[6]), 5)))
            at_last.write(hh*' ' + '    ' +str(round(float(cdft_dipole[int(i)-1].split()[7]), 5)))
            at_last.write(ii*' ' + '    ' +str(round(float(cdft_dipole[int(i)-1].split()[8]), 5)))
        at_last.write('\n\n\n    Atom Indices Other\n')
        at_last.write('\n    Atom         f+         f-           dd          s+          s-        s+/s-       s-/s+         w            N')
        for i in atom_other1:
            jj = 5 + len (i) - 1
            if (float(cdft_other[int(i)-1].split()[0]) > 0):
                aa = 1
            else: aa = 0
            if (float(cdft_other[int(i)-1].split()[1]) > 0):
                bb = 1
            else: bb = 0
            if (float(cdft_other[int(i)-1].split()[2]) > 0):
                cc = 1
            else: cc = 0
            if (float(cdft_other[int(i)-1].split()[3]) > 0):
                dd = 1
            else: dd = 0
            if (float(cdft_other[int(i)-1].split()[4]) > 0):
                ee = 1
            else: ee = 0
            if (float(cdft_other[int(i)-1].split()[5]) > 0):
                ff = 1
            else: ff = 0
            if (float(cdft_other[int(i)-1].split()[6]) > 0):
                gg = 1
            else: gg = 0
            if (float(cdft_other[int(i)-1].split()[7]) > 0):
                hh = 1
            else: hh = 0
            if (float(cdft_other[int(i)-1].split()[8]) > 0):
                ii = 1
            else: ii = 0
            at_last.write('\n' + jj*' ' + i + '       ' + aa*' ' + str(round(float(cdft_other[int(i)-1].split()[0]), 5)))
            at_last.write(bb*' ' + '    ' +str(round(float(cdft_other[int(i)-1].split()[1]), 5)))
            at_last.write(cc*' ' + '    ' +str(round(float(cdft_other[int(i)-1].split()[2]), 5)))
            at_last.write(dd*' ' + '    ' +str(round(float(cdft_other[int(i)-1].split()[3]), 5)))
            at_last.write(ee*' ' + '    ' +str(round(float(cdft_other[int(i)-1].split()[4]), 5)))
            at_last.write(ff*' ' + '    ' +str(round(float(cdft_other[int(i)-1].split()[5]), 5)))
            at_last.write(gg*' ' + '    ' +str(round(float(cdft_other[int(i)-1].split()[6]), 5)))
            at_last.write(hh*' ' + '    ' +str(round(float(cdft_other[int(i)-1].split()[7]), 5)))
            at_last.write(ii*' ' + '    ' +str(round(float(cdft_other[int(i)-1].split()[8]), 5)))
# -------------------------------------------------------------------------- #