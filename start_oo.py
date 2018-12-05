#!/usr/bin/env python
# librerias que utilizaras
import numpy as np
# por si no te lee las tools o functions creadas
import sys
# herramientas para leer pdbs
import read_pdb_tools as rpt
# calculo de distancia
from scipy.spatial.distance import pdist, squareform
# libreria de tablas
import pandas as pd
# funciones de click generadas en pandas
import funciones_CLICK as fc
# iteradores
import itertools as it
# cuenta tiempo de ejecucion
import datetime

# multiprocessing
import multiprocessing
from functools import partial
# filtro distancia
from scipy.spatial import distance


timenow = datetime.datetime.now()

# sys.path.append("../math_tricks/")
# import math_vect_tools as mymath

# assert( len(sys.argv) > 1)
# lectura de archivo
file1 = 'pdbs/1xxa.pdb'  # sys.argv[1]
file2 = 'pdbs/1tig.pdb'  # sys.argv[2]

# numero de cliques, preguntar en el software para generalizarlo...
num_cliques = 3

# outfile = open('hh_%s.txt'%infile.split('.')[0],'w')

# se define la estructura
pdb1 = rpt.PdbStruct("first")
pdb2 = rpt.PdbStruct("second")

# se lee el pdb y se agrega al objeto
pdb1.AddPdbData("%s" % file1)
pdb2.AddPdbData("%s" % file2)

# se obtienen los residuos que perteneces a la cadena de interes por default chain = 'A'
pdb11 = pdb1.GetResChain()
pdb22 = pdb2.GetResChain()

ss1 = pdb1.Get_SS(file1)
ss2 = pdb1.Get_SS(file2)

# se crea atributo a cada residuo
for i, j in zip(pdb11, ss1.structure.values):
    setattr(i, 'structure', j)
for i, j in zip(pdb22, ss2.structure.values):
    setattr(i, 'structure', j)


def get_df_distancias(ref):
    """Funcion para obtener el dataframe de distancias de cada proteina"""
    # se generan listas con coordenadas y numero de atomo
    coord = [res.GetAtom('CA').coord for res in ref]
    index = [res.resi for res in ref]

    # calcula distancia y regresa dataframe
    distancias = []
    # se calcula la distancia euclidiana entre cada atomo de carbon alfalfa
    for v in coord:
        distancia_un_atomo = []
        for av in coord:
            distancia = pdist(np.array([v, av]), metric='euclidean').item()
            distancia_un_atomo.append(distancia)
        distancias.append(distancia_un_atomo)

    # se genera la matriz de adyacencias para la red
    df_da = pd.DataFrame(index=index, columns=index, data=distancias)
    return(df_da, index)

# devuelve tabla e indices de el dataframe de distancias entre atomos de la misma proteina con dth < 10A
df_distancias1, index1 = get_df_distancias(pdb11)
df_distancias2, index2 = get_df_distancias(pdb22)

# se generan cliques, te devuleve dataframe con cliques de 3 y la lista de cliques maximales
df_cliques1, cliques1 = fc.gen_3_cliques(df_distancias1, dth=10, k=num_cliques)
print('**'*50)
df_cliques2, cliques2 = fc.gen_3_cliques(df_distancias2, dth=10, k=num_cliques)
print('**'*50)

# funcion para obtener las propiedades del residuo para los cliques agrupados
def get_df_ca(list_of_residues):
    """Genera dataframe con la informacion necesaria para las siguientes funciones
    FALTA DOCUMENTAR ESTA COSA!!!!"""
    #crear df_ca
    atom_number = []
    atom_name = []
    residue_name = []
    residue_number = []
    coord = []
    for res in list_of_residues:
        for atom in res.atoms:
            atom_number.append(atom.atom_number)
            atom_name.append(atom.name)
            residue_name.append(res.resn)
            residue_number.append(res.resi)
            coord.append(atom.coord)

    df_atoms = pd.DataFrame(columns=['atom_number', 'atom_name', 'residue_name',
                                   'residue_number', 'vector'])
    df_atoms.atom_number = atom_number
    df_atoms.atom_name = atom_name
    df_atoms.residue_name = residue_name
    df_atoms.residue_number = residue_number
    df_atoms.vector = coord

    return(df_atoms)


# CREAR DF_atomos_CA #
df_atoms1 = get_df_ca(pdb11)
df_atoms2 = get_df_ca(pdb22)

# se obtiene la estructura secundaria utilizando dssp
# ss1 = fc.mini_dssp(file1, index1)
# print('**'*50)
# ss2 = fc.mini_dssp(file2, index2)
# ya no por que se obtiene arriba desde read_pdb_tools.py

# se le pega la estructura secundaria al dataframe de los cliques
# esto va a cambiar por que lo tiene que obtener del objeto residuo
# ya se crea en ss1 y no cuesta reevaluar si es mejor desde el residuo
# checar que es mas rapido si desde residuo o desde dataframe ss
df_cliques1 = fc.paste_SS(ss1, df_cliques1, num_cliques=num_cliques)
df_cliques2 = fc.paste_SS(ss2, df_cliques2, num_cliques=num_cliques)

# comparacion SSM #aqui se obtienen los candidatos posibles pasando el filtro de SS
candidatos_ss = fc.compare_SS(df_cliques1,df_cliques2, num_cliques=num_cliques)

# obtienes coordenadas de cada atomo de los cliques.
df_cliques1 = fc.get_coords_clique(df_atoms1, df_cliques1, num_cliques)
df_cliques2 = fc.get_coords_clique(df_atoms2, df_cliques2, num_cliques)

# calculo de baricentro baricentro clique
df_cliques1 = fc.baricenter_clique(df_cliques1, num_cliques)
df_cliques2 = fc.baricenter_clique(df_cliques2, num_cliques)

# calculo de vectores gorro
df_cliques1 = fc.center_vectors(df_cliques1, num_cliques)
df_cliques2 = fc.center_vectors(df_cliques2, num_cliques)

for i, j in enumerate(df_cliques1.columns):
    print(i, j)

# para obtener el vector de columna de interes sin importar el numero de cliques.
idx_rmsd1, idx_rmsd2 = 3*num_cliques, 4*num_cliques+3
# print(list(range(idx_rmsd1,idx_rmsd2)))
# se pasan a numpy arrays para mayor rapidez
array_df_cliques1 = df_cliques1.values[:, range(idx_rmsd1, idx_rmsd2)] #del 9 al 15 columnas de interes
array_df_cliques2 = df_cliques2.values[:, range(idx_rmsd1, idx_rmsd2)]

# Se genera columna del calculo de distancia promedio para posteriormente filtrar por distancia promedio minima (dpm)
df_cliques1 = fc.get_distancia_promedio(num_cliques, df_cliques1)
df_cliques2 = fc.get_distancia_promedio(num_cliques, df_cliques2)

array_dist_promedio1 = df_cliques1.values[:, -1]  # el ultimo valor de distancia.
array_dist_promedio2 = df_cliques2.values[:, -1]

limite_distancia_minima = 0.45
if num_cliques == 4:
    limite_distancia_minima = 0.9
if num_cliques == 5:
    limite_distancia_minima = 1.8
if num_cliques == 6:
    limite_distancia_minima = 3.6
if num_cliques == 7:
    limite_distancia_minima = 4.5
if num_cliques == 8:
    limite_distancia_minima = 8.0

#filtro por distancia minima
candidatos_filter_dist = [(i, j) for i, j in candidatos_ss if (
        array_dist_promedio1[i] - array_dist_promedio2[j] >= -limite_distancia_minima) & (
        array_dist_promedio1[i] - array_dist_promedio2[j] <= limite_distancia_minima)]

print('num candidatos filtro SS', len(candidatos_ss))
print('num candidatos filtro distancia y ss', len(candidatos_filter_dist))

#calculo del RMSD

time = datetime.datetime.now()
print('tiempo pasado en filtros:', time - timenow)

timenow = datetime.datetime.now()

restriccion_rmsd = 0.15
if num_cliques == 4:
    restriccion_rmsd = 0.30
if num_cliques == 5:
    restriccion_rmsd = 0.60
if num_cliques == 6:
    restriccion_rmsd = 0.90
if num_cliques == 7:
    restriccion_rmsd = 1.50
if num_cliques == 8:
    restriccion_rmsd = 1.80

# calculo del RMSD
long_candidatos_ss = len(candidatos_filter_dist)

p = multiprocessing.Pool(multiprocessing.cpu_count()-1)

rmsd_1 = p.map(partial(fc.calculate_rmsd_rot_trans_m,
                    array_cliques1=array_df_cliques1,
                    array_cliques2=array_df_cliques2,
                    num_cliques=num_cliques), candidatos_filter_dist
                     )
p.close()
p.join()

f1 = pd.DataFrame(rmsd_1)
f1 = f1[f1[0] <= restriccion_rmsd]
candidatos = f1[1].values

time = datetime.datetime.now()

print('numero de candidatos:', len(candidatos))
print('tiempo pasado:', time - timenow)
