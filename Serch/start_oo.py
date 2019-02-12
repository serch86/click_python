#!/usr/bin/env python
# librerias que utilizaras
import numpy as np
# por si no te lee las tools o functions creadas
import sys
# herramientas para leer pdbs
import read_pdb_tools as rpt
# calculo de distancia
from scipy.spatial.distance import pdist
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

timenow_bueno = datetime.datetime.now()

timenow = datetime.datetime.now()

# sys.path.append("../math_tricks/")
# import math_vect_tools as mymath

# lectura de archivo
file1 = 'pdbs/1xxa_clean.pdb'  # sys.argv[1]
file2 = 'pdbs/1tig_clean.pdb'  # sys.argv[2]

# numero de cliques, preguntar en el software para generalizarlo INPUT...
number_elements_clique = 3

# se define la estructura
pdb1 = rpt.PdbStruct(file1)
pdb2 = rpt.PdbStruct(file2)

# se lee el pdb y se agrega al objeto
pdb1.AddPdbData("%s" % file1)
pdb2.AddPdbData("%s" % file2)

# Se calculan sus vecinos mas cercanos
pdb1.GetNeighbors()
pdb2.GetNeighbors()

# se obtienen los residuos que perteneces a la cadena de interes por default chain = 'A'
pdb11 = pdb1.GetResChain()
pdb22 = pdb2.GetResChain()

# Siempre la proteina 1 es el que se rota y traslada para embonar en la proteina 2
if len(pdb22) < len(pdb11):

    import copy
    pdb1_temp = copy.copy(pdb1)
    pdb2_temp = copy.copy(pdb2)

    pdb11_temp = copy.copy(pdb11)
    pdb22_temp = copy.copy(pdb22)

    pdb1 = pdb2_temp
    pdb2 = pdb1_temp

    pdb11 = pdb22_temp
    pdb22 = pdb11_temp

    del [pdb1_temp]
    del [pdb2_temp]
    del [pdb11_temp]
    del [pdb22_temp]

    print("Intercambio de nombre ya que la proteina 1 es mas grande que la 2")
    print(pdb1.name, len(pdb11))
    print(pdb2.name, len(pdb22))
    print("No te preocupes ya quedo :)")


pdb1.Set_SS()
pdb2.Set_SS()

ss1 = fc.create_ss_table(pdb11)
ss2 = fc.create_ss_table(pdb22)


def get_df_distancias(ref):
    """Funcion para obtener el dataframe de distancias de cada residuo
    Dudas en codigo pueden revisar fc.distancia_entre_atomos en ese se basa
    esta funcion, la diferencia es que se crea con el objeto residuo"""
    # se generan listas con coordenadas y numero de atomo
    coord = [res.GetAtom('CA').coord for res in ref]
    lista_residuos = [res.resi for res in ref]

    # calcula distancia y regresa dataframe
    distancias = [ [pdist(np.array([v, av]), metric='euclidean').item() for av in coord] for v in coord]

    # se genera la matriz de adyacencias para la red
    df_da = pd.DataFrame(index=lista_residuos, columns=lista_residuos, data=distancias)
    return(df_da, lista_residuos)


# devuelve tabla e indices de el dataframe de distancias entre atomos de la misma proteina con dth < 10A
df_distancias1, lista_residuos_1 = get_df_distancias(pdb11)
df_distancias2, lista_residuos2 = get_df_distancias(pdb22)

# se generan cliques, te devuleve dataframe con cliques de k(numero_de_cliques) y la lista de cliques maximales
df_cliques1, cliques1 = fc.gen_3_cliques(df_distancias1,  dth=10, k=number_elements_clique)  # file1[:-4] pal gexf
print('**'*50)
df_cliques2, cliques2 = fc.gen_3_cliques(df_distancias2,  dth=10, k=number_elements_clique)  # file2[:-4]
print('**'*50)

# REVISAR DONDE GUARDA EL GEXF Y COLOCARLO EN LA CARPETA Grafos

# se agrega filtro de residuos que pertenecen a cliques de 7 elementos
mc1_7 = cliques1[cliques1.numero_elementos == 7].drop('numero_elementos', 1)
mc2_7 = cliques2[cliques2.numero_elementos == 7].drop('numero_elementos', 1)

residuos_unicos_1 = []
for i in [list(mc1_7[i].unique()) for i in mc1_7.dropna(1).columns]:
    for j in i:
        residuos_unicos_1.append(int(j))

residuos_unicos_2 = []
for i in [list(mc2_7[i].unique()) for i in mc2_7.dropna(1).columns]:
    for j in i:
        residuos_unicos_2.append(int(j))

residuos_unicos_1 = set(residuos_unicos_1)
residuos_unicos_2 = set(residuos_unicos_2)


for i in df_cliques1.columns:
    mask = np.where(df_cliques1[i].isin(residuos_unicos_1), True, False)
    df_cliques1 = df_cliques1[mask].reset_index(drop=True)

for i in df_cliques2.columns:
    mask = np.where(df_cliques2[i].isin(residuos_unicos_2), True, False)
    df_cliques2 = df_cliques2[mask].reset_index(drop=True)
# acaba filtro de cliques con 7 elementos

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
    res_ngb = []
    for res in list_of_residues:
        for atom in res.atoms:
            atom_number.append(atom.atom_number)
            atom_name.append(atom.name)
            residue_name.append(res.resn)
            residue_number.append(res.resi)
            coord.append(atom.coord)
            res_ngb.append(res.ngb)

    df_atoms = pd.DataFrame(columns=['atom_number', 'atom_name', 'residue_name',
                                   'residue_number', 'vector', 'ngb'])
    df_atoms.atom_number = atom_number
    df_atoms.atom_name = atom_name
    df_atoms.residue_name = residue_name
    df_atoms.residue_number = residue_number
    df_atoms.vector = coord
    df_atoms.ngb = res_ngb

    return(df_atoms)


# CREAR DF_atomos_CA #
df_atoms1 = get_df_ca(pdb11)
df_atoms2 = get_df_ca(pdb22)

# se genera la columna de ss a la tabla de los cliques
df_cliques1 = fc.paste_SS(ss1, df_cliques1, num_cliques=number_elements_clique)
df_cliques2 = fc.paste_SS(ss2, df_cliques2, num_cliques=number_elements_clique)

# comparacion SSM #aqui se obtienen los candidatos posibles pasando el filtro de SS
candidatos_ss = fc.compare_SS(df_cliques1, df_cliques2, num_cliques=number_elements_clique)

# obtienes coordenadas de cada atomo de los cliques.
df_cliques1 = fc.get_coords_clique(df_atoms1, df_cliques1, number_elements_clique)
df_cliques2 = fc.get_coords_clique(df_atoms2, df_cliques2, number_elements_clique)

# calculo de baricentro baricentro clique
df_cliques1 = fc.baricenter_clique(df_cliques1, number_elements_clique)
df_cliques2 = fc.baricenter_clique(df_cliques2, number_elements_clique)

# calculo de vectores gorro
df_cliques1 = fc.center_vectors(df_cliques1, number_elements_clique)
df_cliques2 = fc.center_vectors(df_cliques2, number_elements_clique)

# para obtener el vector de columna de interes sin importar el numero de cliques.
idx_rmsd1, idx_rmsd2 = 3 * number_elements_clique, 4 * number_elements_clique + 3
# print(list(range(idx_rmsd1,idx_rmsd2)))
# se pasan a numpy arrays para mayor rapidez
array_df_cliques1 = df_cliques1.values[:, range(idx_rmsd1, idx_rmsd2)]  # del 9 al 15 columnas de interes
array_df_cliques2 = df_cliques2.values[:, range(idx_rmsd1, idx_rmsd2)]

# Se genera columna del calculo de distancia promedio para posteriormente filtrar por distancia promedio minima (dpm)
df_cliques1 = fc.get_distancia_promedio(number_elements_clique, df_cliques1)
df_cliques2 = fc.get_distancia_promedio(number_elements_clique, df_cliques2)

array_dist_promedio1 = df_cliques1.values[:, -1]  # el ultimo valor de distancia.
array_dist_promedio2 = df_cliques2.values[:, -1]

limite_distancia_minima = 0.45
if number_elements_clique == 4:
    limite_distancia_minima = 0.9
if number_elements_clique == 5:
    limite_distancia_minima = 1.8
if number_elements_clique == 6:
    limite_distancia_minima = 3.6
if number_elements_clique == 7:
    limite_distancia_minima = 4.5
if number_elements_clique == 8:
    limite_distancia_minima = 8.0

# filtro por distancia minima
candidatos_filter_dist = [(i, j) for i, j in candidatos_ss if (
        array_dist_promedio1[i] - array_dist_promedio2[j] >= -limite_distancia_minima) & (
        array_dist_promedio1[i] - array_dist_promedio2[j] <= limite_distancia_minima)]

print('num candidatos filtro SS', len(candidatos_ss))
print('num candidatos filtro distancia y ss', len(candidatos_filter_dist))
#Termina filtro distancia minima

# calculo del RMSD

time = datetime.datetime.now()
print('tiempo pasado en filtros:', time - timenow)

timenow = datetime.datetime.now()

restriccion_rmsd = 0.15
if number_elements_clique == 4:
    restriccion_rmsd = 0.30
if number_elements_clique == 5:
    restriccion_rmsd = 0.60
if number_elements_clique == 6:
    restriccion_rmsd = 0.90
if number_elements_clique == 7:
    restriccion_rmsd = 1.50
if number_elements_clique == 8:
    restriccion_rmsd = 1.80

# calculo del RMSD
long_candidatos_ss = len(candidatos_filter_dist)

p = multiprocessing.Pool(multiprocessing.cpu_count()-1)

rmsd_1 = p.map(partial(fc.calculate_rmsd_rot_trans_m,
                       array_cliques1=array_df_cliques1,
                       array_cliques2=array_df_cliques2,
                       num_cliques=number_elements_clique), candidatos_filter_dist
               )

p.close()
p.join()

f1 = pd.DataFrame(rmsd_1)
f1 = f1[f1[0] <= restriccion_rmsd]
f1['candidato_clique_1'] = f1[1].str.get(0)
f1['candidato_clique_2'] = f1[1].str.get(1)
candidatos = f1[1].values

time = datetime.datetime.now()

print('numero de candidatos:', len(candidatos))
print('tiempo pasado:', time - timenow)

new_df_cliques1 = fc.add_element_clique(df_cliques1, 'candidato_clique_1', cliques1, f1, number_elements_clique)
new_df_cliques2 = fc.add_element_clique(df_cliques2, 'candidato_clique_2', cliques2, f1, number_elements_clique)
number_elements_clique = number_elements_clique + 1


def iter_rmsd(new_df_cliques1,new_df_cliques2,number_elements_clique):
    # filtro residuos de 7 unicos
    for i in new_df_cliques1.columns:
        mask = np.where(new_df_cliques1[i].isin(residuos_unicos_1), True, False)
        new_df_cliques1 = new_df_cliques1[mask].reset_index(drop=True)

    for i in new_df_cliques2.columns:
        mask = np.where(new_df_cliques2[i].isin(residuos_unicos_2), True, False)
        new_df_cliques2 = new_df_cliques2[mask].reset_index(drop=True)

    # filtro estructura secundaria
    new_df_cliques1 = fc.paste_SS(ss1, new_df_cliques1, num_cliques=number_elements_clique)
    new_df_cliques2 = fc.paste_SS(ss2, new_df_cliques2, num_cliques=number_elements_clique)
    candidatos_ss = fc.compare_SS(new_df_cliques1, new_df_cliques2, num_cliques=number_elements_clique)

    # rotando y trasladando
    df_cliques1 = fc.get_coords_clique(df_atoms1, new_df_cliques1, number_elements_clique)
    df_cliques2 = fc.get_coords_clique(df_atoms2, new_df_cliques2, number_elements_clique)
    df_cliques1 = fc.baricenter_clique(df_cliques1, number_elements_clique)
    df_cliques2 = fc.baricenter_clique(df_cliques2, number_elements_clique)
    df_cliques1 = fc.center_vectors(df_cliques1, number_elements_clique)
    df_cliques2 = fc.center_vectors(df_cliques2, number_elements_clique)
    idx_rmsd1, idx_rmsd2 = 3 * number_elements_clique, 4 * number_elements_clique + 3
    array_df_cliques1 = df_cliques1.values[:, range(idx_rmsd1, idx_rmsd2)]  # del 9 al 15
    array_df_cliques2 = df_cliques2.values[:, range(idx_rmsd1, idx_rmsd2)]

    # filtro de distancia minima
    df_cliques1 = fc.get_distancia_promedio(number_elements_clique, df_cliques1)
    df_cliques2 = fc.get_distancia_promedio(number_elements_clique, df_cliques2)
    array_dist_promedio1 = df_cliques1.values[:, -1]  # el ultimo valor de distancia.
    array_dist_promedio2 = df_cliques2.values[:, -1]
    limite_distancia_minima = 0.45
    if number_elements_clique == 4:
        limite_distancia_minima = 0.9
    if number_elements_clique == 5:
        limite_distancia_minima = 1.8
    if number_elements_clique == 6:
        limite_distancia_minima = 3.6
    if number_elements_clique == 7:
        limite_distancia_minima = 4.5
    if number_elements_clique == 8:
        limite_distancia_minima = 8.0

    # filtro por distancia minima
    candidatos_filter_dist = [(i, j) for i, j in candidatos_ss if (
            array_dist_promedio1[i] - array_dist_promedio2[j] >= -limite_distancia_minima) & (
                                      array_dist_promedio1[i] - array_dist_promedio2[j] <= limite_distancia_minima)]

    #filtro por restriccion de RMSD despues de ajuste 3D
    p = multiprocessing.Pool(multiprocessing.cpu_count() - 1)
    restriccion_rmsd = 0.15
    if number_elements_clique == 4:
        restriccion_rmsd = 0.30
    if number_elements_clique == 5:
        restriccion_rmsd = 0.60
    if number_elements_clique == 6:
        restriccion_rmsd = 0.90
    if number_elements_clique == 7:
        restriccion_rmsd = 1.50
    if number_elements_clique == 8:
        restriccion_rmsd = 1.80

    rmsd_1 = p.map(partial(fc.calculate_rmsd_rot_trans_m,
                           array_cliques1=array_df_cliques1,
                           array_cliques2=array_df_cliques2,
                           num_cliques=number_elements_clique), candidatos_filter_dist)

    p.close()
    p.join()

    df_candidatos = pd.DataFrame(rmsd_1, columns=['rmsd', 'candidatos', 'matriz_rotacion'])
    df_candidatos = df_candidatos[df_candidatos.rmsd <= restriccion_rmsd]
    df_candidatos['candidato_clique_1'] = df_candidatos.candidatos.str.get(0)
    df_candidatos['candidato_clique_2'] = df_candidatos.candidatos.str.get(1)
    candidatos = df_candidatos.candidatos.values
    time = datetime.datetime.now()

    print('numero de candidatos:', len(candidatos))
    print('tiempo pasado:', time - timenow)

    # Se agrega un nuevo elemento a los cliques.
    if number_elements_clique < 7:
        new_df_cliques1 = fc.add_element_clique(df_cliques1, 'candidato_clique_1', cliques1, df_candidatos, number_elements_clique)
        new_df_cliques2 = fc.add_element_clique(df_cliques2, 'candidato_clique_2', cliques2, df_candidatos, number_elements_clique)
        number_elements_clique = number_elements_clique + 1

    return(new_df_cliques1, new_df_cliques2, number_elements_clique, candidatos,df_candidatos)


# Si se empieza en 4-cliques solo itera 3 veces.
for k in range(4):
    new_df_cliques1, new_df_cliques2, number_elements_clique, candidatos, rmsd = iter_rmsd(
        new_df_cliques1, new_df_cliques2, number_elements_clique)
    print("iteracion", k + 1)
    print("===" * 10)


# Obtenemos los elementos para generar el Alineamiento final.
# Utilizando los datos de los cliques, el numero de elementos del clique, la lista de candidatos que sobrevivieron
# y el dataframe del indice de candidatos con el RMSD de cada alineamiento
# quitando el clique seleccionado
rmsd.reset_index(drop=True, inplace=True)


# se obtiene la matriz de rotacion del menor rmsd
# se aplica a todos los vectores gorro de la proteina 1 que ya se le quito el baricentro del candidato 1
# para cada candidato
candidato = []
protein_to_compare = np.array([vectors for vectors in df_atoms2[df_atoms2.atom_name == 'CA'].vector.values])

df_rmsd = rmsd

# CHECK DE ALINUEAMINETO
# df_rmsd.to_pickle('df_rmsd.pkl')
#
# new_df_cliques1.to_pickle('clique1.pkl')
# new_df_cliques2.to_pickle('clique2.pkl')
#
# df_atoms1.to_pickle('clique1_df_atributos.pkl')
# df_atoms2.to_pickle('clique2_df_atributos.pkl')
# vectores_gorro_proteina_1.to_pickle('vectore_gorro_1.pkl')
#
# # df_candidatos.to_pickle('df_alineados.pkl')
# # pd.DataFrame(parejas_clique).to_pickle('parejas.pkl')
# exit()

def align_residues(idx, so):
    df_new_df_clique1 = new_df_cliques1.values
    df_new_df_clique2 = new_df_cliques2.values

    df_rmsd_1 = df_rmsd.values

    clique_in_protein_1 = df_rmsd_1[idx, 3]  # se extrae el indice de los cliques
    clique_in_protein_2 = df_rmsd_1[idx, 4]

    # se obtienen el numero de residuo de los candidatos
    los_del_clique_1 = df_new_df_clique1[clique_in_protein_1, [0, 1, 2, 3, 4, 5, 6]]

    los_del_clique_2 = df_new_df_clique2[clique_in_protein_2, [0, 1, 2, 3, 4, 5, 6]]

    parejas_cliques = [(i, j) for i, j in
                       zip(los_del_clique_1, los_del_clique_2)]  # se crean las parejas de residuo de los cliques

    lista_vectores_gorro = []
    for bari in new_df_cliques1.baricentro_clique.values:
        lista_pre_vectores = [coord - bari for coord in df_atoms1[
            df_atoms1.atom_name == 'CA'][
            df_atoms1.residue_number.isin(los_del_clique_1) == False].vector.values]

        lista_vectores_gorro.append(lista_pre_vectores)

    vectores_gorro_proteina_1 = pd.DataFrame(lista_vectores_gorro)

    df_vectores_gorro_proteina_1 = vectores_gorro_proteina_1.values
    # OJO se extrae por indice no por numero de atomos!!!
    atoms1 = df_atoms1[df_atoms1.atom_name == 'CA'][
        df_atoms1.residue_number.isin(los_del_clique_1) == False].index.values
    atoms2 = df_atoms2[df_atoms2.atom_name == 'CA'][
        df_atoms2.residue_number.isin(los_del_clique_2) == False].index.values
    # valor para SO
    number_of_residues = df_atoms1.iloc[atoms1].shape[0]
    # se guardan los vectores
    protein_to_compare = np.array([vectors for vectors in df_atoms2.iloc[atoms2].vector.values])

    matriz_rotacion = df_rmsd_1[idx, 2]  # se obtiene la matriz de rotacion de los cliques candidatos

    # se obtienen los vectores gorro de la proteina a rotar y trasladar
    vector_gorro = df_vectores_gorro_proteina_1[clique_in_protein_1]  # falta quitar los que no son del clique
    # se aplica la rotacion a vectores gorro
    coord_vectores_rotados = [np.matmul(matriz_rotacion, i.reshape(3, 1)).T[0] for i in vector_gorro]
    # se obtiene el baricentro de la proteina 2
    baricentro_proteina_2 = df_new_df_clique2[clique_in_protein_2, 22]
    # se suma el baricentro a vectores rotados de la proteina 1
    protein_trasladado_rotado = coord_vectores_rotados + baricentro_proteina_2  # nuevas coordendas proteina 1
    # se quitan residuos que pertenezcan al clique candidato para calcular la distancia y emparejar al mejor.

    protein_trasladado_rotado_sin_clique = protein_trasladado_rotado
    protein_to_compare_sin_clique = protein_to_compare

    atom_numbers1 = atoms1
    atom_numbers2 = atoms2

    residuos_match = []  # aqui se guardan las parejas de residuos

    # se itera por cada residuo ya rotado y trasladado

    for residuo_not_in_clique, res_1 in zip(protein_trasladado_rotado_sin_clique, atom_numbers1):
        for residuo, res_2 in zip(protein_to_compare_sin_clique, atom_numbers2):
            distancia = pdist(np.array([residuo_not_in_clique, residuo]), metric='euclidean').item()
            if distancia < 3.5:
                residuos_match.append([distancia, (res_1, res_2)])

    df = pd.DataFrame(residuos_match, columns=['distancia', 'candigatos'])
    df = df.sort_values('distancia')
    df['cand1'] = df.candigatos.str[0]
    df['cand2'] = df.candigatos.str[1]

    c1 = []
    c2 = []
    cand_n = []
    for i in df.values:
        if (i[2] in c1) or (i[3] in c2) or (i[0] > 3.5):
            continue
        else:
            c1.append(i[2])
            c2.append(i[3])
            cand_n.append(i)

    df = pd.DataFrame(cand_n, columns=['distancia', 'candigatos', 'cand_1', 'cand2'])

    vecs1 = np.array([vecs for vecs in df_atoms1.vector.values])
    vecs2 = np.array([vecs for vecs in df_atoms2.vector.values])

    # baricentro
    bari_1 = vecs1.mean(0)
    bari_2 = vecs2.mean(0)

    # vectores centricos
    vecs_center_1 = vecs1 - bari_1
    vecs_center_2 = vecs2 - bari_2

    parejas = df.candigatos.values
    parejas = list(parejas)

    def R_ij(i, j):
        valor = sum([vecs_center_1[pos_res[0] - 1, i] * vecs_center_2[pos_res[1] - 1, j] for pos_res in parejas])
        return valor

    # temp = pd.DataFrame(parejas)[1].values

    def matrix_R():

        """cliques a comparar: i,j
        desde aqui se itera sobre cada i y hay que variar los vectores
        coordenada
        Regresa la matriz gigante (matriz simetrica del articulo)"""
        # primer renglon
        R11R22R33 = (R_ij(0, 0) + R_ij(1, 1) + R_ij(2, 2))
        R23_R32 = (R_ij(1, 2) - R_ij(2, 1))
        R31_R13 = (R_ij(2, 0) - R_ij(0, 2))
        R12_R21 = (R_ij(0, 1) - R_ij(1, 0))
        # segundo renglon
        R11_R22_R33 = (R_ij(0, 0) - R_ij(1, 1) - R_ij(2, 2))
        R12R21 = (R_ij(0, 1) + R_ij(1, 0))
        R13R31 = (R_ij(0, 2) + R_ij(2, 0))
        # tercer renglon
        _R11R22_R33 = (-R_ij(0, 0) + R_ij(1, 1) - R_ij(2, 2))
        R23R32 = (R_ij(1, 2) + R_ij(2, 1))
        # cuarto renglon
        _R11_R22R33 = (-R_ij(0, 0) - R_ij(1, 1) + R_ij(2, 2))

        matriz_R = [
            [R11R22R33, R23_R32, R31_R13, R12_R21],
            [R23_R32, R11_R22_R33, R12R21, R13R31],
            [R31_R13, R12R21, _R11R22_R33, R23R32],
            [R12_R21, R13R31, R23R32, _R11_R22R33]
        ]
        return (matriz_R)

    def rotation_matrix(matriz_R):
        """utilizando la funcion giant_matrix, fijando los valores de i,j
        se calcula la matriz de rotacion con los eigenvectores y eigenvalores
        arroja una matriz de rotacion que depende de la matriz gigante
        """
        eignvalues, eigenvectors = np.linalg.eig(matriz_R)
        q = eigenvectors[:, np.argmax(eignvalues)]
        q0, q1, q2, q3 = q[0], q[1], q[2], q[3]
        # matriz de rotacion con eigenvectores
        matriz_rotacion = np.array([
            [(q0 ** 2 + q1 ** 2 - q2 ** 2 - q3 ** 2), 2 * (q1 * q2 - q0 * q3), 2 * (q1 * q3 + q0 * q2)],
            [2 * (q1 * q2 + q0 * q3), (q0 ** 2 - q1 ** 2 + q2 ** 2 - q3 ** 2), 2 * (q2 * q3 - q0 * q1)],
            [2 * (q1 * q3 - q0 * q2), 2 * (q2 * q3 + q0 * q1), (q0 ** 2 - q1 ** 2 - q2 ** 2 + q3 ** 2)]
        ], dtype=np.float64)
        return (matriz_rotacion)

    def rotation_vectors(vector_gorro, matriz_rotacion):
        """obtencion de vector rotado,
        utilizando la matriz de rotacion
        y los vectores gorro a rotar y trasladar"""

        coord_rotado = [np.matmul(
            matriz_rotacion, i.reshape(3, 1)).T[0] for i in vector_gorro]
        return (coord_rotado)

    matriz_R = matrix_R()
    matriz_rotacion = rotation_matrix(matriz_R)
    vector_rotado = rotation_vectors(vecs_center_1, matriz_rotacion)
    vector_rotado_trasladado_a_clique2 = vector_rotado + np.array(bari_2)
    protein_trasladado_rotado = vector_rotado_trasladado_a_clique2
    protein_to_compare = vecs2

    for residuo_not_in_clique, res_1 in zip(protein_trasladado_rotado_sin_clique, atom_numbers1):
        for residuo, res_2 in zip(protein_to_compare_sin_clique, atom_numbers2):
            distancia = pdist(np.array([residuo_not_in_clique, residuo]), metric='euclidean').item()
            if distancia < 3.5:
                residuos_match.append([distancia, (res_1, res_2)])

    df = pd.DataFrame(residuos_match, columns=['distancia', 'candigatos'])
    df = df.sort_values('distancia')
    df['cand1'] = df.candigatos.str[0]
    df['cand2'] = df.candigatos.str[1]

    c1 = []
    c2 = []
    cand_n = []
    for i in df.values:
        if (i[2] in c1) or (i[3] in c2) or (i[0] > 3.5):
            continue
        else:
            c1.append(i[2])
            c2.append(i[3])
            cand_n.append(i)

    df = pd.DataFrame(cand_n, columns=['distancia', 'candigatos', 'cand_1', 'cand2'])
    so_temp = len(cand_n) / number_of_residues

    return idx, so_temp, df, parejas_cliques

timenow_bueno = datetime.datetime.now()

id_so = 0
so = 0.0
df_candidatos = pd.DataFrame()
parejas_clique = []
print(len(df_rmsd))
for i in np.arange(len(df_rmsd)):
    id_temp, so_temp, df_temp, parejas_clique_temp = (align_residues(i, so))
    if so <= so_temp:
        so = so_temp
        id_so = id_temp
        df_candidatos = df_temp
        parejas_clique = parejas_clique_temp
        print(id_so, so)
    if so == 1:
        break

print('ya termine de alinear residuo a residuo')
print(id_so, so)

new_df_cliques1.to_pickle('clique1.pkl')
new_df_cliques2.to_pickle('clique2.pkl')
df_atoms1.to_pickle('clique1_df_atributos.pkl')
df_atoms2.to_pickle('clique2_df_atributos.pkl')
df_candidatos.to_pickle('df_alineados.pkl')
pd.DataFrame(parejas_clique).to_pickle('parejas.pkl')
df_rmsd.to_pickle('df_rmsd.pkl')

print('falta un alineamiento')
exit()
print('imprimiendo pdbs')


# lista_vectores_gorro = []
# for bari in new_df_cliques1.baricentro_clique.values:
#     lista_pre_vectores = [coord - bari for coord in df_atoms1.vector.values]
#
#     lista_vectores_gorro.append(lista_pre_vectores)
#
# vectores_gorro_proteina_1 = pd.DataFrame(lista_vectores_gorro)
#
# lista_vectores_gorro = []
# for bari in new_df_cliques2.baricentro_clique.values:
#     lista_pre_vectores = [coord - bari for coord in df_atoms2.vector.values]
#
#     lista_vectores_gorro.append(lista_pre_vectores)
#
# vectores_gorro_proteina_2 = pd.DataFrame(lista_vectores_gorro)
#
# matriz_rotacion = rmsd.iloc[id_so].matriz_rotacion
#
# vector_gorro = vectores_gorro_proteina_1.iloc[rmsd.iloc[id_so].candidato_clique_1].values
#
# coord_vectores_rotados = [np.matmul(matriz_rotacion, i.reshape(3, 1)).T[0] for i in vector_gorro]
#
# baricentro_proteina_2 = new_df_cliques2.iloc[rmsd.iloc[id_so].candidato_clique_2].baricentro_clique
#
# protein_trasladado_rotado = coord_vectores_rotados + baricentro_proteina_2  # nuevas coordendas proteina 1
#
# new_df_atom1 = pd.concat([df_atoms1, pd.DataFrame(protein_trasladado_rotado, columns=['x', 'y', 'z'])], 1)
# new_df_atom1['new_vector'] = [
#     [new_df_atom1.iloc[i]['x'], new_df_atom1.iloc[i]['y'], new_df_atom1.iloc[i]['z']] for i in range(new_df_atom1.shape[0])]

# for i in pdb11:
#     mask = np.where(i.resi == new_df_atom1.residue_number, True, False)
#     for j in new_df_atom1[mask].atom_name:
#         mask_2 = np.where(new_df_atom1[mask].atom_name == j, True, False)
#         i.GetAtom(j).UpDateValue('coord', new_df_atom1[mask][mask_2].new_vector.values[0])
#
# pdb1.pdbdata = pdb11
# pdb1.WriteToFile(file_out_name='1xxa_el_bueno')

time_bueno = datetime.datetime.now()
print('iteraciones completas:', time_bueno - timenow_bueno)
