# -*- coding: utf-8 -*-
# Auteurs : Alexis Neven & Valentin Dall'Alba

import numpy as np


##############################
# Fonctions d’écriture des fichiers d’entrée CFP
##############################

def write_geoheight(newX, path_dir):
    """
    Écrit le fichier geoheight.txt.

    Paramètres
    ----------
    newX : list
        Liste des coordonnées de nœuds.
    path_dir : str
        Répertoire de sortie.
    """
    with open(f"{path_dir}/geoheight.txt", 'w') as f:
        for i in range(len(newX)):
            endline = "" if i == len(newX) - 1 else "\n"
            f.write(f"{i+1}\t-1{endline}")


def write_kexch(newX, path_dir):
    """
    Écrit le fichier K_exch.txt (valeur par défaut = 5 pour chaque nœud).

    Paramètres
    ----------
    newX : list
        Liste des coordonnées de nœuds.
    path_dir : str
        Répertoire de sortie.
    """
    with open(f"{path_dir}/K_exch.txt", 'w') as f:
        for i in range(len(newX)):
            endline = "" if i == len(newX) - 1 else "\n"
            f.write(f"{i+1}\t5{endline}")


def write_node_heads(newX, path_dir, node_out=99):
    """
    Écrit le fichier node_head.txt.

    Paramètres
    ----------
    newX : list
        Liste des coordonnées de nœuds.
    path_dir : str
        Répertoire de sortie.
    node_out : int
        Index du nœud de sortie (valeur de tête fixée à 0).
    """
    with open(f"{path_dir}/node_head.txt", 'w') as f:
        for i in range(len(newX)):
            node_id = i + 1
            if node_id == node_out and node_out != len(newX):
                f.write(f"{node_id}\t0\n")
            elif node_id == node_out:
                f.write(f"{node_id}\t0")
            elif node_id == len(newX):
                f.write(f"{node_id}\t-1")
            else:
                f.write(f"{node_id}\t-1\n")


def write_pipe_info(newSeg, path_dir, tortuosity=1,
                    rugosity=0.01, Remin=2000, Remax=3000):
    """
    Écrit le fichier pipe_info.txt.

    Paramètres
    ----------
    newSeg : list
        Liste des segments (chaque segment = [id_start, id_end, diamètre]).
    path_dir : str
        Répertoire de sortie.
    tortuosity : float
        Tortuosité du conduit.
    rugosity : float
        Rugosité (sans dimension).
    Remin, Remax : int
        Bornes de Reynolds.
    """
    with open(f"{path_dir}/pipe_info.txt", 'w') as f:
        for i, seg in enumerate(newSeg):
            diam = seg[2]
            endline = "" if i + 1 == len(newSeg) else "\n"
            f.write(f"{i+1}\t{diam}\t{tortuosity}\t{rugosity}\t{Remin}\t{Remax}{endline}")


def write_network_info(newX, newsegs, path_dir, dx, dy):
    """
    Écrit le fichier network_info.txt avec les connexions entre nœuds.

    Paramètres
    ----------
    newX : list
        Liste des coordonnées des nœuds (x, y).
    newsegs : list
        Liste des segments (liens entre nœuds).
    path_dir : str
        Répertoire de sortie.
    dx, dy : float
        Tailles de maille en x et y.
    """
    with open(f"{path_dir}/network_info.txt", 'w') as f:
        for i, node in enumerate(newX):
            X = max(1, int(node[0] / dx) + 1)
            Y = max(1, int(node[1] / dy) + 1)

            id_connect, pipe_connect = find_voisin(newsegs, i)

            line = f"{i+1}\t{X}\t{Y}\t1\t" + \
                   "\t".join(map(str, id_connect)) + "\t" + \
                   "\t".join(map(str, pipe_connect))

            endline = "" if i + 1 == len(newX) else "\n"
            f.write(line + endline)


def find_voisin(newsegs, id_node):
    """
    Trouve les voisins d’un nœud et les conduits qui y sont connectés.

    Paramètres
    ----------
    newsegs : list
        Liste des segments (id_start, id_end, diamètre).
    id_node : int
        Index du nœud courant.

    Retour
    ------
    id_connect : list
        Liste des indices de nœuds connectés (complétée à 6 voisins max).
    pipe_connect : list
        Liste des indices de segments connectés (complétée à 6).
    """
    id_connect = []
    pipe_connect = []

    for i, seg in enumerate(newsegs):
        if seg[0] == id_node:
            id_connect.append(int(seg[1] + 1))
            pipe_connect.append(i + 1)
        elif seg[1] == id_node:
            id_connect.append(int(seg[0] + 1))
            pipe_connect.append(i + 1)

    # Compléter jusqu’à 6 entrées
    while len(id_connect) < 6:
        id_connect.append(0)
    while len(pipe_connect) < 6:
        pipe_connect.append(0)

    return id_connect, pipe_connect
