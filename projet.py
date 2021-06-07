import numpy as np
import random 
import matplotlib.pyplot as plt

def stochastic(matrix):
    """ Fonction qui retourne True si la matrice est stochastique
    Argument : matrice
    Retroune : True si matrice est une matrice carrée dont chaque élément est un réel positif 
    et dont la somme des éléments de chaque ligne vaut 1.
    Elle retourne False sinon.
    """
    num_rows, num_cols = matrix.shape

    if (num_rows != num_cols):
        return False

    for i in range(0,num_rows):
        sum=0
        for j in range(0,num_cols):
            if(matrix[i][j]<0) :
                print("Valeur dans mat[{}][{}] est negative".format(i,j))
                return False
            sum+=matrix[i][j]
        if sum != 1 :
            print("Somme dans la ligne {} n'est pas egale a 1".format(i))
            return False
    return True

def tirage_aleatoire_un_etat(etat,matrice):
    """Fonction qui tire aleatoirement un etat, a partir d'un etat actuel donné en parametre 
    Elle recupere les etats qu'on peut atteindre a partir de cet etat dans la matrice de transition.
    """
    #cas ou c'est l'etat initial :
    if(etat==-1):
        raw=matrice # on donne en parametre le vecteur initial
    else :
        raw=matrice[etat]
    
    return random.choices([0,1,2],raw)[0]



def tirage_aleatoire(temps,etat_initial,matrice):
    """ Fonction qui tire aléatoirement une sequence de S,I,R pour un individu
    Argument : $1 : temps , $2 vecteur initial , $3 matrice des transitions de la ch de Markov
    Retroune : une chaine de caractere.
    """

    #etat initial
    init=tirage_aleatoire_un_etat(-1,etat_initial)
    
    resultat=""

    for i in range(0,temps):
        resultat+=str(init)+","
        init=tirage_aleatoire_un_etat(init,matrice)

    resultat= resultat[:-1]
    resultat+="\n"
    return resultat



def modelisation_population(matrice,vecteur_initial,individus,temps,nom_fichier):
    """ Fonction qui rempli le fichier donné en parametre en tournant le tirage sur tous les individus
    Argument :$1 matrice de transition, $2 vecteur initial $3 nombre d'individus
    $4 temps , $5 nom du fichier
    """
    outputFile = open(nom_fichier, "w")
    s=""
    for i in range(0,temps):
        s+=str(i)+","
    s=s[:-1]
    outputFile.write(s+"\n")

    for i in range(0,individus):
        outputFile.write(tirage_aleatoire(temps,vecteur_initial,matrice))



def nombre_sains_infectes_instant_t(df,instant):
    """ Fonction qui retourne le nombre d'individus sains et infectes a un instant t
    Argument : $1: dataframe, $2: instant
    Retroune : un dictionnaire contenant pour chaque cle qui est l'etat d'un individu, une valeur qui correspond au nombre
    de personnes dans cet etat a l'instant t donne
    """
    d= df.to_dict('list')
    dic= dict((str(k), v) for k, v in d.items())
    column=dic[instant]
    dic=[]
    dic.append(column.count(0))
    dic.append(column.count(1))
    dic.append(column.count(2))
    return dic
    
    
    

def nombre_sains_infectes_tous_instants(df):
    """ Fonction qui retourne le nombre d'individus sains et infectes pendant la simulation
    Argument : $1: dataframe,
    Retroune : une liste contenant le nombre de personnes infecteees et saines a tout instant t
    dic[0] contient [nb infectes t=0, nb sains t=0, nb retablis t=0]
    """
    total_cols=len(df.axes[1])
    dic=[]
    for i in range(0,total_cols):
        dic.append(nombre_sains_infectes_instant_t(df,str(i)))
    return dic

def draw_graph(df,nom_graphe,labels):
    x = range(0,len(df.axes[1]))
    y = nombre_sains_infectes_tous_instants(df)
    plt.xlabel("Temps")
    plt.ylabel("Nombre de personnes dans chaque categorie")
    plt.title(nom_graphe)
    for i in range(len(y[0])):
        plt.plot(x,[pt[i] for pt in y],label = labels[i])
    plt.legend()
    plt.show()

def pic_epidemie(df):
    """Fonction qui retourne le moment et le nombre d'infectes au moment du pic de l'epidemie.
    Cela correspond au moment ou il y a le plus de personnes infectées pendant l'observation
    On peut aussi le voir sur le graphe.
    """

    total_cols=len(df.axes[1])
    liste=[]
    max_infecte=0
    max_temps=0
    for i in range(0,total_cols):
        nb_infectes = nombre_sains_infectes_instant_t(df,str(i))
        liste.append(nb_infectes[1])
        if(nb_infectes[1] > max_infecte):
            max_infecte=nb_infectes[1]
            max_temps=i
    print("Le pic de l'epidemie a lieu a l'instant {} avec un nombre total de personnes infectées : {}".format(max_temps,max_infecte))

def longueur_infection_un_individu(ligne_df):
    """Fonction qui calcule la longueur d'infection d'un individu, et donc 
    le temps de guerison """
    ligne=str(ligne_df)
    nb_infecte=0
    for etat in ligne :
        if (etat =="1"):
            nb_infecte+=1
    return nb_infecte

def longueur_infection_individus(df):
    """Fonction qui calcule la moyenne de la longueur d'infection des individus
    """
    file = open(df)
    lines = file.readlines()
    len_infection_totale=0
    for index in range(1,len(lines)):
        line = lines[index].strip()
        len_infection_totale += longueur_infection_un_individu(line)

    return len_infection_totale/(len(lines)-1)


##---------------Modele distanciation sociale----------------------- 

def alternance_periodes(matrice_i,vecteur_initial,matrice_s,individus,temps):
    """Fonction pour permettre d'alterner entre les periodes de non distanciation 
    et de distanciation
    Parametre : $1 matrice initiale : matrice des transitions modele sans distanciation
    $2 vecteur initial, $3 matrice secondaire : apres confinement : proba de devenir infectee <<
    $4 nb individus , $5 temps

    Retourne : un dictionnaire donc les cles sont les valeurs du temps,
    valeurs : liste pour les etat de toutes les personnes a ce temps la.
    Nous avons modifié par rapport aux parties precedentes car ici c'est mieux d'avoir les personnes en colonne
    et difficultê d'ecrire dans le fichier.
    """
    dic=dict()
    confinement=False
    matrice=matrice_i
    
    #initialisation des etats initiaux de tous les individus
    dic[0]=[]
    for i in range(0,individus):
        dic[0].append(tirage_aleatoire_un_etat(-1,vecteur_initial))
    
    for i in range(1,temps):
        dic[i]=[]
        nb_infecte=dic[i-1].count(1)
       

        if ((nb_infecte/individus)>=0.3):
            matrice=matrice_s
            confinement=True

        elif((nb_infecte/individus)<=0.15 and confinement):
            matrice=matrice_i

        for j in range(0,individus):
            #on recupere l'etat de la personne au temps precedent 
            liste_temps=dic[i-1]
            etat_temps_precedent=liste_temps[j]
            dic[i].append(tirage_aleatoire_un_etat(etat_temps_precedent,matrice))
    return dic


    

