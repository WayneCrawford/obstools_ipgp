*******************************
plotPSDs
*******************************

#################
Introduction
#################

La logiciel « plotPSDs » crée des Probabilistic Power Spectral Densitities
(PPSDs) pour chaque station-voie et les comparaisons des médian Power
Spectral Densities (med_PSDs) pour chaque nom de voie.
Les med_PSDs permet d’identifier rapidement s’il y a un problème avec un ou
plusieurs voies, tandis que les PPSDs permettent de regarder chaque voie en
détail.
Ce qu’il faut regarder c’est :
1. Si on voit le pique « microsismique » vers 10 seconds (présent dans presque
   toutes les données sismologiques
2. Si cette pique est d’à peu près la même taille sur chaque station.

Par default, « plotPSDs « fait ses calculs sur 1 jour, à partir d’un jour
après le début des données de l’instrument (selon ``<Channel><Start_Date>``
dans le fichier StationXML).
Ça suffit pour une validation rapide, mais c’est préférable d’utiliser
plusieurs jours, voir tout le déploiement pour la validation finale.

Le help (y compris une liste des options) se trouve en tapant :

.. code-block::
	plotPSDs -h

La ligne de commande pour créer les images avec les options standards est :

.. code-block::
	plotPSDs <SDS-dir> <StationXML-file>
 
Une fois le programme tourné, vous aller trouver des répertoires plot_PPSDs/
et plot_medPSDs/.
Voici un exemple des images dans plot_medPSDs/:

.. figure:: images/plotPPSDs_medPSDs_BDH.jpg

   Médian PSDs pour la voie d’hydrophone (BDH), avec les limites « standards »
   indiqués par les lignes pointillés.
   Rien à signaler, à part quelques pics à gauche : pour la reste, on voit bien
   le pic des « microséismes » vers 10 secondes et les courbes se superposent.


.. figure:: images/plotPPSDs_medPSDs_SH3.jpg

    Médian PSDs pour la voie géophone verticale (SH3).
    On voit que certaines voies n’ont pas fonctionnés (pas de pic de microséismes,
    lignes bien en-dessous les autres) et qu’il y a une fort « bipolarité » aux
    longues périodes qui pourrais lui-même masquer le pic des microséismes.

On constate que certaines voies verticales n’ont pas fonctionnés.
On note aussi que, dans certaines bande de périodes (basses pour les hydrophones,
hautes pour les géophones) les med_PSDs de certaines stations oscillent entre
deux limites.

Voici le même chose pour des BBOBS :
 	 
.. figure:: images/plotPPSDs_medPSDs_BDG.jpg

    Médian PSDs pour la voie de jauge de pression différentielle (BDG).
    Un des jauges (sur J6) n’a pas fonctionné.
    On voit aussi que le niveau du bruit est plus faible entre 10 et 100
    secondes que sur les hydrophones
    
    
.. figure:: images/plotPPSDs_medPSDs_BHZ.jpg

    Médian PSDs pour la voie sismometre verticale (BHZ).
    On voit que le niveau du bruit est bien inférieur aux géophones entre
    ~1 et 500 secondes.

Pour regarder en plus de détail ces phénomènes, on regarde les PPSDs
individuelles, dans le répertoire « plot_PPSDs » .
Les PPSDs indiquent la probabilité qu’un spectre sera d’un certain niveau,
tandis que le médian donne simplement la valeur pour lequel 50% des spectres
dont plus haut et 50% plus bas.
Le PPSD fourni donc plus de détails, mais ne permet pas de superposer
plusieurs voies.  Voici quelques exemples :

 	 
.. figure:: images/plotPPSDs_medPSDs_BDH.jpg

    PPSD pour la station-voie K6.BDH.
    On voit qu’aux basse périodes (hautes fréquences) il y a bien deux
    différents niveaux « typiques » du bruit : peut-être liée à la météo
    locale ?
    
    
    
.. figure:: images/plotPPSDs_medPSDs_SH3.jpg

    PPSD pour la station-voie MM1.SH3.  
    On voie une nette bifurcation des spectres sur la voie SH3,
    dont le plus haut est au-dessous du pic des microséismes.
    Il s’agit du bruit du disque dur, ce qui ne gêne pas trop la plupart de
    l’enregistrement parce qu’il ne dure que 1 minute sur 60, mais dans les
    longs spectres comme celle-ci, ce bruit peut se trouver dans la moitié
    des fenêtres.
    C’est pour ça que le niveau du bruit « haut » est aussi « probable »
    que le niveau du bruit bas dans ces spectres, et que les med_PSDs peuvent
    osciller entre ces deux valeurs.
    Cet exemple est pour l’expérience SISMANTILLES, ou nous utilisions des
    disques durs mecaniques : le bruit liée à l’activité de la disque est *à
    priori* plus faible maintenant
