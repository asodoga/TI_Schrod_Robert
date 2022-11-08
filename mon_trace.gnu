reset
set terminal png
set output "potRR.png"  # Nom du fichier de sortie
set title textcolor rgb "red" "Potential"  # Titre du graphique, de couleur rouge
set xlabel "R1"  # Nom de l'axe x
set ylabel textcolor rgb "green" "R2 "  #Nom de l'axe y, de couleur verte
set zlabel  "Pot "  # Nom de l'axe z et repositionnement au-dessus
set pm3d  # Colorisation de la surface
set hidden3d  # Masquage du quadrillage
set isosamples 80,80  # Dimensionnement des entre-axes de la surface
splot "Pot_r.dat" u 1:2:4 with points palette pointsize 2 pointtype 7
 # Création du graphique 3D, avec splot
