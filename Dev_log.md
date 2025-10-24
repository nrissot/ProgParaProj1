### Choix d'implémentation

Matrice de départ :
1. Scatter de ma matrice
   1. Envoi ligne par ligne aux PID correspondants  nb_nodes scatter, avec des communications correspondant à chaque groupe de processus

   2. &#9989; Mise en ligne de la donnée permettant l’envoi de la donnée avec un seul scatter 


2. Adaptation de la matrice
    1. Chaque processus adapte son bloc a la reception
    2. &#9989; La preparation d'un bloc demandant la réécriture d'une matrice entière, effectuer l'adaptation au cours de cette étape est plus efficace. (Cela rend la matrice préparer moins lisible).