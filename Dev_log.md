### Choix d'implémentation

Matrice de départ :
1. Scatter de ma matrice
   1. Envoi ligne par ligne aux PID correspondants  nb_nodes scatter, avec des communications correspondant à chaque groupe de processus

   2. &#9989; Mise en ligne de la donnée permettant l’envoi de la donnée avec un seul scatter 


2. Adaptation de la matrice
   1. Chaque processus adapte son bloc a la reception
   2. &#9989; La preparation d'un bloc demandant la réécriture d'une matrice entière, effectuer l'adaptation au cours de cette étape est plus efficace. (Cela rend la matrice préparer moins lisible).

3. PAM details:
   1. Menoide locaux : On considère que k sera "petit" comparé aux nombre de noeud, la boucle while pour initialisé les candidate (while(candidat[choice != 0])) ne pausera pas de problème. En cas de k plus grand une autre implémentation sera necessaire.
   2. Le choix de déterminer des ménoïdes locaux vient de deux hypothèse. Si le nombre de noeud est suffisament grand par rapport au nombre de processus, et avec une 
   distribution aléatoire dans les scatters des noeuds. Et avec l'hypothèse que les communauté sont de taille comparable, alors le calcule de ménoïde locaux ne devrai pas être trop éloignée des ménoïde globaux. Cela servira de base a l'algorithme global.
   3. Idée pour la "phase d'ébauche, ne pas pouvoir considérer ses propres noeud comme des ménoides"
