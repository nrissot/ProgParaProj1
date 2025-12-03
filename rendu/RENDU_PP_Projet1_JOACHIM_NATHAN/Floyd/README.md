# Partie Floyd

Nous avons inclut trois tailles d'exemples :
.--------+----------+--------+-------------------.
| Taille | nb_nodes | nprocs |       fichier.dot |
+--------+----------+--------+-------------------+
| small  |        8 |      4 |  PetitExemple.dot |
| medium |       20 |      4 |  MoyenExemple.dot |
| big    |      100 |     16 |   GrosExemple.dot |
'--------+----------+--------+-------------------'

> **Remarque**
> Les warnings viennent de la librairie MPI, ils apparaissent en raison de l'inclusion des flags qui nous ont été recommandés lors des cours de C/C++ de l'an dernier. 

### Compilation + Execution + Cleanup des .o (Recommandé)
```bash
make small
make medium
make big
```
### Execution Seule (Pour ne pas tout recompiler à chaque fois)
```bash
make run_small
make run_medium
make run_big
```
ou

```bash
mpirun -np 4 --oversubscribe ./main "exemples/PetitExemple.dot"
mpirun -np 4 --oversubscribe ./main "exemples/MoyenExemple.dot"
mpirun -np 16 --oversubscribe ./main "exemples/GrosExemple.dot"
```

### Compilation Seule
```bash
make main
```
ou
```bash
mpic++ -c main.cpp -o main.o -std=c++11 -Wall -Wextra -pedantic -lcgraph -lgvc -lxdot
mpic++ -c Utils.cpp -o Utils.o -std=c++11 -Wall -Wextra -pedantic -lcgraph -lgvc -lxdot
mpic++ -c FuncP.cpp -o FuncP.o -std=c++11 -Wall -Wextra -pedantic -lcgraph -lgvc -lxdot
mpic++ -c ForGraph.cpp -o ForGraph.o -std=c++11 -Wall -Wextra -pedantic -lcgraph -lgvc -lxdot
mpic++ -o main main.o ForGraph.o FuncP.o Utils.o -std=c++11 -Wall -Wextra -pedantic -lcgraph -lgvc -lxdot
```

### Cleanup

Uniquement les `.o`:
```bash
make cleano
```

`.o` & executable `main`:
```bash
make clean
```