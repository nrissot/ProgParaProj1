# Application ARN

Nous avons inclut 4 exemples :
.------------+----------+--------+--------------------.
| Nom        | nb_nodes | nprocs |         fichier.fa |
+------------+----------+--------+--------------------+
|  arn_500_4 |      500 |      4 |  dataset_500seq.fa |
| arn_500_16 |      500 |     16 |  dataset_500seq.fa |
|   arn_2k_4 |     2000 |      4 | dataset_2000seq.fa |
|  arn_2k_16 |     2000 |     16 | dataset_2000seq.fa |
'------------+----------+--------+--------------------'

> **Remarque** :
> le code génére un fichier dot contenant un graphe dont les noeuds sont colorés en fonction des communautées `./output/arn.dot`. 

> **Remarque2** :
> Les warnings viennent de la librairie MPI, ils apparaissent en raison de l'inclusion des flags qui nous ont été recommandés lors des cours de C/C++ de l'an dernier. 

### Compilation + Execution + Cleanup des .o (Recommandé)
```bash
make arn_500_4
make arn_500_16
make arn_2k_4
make arn_2k_16
```
### Execution Seule (Pour ne pas tout recompiler à chaque fois)
```bash
make run_arn_500_4
make run_arn_500_16
make run_arn_2k_4
make run_arn_2k_16
```
ou

```bash
mpirun -np 4 --oversubscribe ./main_arn "exemples/dataset_500seq.fa" 4 500
mpirun -np 16 --oversubscribe ./main_arn "exemples/dataset_500seq.fa" 4 500
mpirun -np 4 --oversubscribe ./main_arn "exemples/dataset_2000seq.fa" 4 2000
mpirun -np 16 --oversubscribe ./main_arn "exemples/dataset_2000seq.fa" 4 2000
```

usage : `mpi_run -np <nprocs> [--oversubscribe] ./main_arn <fichier_seq.fa> <k> <nb_seq>`

### Compilation Seule
```bash
make main_arn
```
ou
```bash
mpic++ -c main_arn.cpp -o main_arn.o -std=c++11 -Wall -Wextra -pedantic -lcgraph -lgvc -lxdot
mpic++ -c ForARN.cpp -o ForARN.o -std=c++11 -Wall -Wextra -pedantic -lcgraph -lgvc -lxdot
mpic++ -c FuncP.cpp -o FuncP.o -std=c++11 -Wall -Wextra -pedantic -lcgraph -lgvc -lxdot
mpic++ -c Utils.cpp -o Utils.o -std=c++11 -Wall -Wextra -pedantic -lcgraph -lgvc -lxdot
mpic++ -o main_arn main_arn.o ForARN.o FuncP.o Utils.o -std=c++11 -Wall -Wextra -pedantic -lcgraph -lgvc -lxdot

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