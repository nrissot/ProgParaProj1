# Partie Application

Nous avons inclut 3 exemples :

```
.------------+----------+--------+--------------------.
| Nom        | nb_nodes | nprocs |         fichier.fa |
+------------+----------+--------+--------------------+
|    make 20 |       20 |      4 |   dataset_20seq.fa |
|   make 500 |      500 |      4 |  dataset_500seq.fa |
|  make 2000 |     2000 |     16 | dataset_2000seq.fa |
'------------+----------+--------+--------------------'
```

### Compilation + Execution + Cleanup des .o (Recommandé)
```bash
make 20
make 500
make 2000
```

### Compilation Seule
```bash
make main
```

ou
```bash
mpic++ -c main.cpp -o main.o -std=c++11 -Wall -Wextra -pedantic -lcgraph -lgvc -lxdot -fopenmp
mpic++ -c needleman_wunsch.cpp -o needleman_wunsch.o -std=c++11 -Wall -Wextra -pedantic -lcgraph -lgvc -lxdot -fopenmp
mpic++ -c ForARN.cpp -o ForARN.o -std=c++11 -Wall -Wextra -pedantic -lcgraph -lgvc -lxdot -fopenmp
mpic++ -c FuncP.cpp -o FuncP.o -std=c++11 -Wall -Wextra -pedantic -lcgraph -lgvc -lxdot -fopenmp
mpic++ -c Utils.cpp -o Utils.o -std=c++11 -Wall -Wextra -pedantic -lcgraph -lgvc -lxdot -fopenmp
mpic++ -o main main.o needleman_wunsch.o ForARN.o FuncP.o Utils.o -std=c++11 -Wall -Wextra -pedantic -lcgraph -lgvc -lxdot -fopenmp
```

### Execution Seule (Pour ne pas tout recompiler à chaque fois)
```bash
make run_20
make run_500
make run_2000
```

ou
```bash
mpirun -np <1|4|16> <--oversubscribe> ./main "arn/dataset_<20|500|2000>seq.fa" <nb_communauté> <20|500|2000>
```


### Cleanup

Retire uniquement les `.o`
```bash
make cleanox
```

Retire les `.o` et `main`
```bash
make clean
```

ou
```bash
rm -f *.o *~ main
```