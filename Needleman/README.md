# Partie Needleman-Wunsch

Un exemple à été pré-entré

A: "CGGTATTCTTCGAAAAGCGTAGGCGAAGCCCGACGGCATGCGGACATGGGCGAGAGCCTACTTGAGGGTAGCACTACATAGGTTTGATCGTACTGGGGGA"
B: "ACCCGGAAGATGCGGCGGAGAAAACTCGGGTATGCCCCAGACAGTTGAGTGGGTTTCGTATATTCATGCGCCGTGCCAACCGGGGAGATGCTAGGCATGA"

pour tester d'autres chaines, vous pouvez changer les chaines dans `main.cpp`, ~l8 en modifiant les variables SEQ_A, SEQ_B, LENA, LENB

### Compilation + Execution + Cleanup des .o (Recommandé)
```bash
make all
```

### Compilation Seule
```bash
make main
```

ou
```bash
mpic++ -c main.cpp -o main.o -std=c++11 -Wall -Wextra -pedantic -lcgraph -lgvc -lxdot -fopenmp
mpic++ -c needleman_wunsch.cpp -o needleman_wunsch.o -std=c++11 -Wall -Wextra -pedantic -lcgraph -lgvc -lxdot -fopenmp
mpic++ -o main main.o needleman_wunsch.o -std=c++11 -Wall -Wextra -pedantic -lcgraph -lgvc -lxdot -fopenmp
```

### Execution Seule (Pour ne pas tout recompiler à chaque fois)
```bash
make run
```

ou
```bash
./main
```


### Cleanup

Retire uniquement les `.o`
```bash
make cleano
```

Retire les `.o` et `main`
```bash
make clean
```

ou
```bash
rm -f *.o *~ main
```