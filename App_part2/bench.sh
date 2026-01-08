#!/bin/bash

# "> 2k : 16 proc --oversubscribe" 
echo "> 2k : 16 proc --oversubscribe"
echo "______________________________________"
for i in {0..10}
do
    mpirun -np 16 --oversubscribe ./main "arn/dataset_2000seq.fa" 4 2000
    echo "______________________________________"
done
# "> 2k : 4 proc" 
echo "> 2k : 4 proc"
echo "______________________________________"
for i in {0..10}
do
    mpirun -np 4 --oversubscribe ./main "arn/dataset_2000seq.fa" 4 2000
    echo "______________________________________"
done
# "> 2k : 1 proc"
echo "> 2k : 1 proc"
echo "______________________________________"
for i in {0..10}
do
    mpirun -np 1 ./main "arn/dataset_2000seq.fa" 4 2000
    echo "______________________________________"
done
# "> 500 : 16 proc --oversubscribe"
echo "> 500 : 16 proc --oversubscribe"
echo "______________________________________"
for i in {0..10}
do
    mpirun -np 16 --oversubscribe ./main "arn/dataset_500seq.fa" 4 500
    echo "______________________________________"
done
# "> 500 : 4 proc"
echo "> 500 : 4 proc"
echo "______________________________________"
for i in {0..10}
do
    mpirun -np 4 --oversubscribe ./main "arn/dataset_500seq.fa" 4 500
    echo "______________________________________"
done
# "> 500 : 1 proc"
echo "> 500 : 1 proc"
echo "______________________________________"
for i in {0..10}
do
    mpirun -np 1 --oversubscribe ./main "arn/dataset_500seq.fa" 4 500
    echo "______________________________________"
done

echo "> EOF"

