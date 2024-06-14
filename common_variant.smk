onstart:
    print("##### Genotype-Phenotype Map Pipeline #####")

some_result = 'some_result.txt'

rule all:
    input: some_result

onsuccess:
    print('yay')

onerror:
    print('nooooo')
