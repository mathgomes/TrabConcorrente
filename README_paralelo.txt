-------INSTRUÇÕES DE USO DO PROGRAMA PARALELO--------------------------------------------------------------------
	Para executar o programa sequencial.c, devemos inserir um arquivo de entrada contendo 
todos os dados referentes ao sistema linear a ser resolvido.
	No Windows, devemos ir até o local das pasta dos arquivos e colocar no prompt de comando: 
	sequencial.exe < "nome_do_arquivo".txt
	No Linux, devemos ir até o local da pasta dos arquivos e colocar no prompt de comando:
	gcc -o <nome_do_executavel> <nome_do_arquivo.c>
	.\<nome_do_executavel> < <nome_do_arquivo_de_matriz>.txt -lpthreads

	Na ordem de entrada, temos:
	1. A ordem do sistema linear (valor inteiro);
	2. A linha de comparação dos resultados (valor inteiro);
	3. O erro mínimo para a solução do sistema linear (valor decimal);
	4. O número de iterações máximas do método (valor inteiro).
	5. Matriz A do sistema linear;
	6. Vetor B do sistema linear.

Após executado o programa com as entradas apresentadas, será finalizado mostrando as 
segintes informações:
	1. Iterations: %quantidade_da_iteracões
	2. RowTest: %linha_de_comparacao => [%valor_aproximado_da_matriz_na_linha] =? %valor_da_matriz_B_na_linha

Também será gerado um arquivo com o nome "resultado.txt", que conterá as seguintes informações:
	1. Ordem da matriz utilizada
	2. Tempo de execução do algoritmo
	3. Valores finais dos Xi
------------------------------------------------------------------------------------------------------------------
