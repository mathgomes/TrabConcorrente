/************************************************************
* Grupo_27: Matheus Gomes ; Renato José.                    *
* Nº USP: 8532321 ; 8531790.                                *
* Início: 24/08/15                                          *
* Última alteração: 03/09/15                                *
* Função: Método Jacobi-Richardson Sequencial               *
************************************************************/

#include <stdbool.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/**Função de leitura das entradas do programa

    Args:
    int *J_ORDER : ponteiro para a ordem da matriz
    int *J_ROW_TEST : ponteiro para a coluna teste
    double *J_ERROR : ponteiro para o erro a ser avaliado
    int *J_ITE_MAX : ponteiro para o maximo de iteracoes
*/
void leitura(int *J_ORDER, int *J_ROW_TEST, double *J_ERROR, int *J_ITE_MAX){
    scanf("%d %d %lf %d", J_ORDER, J_ROW_TEST, J_ERROR, J_ITE_MAX);
}

/**Função de leitura da matriz A

    Args:
    int J_ORDER : ordem da matriz
    Return:
    Ponteiro de ponteiro para a matriz lida
*/
double **leituraMa(int J_ORDER){
    double **a;
    int i, j;

    a = (double**) malloc(J_ORDER*sizeof(double*));
    for(i=0;i<J_ORDER;++i){
        a[i] = (double*) malloc(J_ORDER*sizeof(double));
        for(j=0;j<J_ORDER;++j){
            scanf("%lf", &a[i][j]);

        }
    }
    return a;
}

/**Função de leitura da matriz B

    Args:
    int J_ORDER : ordem da matriz
    Return:
    Ponteiro para a matriz lida
*/
double *leituraMb(int J_ORDER){
    double *b;
    int j;

    b = (double*) malloc(J_ORDER*sizeof(double));
    for(j=0;j<J_ORDER;++j){
        scanf("%lf", &b[j]);
    }
    return b;
}

/**Função que cria o vetor usado para guardar os valores iniciais dos Xi

    Args:
    double **a : Ponteiro para a matriz A
    double *b : Ponteiro para a matriz B
    int J_ORDER : ordem da matriz

    Return:
    ponteiro para o vetor inicial

*/
double *criaVInicio(double **a, double *b, int J_ORDER){
    double *result;
    int i;

    result = (double*) malloc(J_ORDER*sizeof(double));
    for(i=0;i<J_ORDER;++i){
        result[i] = 0;
    }

    return result;
}

/**Função que executa o método Jacobi-Richardson

    Args:
    double **a : ponteiro para matriz A
    double *b : ponteiro para matriz B
    double *vetorInicio : ponteiro para o vetor inicial
    int J_ORDER : Ordem da matriz
    int J_ROW_TEST : Coluna a ser testada
    int J_ITE_MAX : Maximo de iteracoes
    double J_ERROR : Erro a ser avaliado

    Return:
    Ponteiro para o vetor de solucoes
*/
double *jacobi_richardson(double **a, double *b, double *vetorInicio, int J_ORDER, int J_ROW_TEST, int J_ITE_MAX, double J_ERROR){
    int i=0, j=0, converge=0, k=0;
    double soma, res = 0, maiorDif, maiorValor, dif, *result, *vetAux;

    //Criando o vetor que conterá os resultados.
    result = (double*) malloc(J_ORDER*sizeof(double));
    vetAux = (double*) malloc(J_ORDER*sizeof(double));

    //Verifica se o sistema irá convergir. A variável 'fila' é a linha da matriz em que testará o processo
    for(i=0;i<J_ORDER;++i){
        soma=0;
       for(j=0;j<J_ORDER;++j){
            if(i!=j){
                soma = soma + fabs(a[i][j]);
                a[i][j] = a[i][j]/a[i][i];
            }

       }
       if(soma < a[i][i]){
            converge++;
       }
       b[i] = b[i]/a[i][i];
       a[i][i] = 1;
    }

    //faz result receber o vetorInicio.
    if(converge==J_ORDER){
        while(k<J_ITE_MAX){
            maiorDif = -100000;
            maiorValor = -100000;
            for(i=0;i<J_ORDER;++i){
                soma=0;
                for(j=0;j<J_ORDER;++j){
                    if(i!=j){
                        soma += -(a[i][j]*vetorInicio[j]);
                    }
                }
                soma += b[i];

                //quarda os valores dos Xi.
                result[i] = soma;
                //guarda o maior valor dentre a diferença dos atuais Xi e os da interação anterior.
                if(maiorDif < fabs(result[i]-vetorInicio[i])){
                    maiorDif = fabs(result[i]-vetorInicio[i]);
                }
                //guarda o maior valor dentre os atuais Xi.
                if(maiorValor < fabs(result[i])){
                    maiorValor = fabs(result[i]);
                }
                vetAux[i] = result[i];
            }
            for(i=0;i<J_ORDER;++i){
                vetorInicio[i] = vetAux[i];
            }

            //subtrai em modulo o atual valor da posição result pelo valor anterior. Tudo isso divide pelo modulo do atual valor da posicao result.
            dif = maiorDif / maiorValor;

            //compara para ver se é menor que o erro.
            if(dif <= J_ERROR){
                break;
            }

            k++;
        }
        soma=0;
        for(i=0;i<J_ORDER;++i){
            soma += a[J_ROW_TEST][i]*result[i];
        }
        printf("----------------------------------\n");
        printf("Iterations: %d\n", k);
        printf("RowTest: %d => [%lf] =? %lf\n", J_ROW_TEST, soma, b[J_ROW_TEST]);
        printf("----------------------------------");

    }

    return result;

}

/**Função que escreve no arquivo os resultados finais dos Xi

    Args:
    double *resultado : Ponteiro para o vetor de solucoes
    int J_ORDER : Ordem da matriz
    double tempo : Tempo calculado da execucao do programa
*/
void escreverArquivo(double *resultado, int J_ORDER, double tempo){
    int i;
    FILE *arq;

    arq = fopen("resultado.txt","w+");
    fprintf(arq, "Tempo de execução: %lf\nOrdem do sistema linear: %d \n\n", tempo, J_ORDER);

    for(i=0;i<J_ORDER;++i){
        fprintf(arq, "X%d: %lf\n", i, resultado[i]);
    }

    fclose(arq);
}

int main(){
    double J_ERROR, **a, *b, *vetorInicio, *resultado, tempo;
    int i, j, J_ROW_TEST, J_ORDER, J_ITE_MAX;
    clock_t itime, ftime;

    leitura(&J_ORDER, &J_ROW_TEST, &J_ERROR, &J_ITE_MAX);

    a = leituraMa(J_ORDER);
    b = leituraMb(J_ORDER);

    vetorInicio = criaVInicio(a, b, J_ORDER);

    itime = clock();
    resultado = jacobi_richardson(a, b, vetorInicio, J_ORDER, J_ROW_TEST, J_ITE_MAX, J_ERROR);

    ftime = clock();
    tempo = (ftime-itime) / (CLOCKS_PER_SEC * 1.0);
    printf("\nTempo : %lf\n",tempo);

    escreverArquivo(resultado, J_ORDER, tempo);

    return 0;
}
