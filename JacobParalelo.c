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
#include <pthread.h>
#include <unistd.h>
#include <errno.h>
//Define que checa se uma expressao e verdadeira e encerra a execucao com EXIT_FAILURE
#define fatal(expression)   do {if (expression) {fprintf (stderr, \
                                    "%s : %d : %d\n", __FILE__, __LINE__,  \
                                    strerror(errno)); exit(EXIT_FAILURE) ;}} \
                            while (0)


/**
Estrutura representando os argumentos que cada thread usa para si propria
*/
typedef struct {

    // identificacao
    int idNum;
    int *id;
    // inicio da secao do loop
    int *forBegin;
    // fim da secao do loop
    int *forEnd;

}thread_args;

/**
Estrutura representando os argumentos da funcao utilizada pelas threads
*/
typedef struct {

    int J_ORDER;
    thread_args *Targs;
    double maiorDif,maiorValor;
    double **a;
    double *b;
    double *result;
    double *vetorInicio;
    double *vetAux;

}function_args;

FILE *arq2;
/**
Funcao que sera executava em paralelo pelas threads
*/
void *process( void *ptr ) {

    function_args *args = (function_args*)ptr;
    int i,j;
    double soma;
    // dependendo do numero da thread
    switch(args->Targs->idNum) {

        // executa uma secao do for loop
        case 0 :
            for(i=args->Targs->forBegin[args->Targs->id[0]]; i < args->Targs->forEnd[args->Targs->id[0]]; ++i){
                soma=0;
                for(j=0;j<args->J_ORDER;++j){
                    if(i!=j){
                        soma += -(args->a[i][j]*args->vetorInicio[j]);
                    }
                }
                soma += args->b[i];

                //guarda os valores dos Xi.
                args->result[i] = soma;
            }
            pthread_exit((void *)NULL);
            break;
        case 1:
            for(i=args->Targs->forBegin[args->Targs->id[1]]; i < args->Targs->forEnd[args->Targs->id[1]]; ++i){
                soma=0;
                for(j=0;j<args->J_ORDER;++j){
                    if(i!=j){
                        soma += -(args->a[i][j]*args->vetorInicio[j]);
                    }
                }
                soma += args->b[i];

                //guarda os valores dos Xi.
                args->result[i] = soma;
            }
            pthread_exit((void *)NULL);
            break;
        case 2:
            for(i=args->Targs->forBegin[args->Targs->id[2]]; i < args->Targs->forEnd[args->Targs->id[2]]; ++i){
                soma=0;
                for(j=0;j<args->J_ORDER;++j){
                    if(i!=j){
                        soma += -(args->a[i][j]*args->vetorInicio[j]);
                    }
                }
                soma += args->b[i];

                //guarda os valores dos Xi.
                args->result[i] = soma;
            }
            pthread_exit((void *)NULL);
            break;
        case 3:
            for(i=args->Targs->forBegin[args->Targs->id[3]]; i < args->Targs->forEnd[args->Targs->id[3]]; ++i){
                soma=0;
                for(j=0;j<args->J_ORDER;++j){
                    if(i!=j){
                        soma += -(args->a[i][j]*args->vetorInicio[j]);
                    }
                }
                soma += args->b[i];

                //guarda os valores dos Xi.
                args->result[i] = soma;
            }
            pthread_exit((void *)NULL);
            break;
        case 4:
            for(i=args->Targs->forBegin[args->Targs->id[4]]; i < args->Targs->forEnd[args->Targs->id[4]]; ++i){
                soma=0;
                for(j=0;j<args->J_ORDER;++j){
                    if(i!=j){
                        soma += -(args->a[i][j]*args->vetorInicio[j]);
                    }
                }
                soma += args->b[i];

                //guarda os valores dos Xi.
                args->result[i] = soma;
            }
            pthread_exit((void *)NULL);
            break;
    }
}
/**Função de leitura das entradas do programa

    Args:
    int *J_ORDER : ponteiro para a ordem da matriz
    int *J_ROW_TEST : ponteiro para a coluna teste
    double *J_ERROR : ponteiro para o erro a ser avaliado
    int *J_ITE_MAX : ponteiro para o maximo de iteracoes
*/
void leitura(int *J_ORDER, int *J_ROW_TEST, double *J_ERROR, int *J_ITE_MAX){
    fscanf(arq2,"%d %d %lf %d", J_ORDER, J_ROW_TEST, J_ERROR, J_ITE_MAX);
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
            fscanf(arq2,"%lf", &a[i][j]);

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
        fscanf(arq2,"%lf", &b[j]);
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
    pthread_t thread* : ponteiro para as threads
    funtion_args *args : argumentos da funcao paralela
    int numThreads : numero de threads a ser utilizado
    int J_ORDER : Ordem da matriz
    int J_ROW_TEST : Coluna a ser testada
    int J_ITE_MAX : Maximo de iteracoes
    double J_ERROR : Erro a ser avaliado

    Return:
    Ponteiro para o vetor de solucoes
*/
double *jacobi_richardson(pthread_t *thread, function_args *args, int numThreads, int J_ORDER, int J_ROW_TEST, int J_ITE_MAX, double J_ERROR){
    int i=0, j=0, converge=0, k=0,return1;
    double soma, res = 0, dif;
    numThreads = 5;
    args->J_ORDER = J_ORDER;
    //Criando o vetor que conterá os resultados.
    args->result = (double*) malloc(J_ORDER*sizeof(double));
    args->vetAux = (double*) malloc(J_ORDER*sizeof(double));


    //Verifica se o sistema irá convergir. A variável 'fila' é a linha da matriz em que testará o processo
    for(i=0;i<J_ORDER;++i){
        soma=0;
       for(j=0;j<J_ORDER;++j){
            if(i!=j){
                //printf("%d %d %lf %lf\n",i,j,soma,args->a[i][j]);
                soma = soma + fabs(args->a[i][j]);
                args->a[i][j] = args->a[i][j]/args->a[i][i];
            }

       }
       if(soma < args->a[i][i]){
            converge++;
       }
       args->b[i] = args->b[i]/args->a[i][i];
       args->a[i][i] = 1;
    }


    args->Targs = (thread_args*)malloc(numThreads*sizeof(thread_args));

    args->Targs->id = (int*)malloc(numThreads*sizeof(int));
    args->Targs->forBegin = (int*)malloc(numThreads*sizeof(int));
    args->Targs->forEnd = (int*)malloc(numThreads*sizeof(int));

    args->Targs->forBegin[0] = 0;
    args->Targs->forEnd[0] = args->J_ORDER/numThreads;
    args->Targs->id[0] = 0;

    for(i=1;i < numThreads;i ++) {
        args->Targs->id[i] = i;
        args->Targs->forBegin[i] = args->Targs->forBegin[i-1] + args->J_ORDER/numThreads;
        args->Targs->forEnd[i] = args->Targs->forEnd[i-1] + args->J_ORDER/numThreads;
    }


    if(converge==J_ORDER){
        while(k<J_ITE_MAX){
            args->maiorDif = -100000;
            args->maiorValor = -100000;

            for(i = 0; i < numThreads; ++i ) {
                args->Targs->idNum = i;
                return1 = pthread_create( &thread[i], NULL, process, (void*)args);
            }

            for(i=0; i< numThreads; i++){
                pthread_join(thread[i], NULL);
            }
            for(i=0;i<J_ORDER;++i){
                //guarda o maior valor dentre a diferença dos atuais Xi e os da interação anterior.
                if(args->maiorDif < fabs(args->result[i]-args->vetorInicio[i])){
                    args->maiorDif = fabs(args->result[i]-args->vetorInicio[i]);
                }
                //guarda o maior valor dentre os atuais Xi.
                if(args->maiorValor < fabs(args->result[i])){
                    args->maiorValor = fabs(args->result[i]);
                }
                args->vetAux[i] = args->result[i];
                // trasnfere os valores para a proxima iteracao
                args->vetorInicio[i] = args->vetAux[i];
            }

            //subtrai em modulo o atual valor da posição result pelo valor anterior. Tudo isso divide pelo modulo do atual valor da posicao result.
            dif = args->maiorDif / args->maiorValor;
            //compara para ver se é menor que o erro.
            if(dif <= J_ERROR){
                break;
            }

            k++;
        }
        soma=0;
        for(i=0;i<J_ORDER;++i){
            soma += args->a[J_ROW_TEST][i]*args->result[i];
        }
        printf("----------------------------------\n");
        printf("Iterations: %d\n", k);
        printf("RowTest: %d => [%lf] =? %lf\n", J_ROW_TEST, soma, args->b[J_ROW_TEST]);
        printf("----------------------------------");

    }
    return args->result;

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

    //int numCPU = sysconf(_SC_NPROCESSORS_ONLN);
    int threadNumber = 5;
    pthread_t *thread = (pthread_t*)malloc(threadNumber*sizeof(pthread_t));
    function_args *args = (function_args*)malloc(sizeof(function_args));

    double J_ERROR, **a, *b, *vetorInicio, *resultado, tempo;
    int i, j, J_ROW_TEST, J_ORDER, J_ITE_MAX;
    clock_t itime, ftime;

    arq2 = fopen("matriz250.txt","r");

    leitura(&J_ORDER, &J_ROW_TEST, &J_ERROR, &J_ITE_MAX);

    args->a = leituraMa(J_ORDER);
    args->b = leituraMb(J_ORDER);

    args->vetorInicio = criaVInicio(args->a, args->b, J_ORDER);
    fclose(arq2);
    arq2 = fopen("result.txt","w+");
    itime = clock();
    resultado = jacobi_richardson(thread,args,threadNumber, J_ORDER, J_ROW_TEST, J_ITE_MAX, J_ERROR);

    ftime = clock();
    tempo = (ftime-itime) / (CLOCKS_PER_SEC * 1.0);

    escreverArquivo(resultado, J_ORDER, tempo);

    return 0;
}
