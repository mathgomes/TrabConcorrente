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
typedef struct {

    int id,J_ORDER;
    int *forBegin;
    int *forEnd;
    double maiorDif,maiorValor;
    double **a;
    double *b;
    double *result;
    double *vetorInicio;
    double *vetAux;

}thread_args;

void *process( void *ptr ) {


    thread_args *args = (thread_args*)ptr;
    int i,j;
    double soma;
    for(i=0;i<args->J_ORDER;++i){
        soma=0;
        for(j=0;j<args->J_ORDER;++j){
            if(i!=j){
                soma += -(args->a[i][j]*args->vetorInicio[j]);
            }
        }
        soma += args->b[i];

        //quarda os valores dos Xi.
        args->result[i] = soma;
        //guarda o maior valor dentre a diferença dos atuais Xi e os da interação anterior.
        if(args->maiorDif < fabs(args->result[i]-args->vetorInicio[i])){
            args->maiorDif = fabs(args->result[i]-args->vetorInicio[i]);
        }
        //guarda o maior valor dentre os atuais Xi.
        if(args->maiorValor < fabs(args->result[i])){
            args->maiorValor = fabs(args->result[i]);
        }
        args->vetAux[i] = args->result[i];
    }
    pthread_exit((void *)NULL);
}

/*Função de leitura das entradas do programa*/
void leitura(int *J_ORDER, int *J_ROW_TEST, double *J_ERROR, int *J_ITE_MAX){
    scanf("%d %d %lf %d", J_ORDER, J_ROW_TEST, J_ERROR, J_ITE_MAX);
}

/*Função de leitura da matriz A*/
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

/*Função de leitura do vetor B*/
double *leituraMb(int J_ORDER){
    double *b;
    int j;

    b = (double*) malloc(J_ORDER*sizeof(double));
    for(j=0;j<J_ORDER;++j){
        scanf("%lf", &b[j]);
    }
    return b;
}

/*Função que cria o vetor usado para guardar os valores iniciais dos Xi*/
double *criaVInicio(double **a, double *b, int J_ORDER){
    double *result;
    int i;

    result = (double*) malloc(J_ORDER*sizeof(double));
    for(i=0;i<J_ORDER;++i){
        result[i] = 0;
    }

    return result;
}

/*Função que executa o método Jacobi-Richardson*/
double *jacobi_richardson(pthread_t *thread, thread_args *args, int numThreads, int J_ORDER, int J_ROW_TEST, int J_ITE_MAX, double J_ERROR){
    int i=0, j=0, converge=0, k=0,return1;
    double soma, res = 0, dif;

    //Criando o vetor que conterá os resultados.
    args->result = (double*) malloc(J_ORDER*sizeof(double));
    args->vetAux = (double*) malloc(J_ORDER*sizeof(double));

    //Verifica se o sistema irá convergir. A variável 'fila' é a linha da matriz em que testará o processo
    for(i=0;i<J_ORDER;++i){
        soma=0;
       for(j=0;j<J_ORDER;++j){
            if(i!=j){
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

    //faz result receber o vetorInicio.
    if(converge==J_ORDER){
        while(k<J_ITE_MAX){
            args->maiorDif = -100000;
            args->maiorValor = -100000;
            for(i =0; i < numThreads; i ++) {
                args->id = i;
                return1 = pthread_create( &thread[i], NULL, process, (void*)args);
                fatal(return1 != 0);

            }

            for(i=0;i<J_ORDER;++i){
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

/*Função que escreve no arquivo os resultados finais dos Xi*/
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

    int numCPU = sysconf(_SC_NPROCESSORS_ONLN);
    pthread_t *thread = (pthread_t*)malloc(numCPU*sizeof(pthread_t));
    thread_args *args = (thread_args*)malloc(sizeof(thread_args));

    double J_ERROR, **a, *b, *vetorInicio, *resultado, tempo;
    int i, j, J_ROW_TEST, J_ORDER, J_ITE_MAX;
    clock_t itime, ftime;

    leitura(&J_ORDER, &J_ROW_TEST, &J_ERROR, &J_ITE_MAX);

    a = leituraMa(J_ORDER);
    b = leituraMb(J_ORDER);

    vetorInicio = criaVInicio(a, b, J_ORDER);

    itime = clock();
    resultado = jacobi_richardson(thread,args,numCPU, J_ORDER, J_ROW_TEST, J_ITE_MAX, J_ERROR);

    ftime = clock();
    tempo = (ftime-itime) / (CLOCKS_PER_SEC * 1.0);

    escreverArquivo(resultado, J_ORDER, tempo);

    return 0;
}
