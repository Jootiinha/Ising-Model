#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#define L 100
#include <gsl/gsl_rng.h>
using namespace std;

ofstream DATA("DATA.1.dat",ios::out);

// Malha de spins
struct malha{
    int x;
    int y;
};


//Total de spins
const int n = L*L;
//Temperatura final
double tempFinal = 5; //T
//Temperatura inicial
const float tempInicial = 0.5; //Tmin
//Passo para temperatura no loop de temperatura
const float dT = 0.1; //dT
//Rede bidim
int S[L][L];
//Número de passos p/ Monte Carlo
long unsigned int MCSteps = 10000; //MS
//Número de passos descartados para evitar a parte fora do equilíbrio.
int transiente = 1000; 
//Seed
const gsl_rng *r;


//Inicializar a rede
void inicializarRede(){
    for(int y=0;y<L;y++){ 
        for(int x=0;x<L;x++){
            if(gsl_rng_uniform(r)>0.5)
                S[x][y]=1;
            else
                S[x][y]=-1;
        }
    }
}
//Escolha aleatória do sítio
void spinAleatorio(malha &pos){
    pos.x=(gsl_rng_uniform(r)*(L));// (0,L-1)
    pos.y=(gsl_rng_uniform(r)*(L));
}
//Calculo da energia associada a posição
int energiaPosicao(malha &pos){
    //condições periódicas de contorno
    int up, down, left, rigth,e;
    rigth=(pos.x+1)%L;
    left=(pos.x+L-1)%L;
    up=(pos.y+1)%L;
    down=(pos.y+L-1)%L;
    e=-1*S[pos.x][pos.y]*(S[left][pos.y]+S[rigth][pos.y]+S[pos.x][up]+S[pos.x][down]);//J=1
    return e;
}
//função decisão
bool testeFlip(malha pos, int &de){
    de=-2*energiaPosicao(pos);
    if(de<0)
        return true; //flip pela energia
    else{
        if(gsl_rng_uniform(r)<exp(-de/tempFinal))
            return true; //flip pelo banho térmico
        else
            return false; //não flip
    }
}
//função de flip
void flip(malha pos){
    S[pos.x][pos.y]=-S[pos.x][pos.y];
}
//função que remove o transiente
void transienteFinal(){
    malha pos;
    int de=0;
    for(int a=1;a<=transiente;a++){
        for(int b=1;b<=n;b++){
            spinAleatorio(pos);
            if(testeFlip(pos,de)){
                flip(pos);
            }
        }
    }
}
//Magnetização total da rede (soma dos spins)
int magnetizacao(){
    int m=0;
    for(int y=L-1;y>=0;y--){
        for(int x=0;x<L;x++){
            m+=S[x][y];
        }
    }
    return m;
}
//Energia total da rede
int energiaTotal(){
    malha pos;
    int e=0;
    for (int y=L-1;y>=0;y--){
        pos.y=y;
        for(int x=0;x<L;x++){
            pos.x=x;
            e+=energiaPosicao(pos);
        }
    }
    return e;
}
//Programa
int main(){
    //variáveis para o cálculo dos observáveis
    r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r, 158235);
    double E=0,Esq=0,Esq_avg=0,E_avg=0,etotal=0, etotalsq=0;
    double M=0, Msq=0, Msq_avg=0,M_avg=0, mtotal=0, mtotalsq=0;
    double Mabs=0,Mabs_avg=0,Mq_avg=0, mabstotal=0, mqtotal=0;
    int de=0;
    malha pos;
    //inicia a rede de maneira aleatória
    inicializarRede();
    //loop de temperatura
    for(;tempFinal>=tempInicial;tempFinal-=dT){
        //transiente
        transienteFinal();
        //observáveis com valores no equilíbrio (após o transiente)
        M=magnetizacao();
        // Mabs=fabs(magnetizacao());
        E=energiaTotal();
        //reseta as variáveis para cada loop
        etotal=0;
        etotalsq=0;
        mtotal=0;
        mtotalsq=0;
        mabstotal=0;
        mqtotal=0;
        //loop do Monte Carlo (MS=tempo)
        for(int a=1;a<=MCSteps;a++){
            //Metropolis
            for(int b=1;b<=n;b++){
                spinAleatorio(pos);
                if(testeFlip(pos,de)){
                    flip(pos);
                    //ajusta os observáveis
                    E+=2*de;
                    M+=2*S[pos.x][pos.y];
                    //Mabs+=abs(S[pos.x][pos.y]);
                }
            }
            //soma dos obsevavéis
            etotal+=E/2.0;
            etotalsq+=(E/2.0)*(E/2.0);
            mtotal+=M;
            mtotalsq+=M*M;
            mqtotal+=M*M*M*M;
            mabstotal+=sqrt(M*M);
        }
        //Médias:
        E_avg=etotal/(MCSteps*n);
        Esq_avg=etotalsq/(MCSteps*n);
        M_avg=mtotal/(MCSteps*n);
        Mabs_avg=mabstotal/(MCSteps*n);
        Msq_avg=mtotalsq/(MCSteps*n);
        Mq_avg=mqtotal/(MCSteps*n);

        //saída para arquivo
        DATA<<tempFinal<< //temperatura
        "\t"<<M_avg<<"\t"<<Mabs_avg<<"\t"<<Msq_avg<< //<M>;<|M|>;<M^2> por spin
        "\t"<<(Msq_avg-(M_avg*M_avg*n))/(tempFinal)<< //suceptibilidade por spin (X)
        "\t"<<(Msq_avg-(Mabs_avg*Mabs_avg*n))/(tempFinal)<< //susceptibilidade por spin (X’)
        "\t"<<E_avg<<"\t"<<Esq_avg<< //<E>;<E^2> por spin
        "\t"<<(Esq_avg-(E_avg*E_avg*n))/(tempFinal*tempFinal)<< // capacidade termica por spin
        "\t"<<1-((Mq_avg)/(3*Msq_avg))<<endl; //cumulante(U_L)
    }
    return 0;
}