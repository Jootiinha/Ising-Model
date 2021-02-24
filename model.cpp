#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#define L 8
#include <gsl/gsl_rng.h> /*Compilar com: g++ model.cpp -lgsl -lgslcblas -lm*/
using namespace std;
ofstream DATA("DATA.1.dat",ios::out);
//Estrutura para a Matriz 2D
struct matrix{
    int x;
    int y;
};
const int n=L*L; //# de spins
double T=5; //Temperatura Final
const float Tmin=0.5; //temperatura inicial
const float dT=0.1; //passo para a temperatura no loop de temperatura
int S[L][L]; //rede bidimensional
long unsigned int MS=10000; //numero de steps de MC
int transiente=1000; //numero de passos descartados (para evitar a parte fora do equilíbrio.
//int seed=1243; //semente para o gerador
const gsl_rng *r;
//Inicializar a rede
void inicia_rede(){
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
void escolha_spin(matrix &pos){
    pos.x=(gsl_rng_uniform(r)*(L));// (0,L-1)
    pos.y=(gsl_rng_uniform(r)*(L));
}
//Calculo da energia associada a posição
int energy_pos(matrix &pos)
{
    //condições periódicas de contorno
    int up, down, left, rigth,e;
    rigth=(pos.x+1)%L;
    left=(pos.x+L-1)%L;
    up=(pos.y+1)%L;
    down=(pos.y+L-1)%L;
    //energia da posição
    e=-1*S[pos.x][pos.y]*(S[left][pos.y]+S[rigth][pos.y]+S[pos.x][up]+S[pos.x][down]);//J=1
    //printf("S[x=%d][y=%d], S[L=%d][y=%d], S[R=%d][y=%d], S[x=%d][U=%d], S[x=%d][D=%d] \n", pos.x,pos.y,
    return e;
}
//função decisão
bool teste_flip(matrix pos, int &de)
{
    de=-2*energy_pos(pos);
    if(de<0)
        return true; //flip pela energia
    else{
        if(gsl_rng_uniform(r)<exp(-de/T))
            return true; //flip pelo banho térmico
        else
            return false; //não flip
    }
}
//função de flip
void flip(matrix pos){
    S[pos.x][pos.y]=-S[pos.x][pos.y];
}
//função que remove o transiente
void transiente_final(){
    matrix pos;
    int de=0;
    for(int a=1;a<=transiente;a++){
        for(int b=1;b<=n;b++){
            escolha_spin(pos);
            if(teste_flip(pos,de)){
                flip(pos);
            }
        }
    }
}
//Magnetização total da rede (soma dos spins)
int magnetization(){
int m=0;
    for(int y=L-1;y>=0;y--){
        for(int x=0;x<L;x++){
            m+=S[x][y];
        }
    }
    return m;
}
//Energia total da rede
int energia_total()
{
    matrix pos;
    int e=0;
    for (int y=L-1;y>=0;y--){
        pos.y=y;
        for(int x=0;x<L;x++){
            pos.x=x;
            e+=energy_pos(pos);
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
    matrix pos;
    //inicia a rede de maneira aleatória
    inicia_rede();
    //loop de temperatura
    for(;T>=Tmin;T-=dT){
        //transiente
        transiente_final();
        //observáveis com valores no equilíbrio (após o transiente)
        M=magnetization();
        // Mabs=fabs(magnetization());
        E=energia_total();
        //reseta as variáveis para cada loop
        etotal=0;
        etotalsq=0;
        mtotal=0;
        mtotalsq=0;
        mabstotal=0;
        mqtotal=0;
        //loop do Monte Carlo (MS=tempo)
        for(int a=1;a<=MS;a++){
            //Metropolis
            for(int b=1;b<=n;b++){
                escolha_spin(pos);
                if(teste_flip(pos,de)){
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
        E_avg=etotal/(MS*n);
        Esq_avg=etotalsq/(MS*n);
        M_avg=mtotal/(MS*n);
        Mabs_avg=mabstotal/(MS*n);
        Msq_avg=mtotalsq/(MS*n);
        Mq_avg=mqtotal/(MS*n);
        //printf("%lf %Lf %Lf \n",T,Mabs_avg,M_avg);
        //saída para arquivo
        DATA<<T<< //temperatura
        "\t"<<M_avg<<"\t"<<Mabs_avg<<"\t"<<Msq_avg<< //<M>;<|M|>;<M^2> por spin
        "\t"<<(Msq_avg-(M_avg*M_avg*n))/(T)<< //suceptibilidade por spin (X)
        "\t"<<(Msq_avg-(Mabs_avg*Mabs_avg*n))/(T)<< //susceptibilidade por spin (X’)
        "\t"<<E_avg<<"\t"<<Esq_avg<< //<E>;<E^2> por spin
        "\t"<<(Esq_avg-(E_avg*E_avg*n))/(T*T)<< // capacidade termica por spin
        "\t"<<1-((Mq_avg)/(3*Msq_avg))<<endl; //cumulante(U_L)
    }
    return 0;
}