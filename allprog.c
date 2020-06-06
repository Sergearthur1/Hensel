
#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdlib.h>
#define LIN(i,j,m,n) (i*n+j)
#define TAILLEM 3

typedef struct _listeMod{
	int64_t* val;
	int64_t mod;
	struct _listeMod *suiv;
}listeMod;


void ext_eucl_div(int64_t *u, int64_t *v, int64_t *g, int64_t *a, int64_t *b){
	//renvoie pour a,b 3entiers u v g tq ua+bv=g=pgcd(a,b)
	int64_t u1=0,u2=0,u3=0,v1=1,v2=0,v3=*a,q=0,t1=0,t2=1,t3=*b; //algo du prof
	while(t3!=0){
		u1 = v1;
		//printf("%d\n",u1);
		
		u2 = v2;
		//printf("%d\n",u2);
		
		u3 = v3;
		//printf("%d\n",u3);
		
		v1 = t1;
		//printf("%d\n",v1);
		v2 = t2;
		//printf("%d\n",v2);
		v3 = t3;
		//printf("%d\n",v3);
		q = u3/v3;
		//printf("%d\n",q);
		t1 = u1-q*v1;
		//printf("%d\n",t1);
		t2 = u2-q*v2;
		//printf("%d\n",t2);
		t3 = u3%v3;
		//printf("%d\n\n\n",t3);
	}
	*u = v1;
	//printf("%d\n",*u);
	*v = v2;//printf("%d",*v);
	*g = v3; //marche lorsque b=0 : g=a
	//printf("%d",*g);
	if(v3<0){
		*g=(-1)*v3;
	}
}

int64_t gcd(int64_t a, int64_t b, int64_t *u, int64_t *v){
	//renvoie g=pgcd(a,b) en utilisant ext_eucl_div
	int64_t g=0;
	ext_eucl_div(u,v,&g,&a,&b);
	printf("pgcd de %ld et %ld est %ld\n",a,b, g);
	return g;
}

void thchinois(int64_t a1 ,int64_t a2 , int64_t m1, int64_t m2, int64_t *mod, int64_t *x){
    int64_t   u,v,pgcd ; 
   // printf("0\n");
    
  /*  *x = a1;
    printf("1\n");
    
    pgcd = gcd(M, m2,&u,&v);		 
    printf("%d   %d   %d  \n",u,v,pgcd);
    
    *x = (u)*m2*(*x) + (v)*M*a2; 	 
    printf("%d\n",*x);
    
    M= m2*M;
    printf("%d\n",M);
    
    *x = (*x)%M;
    printf("%d\n",*x);
    
    *mod = m1*m2; 
    printf("%d\n",*mod);
*/



	//methode arthur
	
	
	pgcd = gcd(m1, m2,&u,&v);
	u=u*(a1-a2);
	v=v*(a1-a2);
	//*x=u-a1;
	*x=m2*v+a2;
	printf("%ld %ld\n",*x , -m1*u+a1);
	
	*mod = m1*m2; 
	(*x)=(*x)%(*mod);
	if((*x)<0){
		(*x)=*mod+*x;
	}
	printf("%ld\n",*mod);
	
	
	
       
}









int64_t *allocateVector(unsigned int m){
	int64_t *v;
	v=calloc(m,sizeof(int64_t));
	if(v==NULL){
		printf("erreur allocation mémoire");
		exit(1);
	}
	return v;
}
int64_t *allocateMatrix(unsigned int m,unsigned int n){
	int64_t *M;
	M=calloc(m*n,sizeof(int64_t));
	if(M==NULL){
		printf("erreur allocation mémoire");
		exit(1);
	}
	return M;
}


void readVector(int64_t *v, unsigned int m){
	unsigned int i;
	for(i=0;i<m;i++){
		printf("saisir v[%d]\n",i);
		if(!(scanf("%d",&v[i]))){
			exit(1);
		}
	}
}

void readMatrix(int64_t *A, unsigned int m, unsigned int n){
	unsigned int i,j;
	for(i=0;i<m;i++){
		for(j=0;j<n;j++){
			printf("saisir [%d][%d]\n", i,j);
			if(!(scanf("%d",&A[LIN(i,j,m,n)]))){
				exit(1);
			}
		}
	}printf("fin de la saisi de la matrice\n\n");
}



void printVector(int64_t *v,unsigned int m){
	unsigned int i;
	for(i=0;i<m;i++){
		printf(" %d ",v[i]);
	}
	printf("\n");
}
void printMatrix(int64_t *A, unsigned int m, unsigned int n){
	unsigned int i,j;
	for(i=0;i<m;i++){
		printf("\n");
		for (j=0;j<n;j++){
			printf("%d\t",A[LIN(i,j,m,n)]);
		}
	printf("\n\n");
	}
}

void freeVector(int64_t *v){
	free(v);
}
void freeMatrix(int64_t *A){
	free(A);
}


void copyVector(int64_t *v,int64_t *u,unsigned int m){
	unsigned int i;
	for(i=0;i<m;i++){
		v[i]=u[i];
	}

}

void copyMatrix(int64_t *B,int64_t *A, unsigned int m, unsigned int n){
	unsigned int i,j;
	for(i=0;i<m;i++){
		for (j=0;j<n;j++){
			B[LIN(i,j,m,n)]=A[LIN(i,j,m,n)];
		}
	}
}




	

void gaussianElimination(int64_t *A,int64_t *b, unsigned int n, unsigned int mod){
	unsigned int i,j,k;
	int64_t coeff,u,v; // coeff l'element collonne du pivot

	for(i=0;i<n;i++){ // on parcours les pivots.

		if(gcd(A[i*n+i] , mod, &u, &v)!=1){ //on verifie qu'il existe un inverse modulaire pour le pivot avec le mod. 
				printf("Calcul de l'inverse Impossible\n");
				u=0;
			}		

		for(j=1;j<n-i;j++){ // on parcours toutes les lignes en dessous du pivot 

			

			coeff=(A[(i+j)*n+i])*u; //on multiplie le coeff par l'inverse modulaire
			coeff=coeff%mod; // on ajuste par rapport au modulo
			
			for(k=i;k<n;k++){ // on parcours les colonnes de la ligne en question, en partant de l'indice de la colonne du pivot 
				//on modifie la matrice A						
				A[(i+j)*n+k]= A[(i+j)*n+k] - coeff*A[i*n+k];
				A[(i+j)*n+k]= A[(i+j)*n+k] % mod;
				if(A[(i+j)*n+k]<0){//on normalise
					A[(i+j)*n+k]= A[(i+j)*n+k] +mod;
				}			
			}
				//On fait pareil avec b
			b[i+j]= b[i+j] - coeff*b[i];
			b[i+j]= b[i+j] % mod;
			if(b[i+j]<0){//on normalise
				b[i+j]= b[i+j] + mod;
			}
		}
	}
}

void solveTriangularUpper(int64_t *x,int64_t *A,int64_t *b,unsigned int n, unsigned int mod){
	/* On resoud le systeme a partir de la matrice triangulaire sup*/

	int i,j;
	int64_t u,v; 
	
	for(i=n-1;i>=0;i--){ //Parcours des lignes en partant du bas
		x[i]=b[i]%mod; 

		for(j=n-1;j>i;j--){ // Parcours des colonnes 
			x[i]=x[i]-(A[i*n+j]*x[j]);
			x[i]= x[i] % mod;
		}
		if(gcd(A[i*n+i] , mod, &u, &v)!=1){
				printf("Calcul de l'inverse Impossible\n");
				u=0;
			}
		x[i]=x[i] * u; 
		x[i]=x[i] % mod;
		if(x[i]<0){//on normalise
			x[i]= x[i] + mod;
		}
	}
}

void solveSystemGauss(int64_t *x,int64_t *A,int64_t *b,unsigned int n, unsigned int mod){

	int64_t *tmpMatrix=allocateMatrix(n,n);
	int64_t *tmpVector=allocateVector(n);

	copyMatrix(tmpMatrix,A,n,n);
	copyVector(tmpVector,b,n);
	gaussianElimination(tmpMatrix,tmpVector,n, mod);

	solveTriangularUpper(x,tmpMatrix,tmpVector,n, mod);
	freeVector(tmpVector);
	freeMatrix(tmpMatrix);
}

listeMod* creerElem(int64_t mod){
	listeMod *nv=malloc(sizeof(listeMod));
	nv->val=allocateVector(TAILLEM);
	nv->mod=mod;
	nv->suiv=NULL;
	return nv;
}

listeMod* addElem(listeMod *l,int64_t mod){
	listeMod *nv=creerElem(mod);
	nv->suiv=l;
	return nv;
}

void freeElem(listeMod *l){
	while(l){
		listeMod *tmp=l;
		l=l->suiv;
		free(tmp);
	}
	
}







int calculMod(listeMod *lm){
	listeMod *l=lm;
	int res=l->mod;
	l=l->suiv;
	while(l){
	res=res*l->mod;
	l=l->suiv;
	}
}


void thchinoisAmeliore(listeMod *lm,int i){
	listeMod *l=lm;
	int64_t modtmp=l->mod;
	printf("%d\n",modtmp);
	int64_t valtmp;
	while(l->suiv!=NULL){
		thchinois(l->val[i], l->suiv->val[i], modtmp, l->suiv->mod, &modtmp, &valtmp);
		printf("%d\n",modtmp);
		l=l->suiv;
		l->val[i]=valtmp;
		
		
	}
	printf("le resultat est %ld mod(%ld)\n",l->val[i],modtmp);
}
void thchinoisvect(listeMod *lm){
	int i;
	for(i=0;i<TAILLEM;i++){
		thchinoisAmeliore(lm,i);
	}
	while(lm->suiv){
		lm=lm->suiv;
	}
	printf("le resultat final est :     ");
	printVector(lm->val,TAILLEM);
		

}


int main(){
	int taille=TAILLEM;
	int i,nbsystem,mod;
	printf("resolution de systeme modulaire de taille %d, pour chnger la taille modifier la macro TAILLEM\n",taille);
	int64_t *A=allocateMatrix(TAILLEM, TAILLEM);
	int64_t *B=allocateVector(TAILLEM);
	printf("rentrer la matrice associer au systeme\n");
	readMatrix(A,TAILLEM,TAILLEM);
	printf("rentrer le vecteur resultat\n");
	readVector(B,TAILLEM);
	printf("rentrer le nombre de systeme a resoudre\n");
	scanf("%d",&nbsystem);
	printf("rentrer le modulo du 1er systeme\n");
	scanf("%d",&mod);
	listeMod *lm=creerElem(mod);
	solveSystemGauss(lm->val,A,B,TAILLEM, mod);
	printVector(lm->val,TAILLEM);
	for(i=1;i<nbsystem;i++){
		printf("rentrer le modulo du %deme systeme\n",i+1);
		scanf("%d",&mod);
		lm=addElem(lm,mod);
		solveSystemGauss(lm->val,A,B,TAILLEM, mod);
		printVector(lm->val,TAILLEM);
	}
	thchinoisvect(lm);
	
	mod=calculMod(lm);
	
	

	freeElem(lm);
	
		

	return 0;
}

	
