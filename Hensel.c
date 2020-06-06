#include <stdio.h>  //this is it
#include <stdint.h>
#include <inttypes.h>
#include <stdlib.h>
#include <gmp.h>


//ce programme contient toutes les fonctions nécessaires pour déterminer les racines rationnelles d'un pokynome grace à la remontée d'Hensel
//le programme utilise la bibliothèque gmp donc on n'est pas limiter dans la taille des entiers.
	
void ext_eucl_div2(mpz_t *u, mpz_t *v, mpz_t *g, mpz_t a, mpz_t b, mpz_t l ){
	//renvoie pour a,b 3entiers u v g tq ua+bv=g=pgcd(a,b)
	mpz_t u1,u2,u3,v1,v2,v3,q,t1,t2,t3,tmp1;
	mpz_t tmp; mpz_init(tmp);

	mpz_init(u1); // initialise a 0
	mpz_init(u2); 
	mpz_init(u3); 
	mpz_init_set_ui(v1,1);
	mpz_init(v2); 
	mpz_init_set(v3,a);
	mpz_init(q);
	mpz_init(t1);
	mpz_init_set_ui(t2,1);
	mpz_init_set(t3,b);
	mpz_init_set(tmp1, v3);
	mpz_sub(tmp1, v3, l);
	
			gmp_printf("-*-- V3 = %Zd -----", v3);  //PB v3 prend direct la valeur de b et on verifie pas que c'est ou pas inferieur a l 
	while((mpz_cmp_si(t3,0)) && (mpz_cmp_si(tmp1,0)<0)  ){ //v3<l
		mpz_set(u1,v1);

		mpz_set(u2,v2);
		
		mpz_set(u3,v3);
		
		mpz_set(v1,t1);
		mpz_set(v2,t2);
		mpz_set(v3,t3);
		mpz_cdiv_q(q,u3,v3);
		

		mpz_mul(tmp, q, v1);
		mpz_sub(t1,u1,tmp);
		
		mpz_mul(tmp, q, v2);
		mpz_sub(t2,u2,tmp);

		mpz_mod(t3, u3,v3);
		
		
		mpz_sub(tmp1, v3, l);
		gmp_printf("---- V3 = %Zd -----", v3);

	}
	mpz_set(*u,v1);

	mpz_set(*v,v2);
	mpz_set(*g, v3);
	if(mpz_cmp_si(v3,0)<0){
		mpz_mul_si(*g,v3, -1 );

	}
}

void ext_eucl_div1(mpz_t *u, mpz_t *v, mpz_t *g, mpz_t a, mpz_t b){//OK
	//renvoie pour a,b 3entiers u v g tq ua+bv=g=pgcd(a,b)
	mpz_t u1,u2,u3,v1,v2,v3,q,t1,t2,t3; 
	
	mpz_t tmp; mpz_init(tmp);
	mpz_init(u1); // initialise a 0
	mpz_init(u2); 
	mpz_init(u3); 
	mpz_init_set_ui(v1,1);
	mpz_init(v2); 
	mpz_init_set(v3,a);
	mpz_init(q);
	mpz_init(t1);
	mpz_init_set_ui(t2,1);
	mpz_init_set(t3,b);
	
	
	while(mpz_cmp_ui(t3,0)){
	
	
		mpz_set(u1,v1);
		
		mpz_set(u2,v2);
		
		mpz_set(u3,v3);
		
		mpz_set(v1,t1);
		
		mpz_set(v2,t2);
		
		mpz_set(v3,t3);
		
		mpz_fdiv_q(q,u3,v3);
		
		mpz_mul(tmp, q, v1);
		
		mpz_sub(t1,u1,tmp);
		
		mpz_mul(tmp, q, v2);
		
		mpz_sub(t2,u2,tmp);
		
		mpz_mod(t3, u3,v3);
		
		

	}
	
	mpz_set(*u,v1);
	
	mpz_set(*v,v2);
	
	mpz_set(*g, v3);
	
	if(mpz_cmp_si(v3,0)<0){
		mpz_mul_si(*g,v3, -1 );
	}
}

void gcd(mpz_t a, mpz_t b, mpz_t *u, mpz_t *v, mpz_t *g){//OK
	//renvoie g=pgcd(a,b) en utilisant ext_eucl_div

	mpz_init(*g);
	ext_eucl_div1(u,v,g,a,b);
	gmp_printf("la relation de bezout est %Zd*%Zd + %Zd*%Zd = %Zd\n",a,*u,b,*v,*g);
}
int64_t *allocateVector(unsigned int m){//OK
	int64_t *v;
	v=calloc(m,sizeof(int64_t));
	if(v==NULL){
		printf("erreur allocation mémoire");
		exit(1);
	}
	return v;
}
void copyVector(int64_t *v,int64_t *u,unsigned int m){//OK
	unsigned int i;
	for(i=0;i<m;i++){
		v[i]=u[i];
	}

}

void printVector(int64_t *v,unsigned int m){//OK
	unsigned int i;
	for(i=0;i<m;i++){
		printf(" %d ",v[i]);
	}
	printf("\n");
}

void randomVector(int64_t *v, unsigned int m){//OK
	unsigned int i;
	for(i=0;i<m;i++){
		v[i]=(int64_t)rand();
	}
	
}
void readVector(int64_t *v, unsigned int m){//OK
	unsigned int i;
	for(i=0;i<m;i++){
		printf("saisir v[%d]\n",i);
		if(!(scanf("%d",&v[i]))){
			exit(1);
		}
	}
}



void evaluatePolynomial(int64_t *c,mpz_t x,mpz_t mod, unsigned int n, mpz_t* q){//OK
	int i;
	for(i=n-1;i>0;i--){
		mpz_add_ui(*q,*q,c[i]); // q= q+c[i] 
		mpz_mul(*q,*q,x);// q= q*x
		
		
	}
	
	mpz_add_ui(*q,*q,c[0]);
	

	mpz_mod(*q,*q,mod);
	gmp_printf("evaluation polynomiale donne %Zd\n",*q);
	
	 //return mpz_mod(p,p,mod);
}

void derivePolynome(int64_t *c,unsigned int n){//OK
	int i;
	for(i=0;i<n;i++){
		c[i]=c[i+1]*(i+1);
	}
	c[n]=0;
}


void racineMod(mpz_t* x0, int64_t *c, mpz_t *p, int n,unsigned int taillepoly){   // Newton , x doit etre la racine (search x)
	mpz_t  x0tmp;
	int64_t *cprim=allocateVector(taillepoly); // rentrer en parametre
	mpz_init_set(x0tmp,*x0); 
	
	mpz_t xprim,u,v,g;
	mpz_init(u);
	mpz_init(v);
	mpz_t mod;
	mpz_t q;
	mpz_init(q);
	mpz_t m;
	mpz_init(m);
	printf("allocation reussi \n");
	mpz_init(xprim);
	mpz_init_set(mod, *p);
	mpz_mul(mod,mod,mod);
	gmp_printf("mod=%Zd, x0tmp=%Zd\n",mod,x0tmp);
	int i;
	copyVector(cprim,c,taillepoly);
	derivePolynome(cprim, taillepoly);
	printf("polynome derivé:\t");
	printVector(cprim,taillepoly);
	
	for(i=0;i<n;i++){
		printf("rentre dans le for\n");
		gmp_printf("q=%Zd, mod=%Zd, taillepoly=%d\n", q,mod,taillepoly);
		evaluatePolynomial(cprim,x0tmp,mod,taillepoly,&q);//"renvoi" q//SA MARCHE
		printf("evaluation effectué\n");
		mpz_set(xprim,q);// on met q ds xprim ->> x' = c'(x)
		mpz_set_ui(q,0);
		
		gmp_printf("xprim=%Zd, q=%Zd\n",xprim, q);
		if((mpz_cmp_ui(xprim,0)==0)){  // c'[x] ==0
			printf("Racine double trouvé \n");
			
			break;
		}
		
		gcd(xprim,mod,&u,&v,&g);   // g = pgcd(x', mod) on est interressé par l'inverse de x' modulo mod  (x' = c'[x] != 0)  1/c' vrai que si g=1!!!
//x*u+mod*v=g => x*u=1[mod]
		
		mpz_set(xprim,u);  // x' = u = inverse de x'
		gmp_printf("xprim=%Zd\n",xprim);
		evaluatePolynomial(c,x0tmp,mod,taillepoly,&m);  
		mpz_mul(xprim,m,xprim);  //u= c(x)/c'(x)
		mpz_set_ui(m,0);
		gmp_printf("xprim=%Zd\n",xprim);
		mpz_sub(x0tmp,x0tmp,xprim);
		mpz_mod(x0tmp,x0tmp,mod);
		gmp_printf("x0tmp=%Zd\n",x0tmp);
		mpz_mul(mod,mod,mod); //res = mod carre
		gmp_printf("mod=%Zd\n",mod);
		printf("%deme tour de boucle effectué\n",i);
	}
	mpz_set(*p,mod);
	
	mpz_set(*x0,x0tmp);
	//return x;
}

void searchX0(int64_t *c, mpz_t p, unsigned int n, mpz_t* x0){//OK
	
	mpz_t q, x0tmp;
	mpz_init_set_ui(q,0);
	
	mpz_init_set(x0tmp,*x0);
	evaluatePolynomial(c,x0tmp,p,n,&q); //q = polynome(x0=0)
	
	
	while((mpz_cmp_ui(q,0))&&(mpz_cmp(x0tmp,p)<0)){ // on regarde les entier inferieur au modulo p et des qu'on trouve q=0 cad on trouve la racine
		 mpz_add_ui(x0tmp,x0tmp,1); // on incremente jusqu'a trouver la racine
		mpz_init_set_ui(q,0);
		evaluatePolynomial(c,x0tmp,p,n,&q); // on teste avec aec le nouvel abscisse x0
	
	}
	

	mpz_set(*x0,x0tmp);
	if(mpz_cmp(x0tmp,p)==0){
		printf("pas de zero trouvé\n");
		//return  0;  //return mpz_t-1;

	}
	//return x0;
}


/*int main(){
	mpz_t a,b,u,v,g,l;
	mpz_init_set_ui(l,150);
	mpz_t tab[2];
	mpz_init(tab[0]);
	mpz_init(tab[1]);
	mpz_set_si(tab[0],1);
	//gmp_scanf("%Zd", tab[0]);
	gmp_printf("*******************%Zd", tab[0]);
	//gmp_scanf("%Zd", tab[1]);
	//gmp_printf("*******************%Zd", tab[1]);
	mpz_init_set_si(a, 155);
	mpz_init_set_si(b, 5);


	mpz_init(u);mpz_init(v);mpz_init(g);
	
	ext_eucl_div1( &u, &v, &g,a,b	);
	gmp_printf("////// %Zd /////// \n", g);
return 0;
}*/

int main(){
	mpz_t x0,p;
	int64_t* vect=allocateVector(4);
	int64_t* vect2=allocateVector(4);
	readVector(vect, 4);
	
	mpz_init(x0);
	mpz_init_set_ui(p,7);
	searchX0(vect,p ,4,&x0);
	gmp_printf("%Zd est une racine du polynome modulo %Zd \n",x0,p);
	
	racineMod(&x0, vect, &p,5, 4);
	gmp_printf("%Zd est une racine du polynome modulo %Zd \n",x0,p);
return 0;}

/*
void ext_eucl_div2(mpz_t *u, mpz_t *v, mpz_t *g, mpz_t a, mpz_t b, mpz_t l ){
	//renvoie pour a,b 3entiers u v g tq ua+bv=g=pgcd(a,b)
	mpz_t u1,u2,u3,v1,v2,v3,q,t1,t2,t3,tmp1;
	mpz_t tmp; mpz_init(tmp);

	mpz_init(u1); // initialise a 0
	mpz_init(u2); 
	mpz_init(u3); 
	mpz_init_set_ui(v1,1);
	mpz_init(v2); 
	mpz_init_set(v3,a);
	mpz_init(q);
	mpz_init(t1);
	mpz_init_set_ui(t2,1);
	mpz_init_set(t3,b);
	mpz_init_set(tmp1, v3);
	mpz_sub(tmp1, v3, l);
	while(t3!=0 && (mpz_cmp_si(tmp1,0)<0)  ){ //v3<l
		mpz_set(u1,v1);

		mpz_set(u2,v2);
		
		mpz_set(u3,v3);
		
		mpz_set(v1,t1);
		mpz_set(v2,t2);
		mpz_set(v3,t3);
		mpz_cdiv_q(q,u3,v3);
		

		mpz_mul(tmp, q, v1);
		mpz_sub(t1,u1,tmp);
		
		mpz_mul(tmp, q, v2);
		mpz_sub(t2,u2,tmp);

		mpz_mod(t3, u3,v3);
		
		
		mpz_sub(tmp1, v3, l);

	}
	mpz_set(*u,v1);

	mpz_set(*v,v2);
	mpz_set(*g, v3);
	if(mpz_cmp_si(v3,0)<0){
		mpz_mul_si(*g,v3, -1 );

	}
}

void gcd(mpz_t a, mpz_t b, mpz_t *u, mpz_t *v, mpz_t *g){
	//renvoie g=pgcd(a,b) en utilisant ext_eucl_div

	mpz_init(*g);
	ext_eucl_div1(u,v,g,a,b);
	gmp_printf("pgcd de %Zd et %Zd est %Zd\n",a,b,*g);
}

void allocateVector(mpz_t *v,unsigned int m){
	mpz_realloc2(*v, m*sizeof(mpz_t));
	if(v==NULL){
		printf("erreur allocation mémoire");
		exit(1);
	}
}




void copyVector(mpz_t *v,mpz_t *u,unsigned int m){
	unsigned int i;
	for(i=0;i<m;i++){
		mpz_set(v[i],u[i]);
	}

}  // erreur de segmentation avec mpz_set (pointeur) !!!!!!!!!!


void printVector(mpz_t *v,unsigned int m){
	unsigned int i;
	for(i=0;i<m;i++){
		gmp_printf(" %Zd ",v[i]);
	}
	gmp_printf("\n");
}

void randomVector(mpz_t *v, unsigned int m){
	unsigned int i;
	gmp_randstate_t state ;
	gmp_randinit_default (state);
	for(i=0;i<m;i++){
		mpz_urandomb (v[i],state,64) ;// renvoi nb entre 0 et 2^64-1
	}
	
}




void readVector(mpz_t *v, unsigned int m){
	unsigned int i;
	for(i=0;i<m;i++){
		gmp_printf("saisir v[%Zd]\n",i);
		if(!(gmp_vscanf("%Zd",&v[i]))){
			exit(1);
		}
	}
}






void evaluatePolynomial(mpz_t *c,mpz_t x,mpz_t mod, unsigned int n, mpz_t q){
	int i;
	for(i=n-1;i>=0;i--){
		mpz_add(q,q,c[i]); // q= q+c[i] 
		mpz_mul(q,q,x);// q= q*x
		
	}
	
	mpz_mod(q,q,mod);
	 //return mpz_mod(p,p,mod);
}

void derivePolynome(mpz_t *c,unsigned int n){

	int i;
	for(i=0;i<n;i++){
		mpz_t k;
		mpz_set_ui(k,(i+1));	
		mpz_mul(c[i],c[i+1],k);

	}
	mpz_set_ui(c[0],0);
}

void racineMod(mpz_t x, mpz_t *c, mpz_t *p, int n,int deg){   // Newton , x doit etre la racine (search x)
	mpz_t *cprim; allocateVector(cprim,deg); // rentrer en parametre 
	mpz_t xprim,u,v,g;
	mpz_t mod;
	mpz_t q;
	mpz_init(q);
	mpz_t m;
	mpz_init(m);
	mpz_t res;
	mpz_init(res);

	mpz_t zero;
	mpz_init(zero);
	
	mpz_t test;
	mpz_init(test);
	mpz_set(mod, *p) ;//mpz_init(mod, *p);
	int i;
	copyVector(cprim,c,deg);
	derivePolynome(cprim, deg);
	for(i=0;i<n;i++){
		evaluatePolynomial(cprim,x,mod,deg,q) ;//"renvoi" q
		mpz_set(xprim,q);// on met q ds xprim ->> x' = c'(x)
		mpz_mod(test,xprim,mod); //test = x'% mod **************************************************************************************** peut etre le mod a deja ete fait  
		
		if(mpz_cmp_ui(test,0)==0){  // c'[x] ==0
			gmp_printf("Racine double trouvé \n");
			//exit(1);
			//continue;
			break;
		}
		
		gcd(xprim,mod,&u,&v,&g);   // g = pgcd(x', mod) on est interressé par l'inverse de x' modulo mod  (x' = c'[x] != 0)  1/c'
		mpz_set(xprim,u);  // x' = u = inverse de x'
		evaluatePolynomial(c,x,mod,deg,m);  
		mpz_mul(u,m,u);  //u= c(x)/c'(x)
		mpz_sub(x,x,u);
		mpz_mul(mod,mod,mod); //res = mod carre
	}
	mpz_set(*p,mod);
	//return x;
}

void searchX0(mpz_t *c, mpz_t p, unsigned int n, mpz_t x0){

	mpz_t q;
	mpz_init(q);
	mpz_t res;
	mpz_init(res);
	mpz_t p_x0;
	evaluatePolynomial(c,x0,p,n,q);
	mpz_set(p_x0,q); 
	while((p_x0!=0)&&x0<p){
		 mpz_add_ui(res,x0,1);
		evaluatePolynomial(c,x0,p,n,q);
		 mpz_set(p_x0,q);
	}
	if(x0==p){
		printf("pas de zero trouvé\n");
		//return  0;  //return mpz_t-1;
	}
	//return x0;
}




/*
int main(){
mpz_t u,v,a,b,g;
mpz_init_set_ui(u,2);
mpz_init_set_ui(v,5);
mpz_init(a);
mpz_init(b);
mpz_init(g);
ext_eucl_div1(u,v, g, &a, &b);
gmp_printf(" pgcd est %Zd et u %Zd et v: %Zd ", g,u,v);


return 0;
}
		
int main() {
	int i;
	mpz_t q;
	mpz_init(q);
	mpz_t *c;
	//mpz_set(c,allocateVector(4));   //(c,allocateVector(4));
	mpz_t p;
	for(i=0;i<20;i++){
		randomVector(c,4);
		printVector(c,4);
		mpz_init_set_ui(p,5);
		mpz_t x0;
		searchX0(c,p,3,x0);
		if(mpz_cmp_ui(x0,-1)<0 || mpz_cmp_ui(x0,-1)>0 ){
			racineMod(x0,c,&p,10,3);
			mpz_set(x0,p);
			evaluatePolynomial(c,x0,p,3,q)	;		
			if(q==0)
				gmp_printf("Test reussi");
			
			else
				gmp_printf("Test échoué");
		}
	}
	free(c);NIsansMpz.c: In function ‘racineMod’:
NIsansMpz.c:217:38: warning: passing argument 5 of ‘evaluatePolynomial’ from incompatible pointer type
   evaluatePolynomial(cprim,x,mod,deg,q) ;//"renvoi" q
                                      ^
NIsansMpz.c:175:6: note: expected ‘struct __mpz_struct (*)[1]’ but argument is of type ‘struct __mpz_struct *’
 void evaluatePolynomial(int64_t *c,mpz_t x,mpz_t mod, unsigned int n, mpz_t* q){
      ^
NIsansMpz.c:230:34: warning: passing argument 5 of ‘evaluatePolynomial’ from incompatible pointer type
   evaluatePolynomial(c,x,mod,deg,m);  
                                  ^
NIsansMpz.c:175:6: note: expected ‘struct __mpz_struct (*)[1]’ but argument is of type ‘struct __mpz_struct *’
 void evaluatePolynomial(int64_t *c,mpz_t x,mpz_t mod, unsigned int n, mpz_t* q){
      ^
3670252@ppti-1
	return 0;
}


////////////////////////euclide
int main(){
	printf("test -1");
	mpz_t a,b;
	mpz_t u,v,g;
	printf("test0");
	mpz_set_ui(a,12);
	mpz_set_ui(b,3);
	mpz_init(u); 
	mpz_init(v);
	mpz_init(g); 
	printf("test1");
	ext_eucl_div(u, v, g, a, b);
	gmp_printf("pgcd de %Zd et %Zd est %Zd et ses coefs de bezout sont %Zd et %Zd \n",a,b,g, u, v);
	return 0;
}



int main (){
	//int *i;
	mpz_t v[2] ;
	mpz_init(v[0]); // Comment affecter une valeur autre que 0 à v[i]???????? mpz_set ne marche pas provoque un seg fault
	mpz_init(v[1]);
	gmp_printf("%Zd",v[0]);
	mpz_t u[2];
	copyVector(u,v,2);
	for(i=0;i<2;i++){
		gmp_printf(" %Zd ",u[i]);
	}
	printf("test\n");
	mpz_t v[2];
	mpz_inits(v[0],v[1]);
	randomVector(v,2);
	//readVector(v,2);
	printVector(v,2);
	//gmp_printf("i\n");
	mpz_t x;
	mpz_t c[3];
	mpz_inits(c[0],c[1],c[2]);
	mpz_t q;
	mpz_t mod;
	mpz_inits(mod,x);
	printf("test");
	evaluatePolynomial(c,x,mod,3,q);
	gmp_printf("%Zd",q);
	

return 0;
}
*/
