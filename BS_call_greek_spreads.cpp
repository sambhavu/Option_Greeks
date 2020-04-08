#include <iostream>
#include <algorithm>
#include <cmath>

using namespace std;

#ifndef Pi 
#define Pi 3.141592653589793238462643 
#endif 


class normal{
public: 
double normdist(double x);
double cnd(double x);
};


double normal::normdist(double x){
	double num = 0;
	double w1 = 0;
	double w2 = 0;
	
	w1 = (1/sqrt(2*Pi));
	w2 = exp(-.5*x*x);
	
	num = w1*w2; 
	
	return num;
}

double normal::cnd(double x){
	
	  double L, K, w ;
	  double const a1 = 0.31938153, 
	  a2 = -0.356563782, a3 = 1.781477937;
  double const a4 = -1.821255978, a5 = 1.330274429;

	   L = fabs(x);
	   K = 1.0 / (1.0 + 0.2316419 * L);
	   w = 1.0 - 1.0 / sqrt(2 * Pi) * exp(-L *L / 2) * (a1 * K + a2 * K *K + a3 * pow(K,3) + a4 * pow(K,4) + a5 * pow(K,5));

  if (x < 0 ){
	w= 1.0 - w;
  }
  return w;
	
}



class option{
	public: 
	
	normal N; 
	
	double call;
	
	double delta;
	double gamma;
	double theta;
	double vega;
	double rho;
	
	double value(double k, double iv, double t, double p, double r);
	
	double d(double k, double iv, double t, double p, double r);
	
	   double g(double k, double iv, double t, double p, double r);

   double t(double k, double iv, double t, double p, double r);

   double v(double k, double iv, double t, double p, double r);

   double r(double k, double iv, double t, double p, double r);
	
	int type; 
	
	virtual void display(){ 
		cout<<"option";
	}
	
};

double option::value(double k, double iv, double t, double p, double r){

double d1; 
double d2;
double roott;
roott=sqrt(t);

d1=(log(p/k)+(r+(iv*iv/2)*iv*roott))/(iv*roott);

d2=d1-(iv*roott);

call=p*N.normdist(d1)-k*exp(-r*t)*N.cnd(d2);

return call;
}

double option::d(double k, double iv, double t, double p, double r){
   double d1; 
   double d;
   double roott;
   roott=sqrt(t);

d1=(log(p/k)+(r+(iv*iv/2)*iv*roott))/(iv*roott);

   d=N.cnd(d1);
	return d;
	
}

double option::g(double k, double iv, double t, double p, double r){
	 double d1; 
   double d;
   double roott;
   roott=sqrt(t);
   
d1=(log(p/k)+(r+(iv*iv/2)*iv*roott))/(iv*roott);

	d=N.normdist(d1)/(p*iv*roott);
	return d;
	
}

double option::t(double k, double iv, double t, double p, double r){
   double d1; 
   double d2;
   double d;
   double roott;
   roott=sqrt(t);

d1=(log(p/k)+(r+(iv*iv/2)*iv*roott))/(iv*roott);

d2=d1-iv*roott;

d=-(p*iv*N.normdist(d1))/(2*roott) - r*k*exp(-r*t)*N.cnd(d2);

return d;
	
}

double option::v(double k, double iv, double t, double p, double r){
	 double d1; 
   double d;
   double roott;
   roott=sqrt(t);

   d1=(log(p/k)+(r+(iv*iv/2)*iv*roott))/(iv*roott);

  d=p*N.normdist(d1)*roott;
  
  return d;
	
}

double option::r(double k, double iv, double t, double p, double r){
   double d1; 
   double d2;
   double d;
   double roott;
   roott=sqrt(t);

d1=(log(p/k)+(r+(iv*iv/2)*iv*roott))/(iv*roott);

d2=d1-iv*roott;


d=k*t*exp(-r*t)*N.cnd(d2);
return d;
	
}


class underlying: public option{
	public:
	
	double strike;
	double time;
	double impliedvol;
	double price;
	double intrate; 
		
};






class stock: public underlying {
	public:
	
	stock(double k, double t, double iv, double r, double p);

  void calculate();
  void display();
	

};


void stock::calculate(){

call=value(strike,impliedvol,time,price,intrate);

delta=d(strike,impliedvol,time,price,intrate);

gamma=g(strike,impliedvol,time,price,intrate);

vega=v(strike,impliedvol,time,price,intrate);
		
theta=t(strike,impliedvol,time,price,intrate);
		
rho=r(strike,impliedvol,time,price,intrate);
	
}

void stock::display(){ 


cout<<"\nStrike	\t\t\t\t:	"<<strike; 

cout<<"\nCall Price		   :\t"<<call;
cout<<"\n\nDelta			\t\t:\t"<<delta;
cout<<"\nGamma			\t\t:\t"<<gamma;
cout<<"\nVega			 \t\t:\t"<<vega;
cout<<"\nTheta			\t\t:\t"<<theta;
cout<<"\nRho			  \t\t:\t"<<rho;
	
cout<<"\n================================";
	
}

stock::stock(double k, double t, double iv, double r, double p)
{ 
	strike=k;
	time=t;
	impliedvol=iv;
	intrate=r;
	price=p; 
} 




class bond: public underlying {
	public:
	bond(double k, double t, 
	double iv, double r, double p); 
	
	
		void display() {
			cout << intrate;
		}
		
};

bond::bond(double k, double t, double iv, double r, double p){
	strike=k;
	time=t;
	impliedvol=iv;
	intrate=r;
	price=p;
}




int main() {

double price; 
double strike;
double t;
double intrate;
double vol; 

//double ytm;
//double div_yield; 

price=236.7;
t=5;
vol=.25;
strike =237; 
intrate=.03;

double step = 1;
stock apple(
	strike, t,
	vol,
	intrate,
   price); 

option *o1=&apple; 

cout<<"\n\nEquity Price\t\t\t:	"<<apple.price;
cout<<"\nDays\t\t\t\t\t:	"<<apple.time;
cout<<"\nInterest Rate\t\t\t:   "<<apple.intrate;
cout<<"\nImplied Vol\t\t\t: 	"<<apple.impliedvol;

for(int i=0;i<10;++i){ 
    

apple.calculate();

o1->display();

apple.strike+=step;

}

step=.01;

apple.strike=240;

for(int i=0;i<10;++i){ 
    

apple.calculate();

o1->display();

apple.impliedvol+=step;

}
	
	
	
return 0;
}
