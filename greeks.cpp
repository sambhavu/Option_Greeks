#include <vector>
#include <iostream>
#include <string>
#include<cmath>
using namespace std;

#ifndef Pi 
#define Pi 3.141592653589793238462643 
#endif 

class option{ 
    public: 
    double S; 
    double K;
    double r; 
    double r_f;
    double sigma; 
    double time; 
    int steps;  
    
    option(); 
    
    double cnd(double x);
    double normdist(double x);
    
    double price(        
                            double S,
                            double K,
                            double r,
                            double sigma,
                            double time); 
                            
    void greeks(            double S,
                            double K,
                            double r,
                            double sigma,
                            double time,
                            double delta,
                            double gamma,
                            double theta,
                            double vega, 
                            double rho);
                            
    double currency_option(
                            double S,
                            double K,
                            double r,
                            double r_f,
                            double sigma,
                            double time, 
                            int steps);
};

option :: option(){
 S =    50;
 K = 50;
 r = .1;
 r_f = 0.05;
 sigma = .3;
 time = .5;
 
 steps =100;
 
}

double option:: normdist(double x)
{
    double num = 0;
    double w1 = 0;
    double w2 = 0;
    
    w1 = (1/sqrt(2*Pi));
    w2 = exp(-.5*x*x);
    
    num = w1*w2; 
    
    return num;
}

    
double option:: cnd(double x)
{
  double L, K, w ;
  /* constants */
  double const a1 = 0.31938153, a2 = -0.356563782, a3 = 1.781477937;
  double const a4 = -1.821255978, a5 = 1.330274429;

  L = fabs(x);
  K = 1.0 / (1.0 + 0.2316419 * L);
  w = 1.0 - 1.0 / sqrt(2 * Pi) * exp(-L *L / 2) * (a1 * K + a2 * K *K + a3 * pow(K,3) + a4 * pow(K,4) + a5 * pow(K,5));

  if (x < 0 ){
    w= 1.0 - w;
  }
  return w;
}



double option :: price(
                            double S,
                            double K,
                            double r,
                            double sigma,
                            double time)
{
        
        double time_sqrt = sqrt(time);
        double d1 = (log(S/K)+r*time)/(sigma*time_sqrt)+.5*sigma*time_sqrt;
        double d2 = d1-(sigma*time_sqrt);
        
        double price = S*cnd(d1)-K*exp(-r*time)*cnd(d2);
        
        return price; 
    }
    
    
void option :: greeks(
                        double S,
                        double K,
                        double r, 
                        double sigma,
                        double time,
                        double delta, 
                        double gamma, 
                        double theta,
                        double vega, 
                        double rho) 
                        
    {
        double time_sqrt = sqrt(time);
        double d1 = ((log(S/K)+r*time)/(sigma*time_sqrt)) + 0.5*sigma*time_sqrt;
        double d2 = d1-(sigma*time_sqrt);
        delta = cnd(d1);
        gamma = normdist(d1)/(S*sigma*time_sqrt);
        theta =- (S*sigma*normdist(d1))/(2*time_sqrt) - r*K*exp(-r*time)*cnd(d2);
        vega = S*time_sqrt*normdist(d1);
        rho = K*time*exp(-r*time)*cnd(d2);
        
        
        cout<<delta<<"\n";
        cout<<gamma<<"\n";
        cout<<theta<<"\n";
        cout<<vega<<"\n";
        cout<<rho<<"\n";
    }
    
double option :: currency_option( 
                                double S, 
                                double K,
                                double r,
                                double r_f,
                                double sigma,
                                double time,
                                int steps)
    {
        vector<double> exchange_rates(steps+1);
        vector<double> call_values(steps+1);
        
        double t_delta = time/steps;
        double Rinv = exp(-r*(t_delta));
        double u = exp(sigma*sqrt(t_delta));
        double d = 1.0/u;
        double uu = u*u;
        double pUp = (exp((r-r_f)*t_delta)-d)/(u-d);
        double pDown = 1.0 - pUp;
        
        exchange_rates[0] = S*pow(d,steps);
        int i;
        for(i = 1;i<=steps;++i){
            exchange_rates[i] = uu*exchange_rates[i-1];         //terminal tree nodes
        };
        for(i=0;i<=steps;++i) call_values[i]=max(0.0,(exchange_rates[i]-K));
        for(int step = steps-1;step>=0;--step){
            exchange_rates[i] = d*exchange_rates[i+1];
            call_values[i] = (pDown*call_values[i] + pUp*call_values[i+1])*Rinv;
            call_values[i] = max(call_values[i], exchange_rates[i]-K);
        };
    return call_values[0];
};
    


int main()
{
    option call;
    
    double delta;
    double gamma;
    double theta;
    double vega;
    double rho;
    
    double price = 0;
    price = call.price(call.S,call.K,call.r,call.sigma,call.time);
    
    call.greeks(call.S,call.K,call.r,call.sigma,call.time,
                        delta,gamma,theta,vega,rho);
                        
                        
    
    cout<<price<<"\n";
    
    
    
    /*
    
    double fx_option_price;
    
    fx_option_price = call.currency_option(call.S,call.K,call.r,call.r_f,call.sigma,call.time,call.steps);
    
    cout<<fx_option_price;
    */
    
    
    
return 0;
}



