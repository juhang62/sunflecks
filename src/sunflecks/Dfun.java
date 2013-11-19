/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package sunflecks;

/**
 *
 * @author lfhan
 */
public class Dfun {
    double[] tau;
    double seq;
    double rmax;
    double tmax;
    double kr;
    double kt;
    double gamma;
    double psi;
    double vj;
    double vfeq;
    double vceq;
    double ci;
    double wc;
    boolean isgs;
    
    //constructor for stomata use
    public Dfun(double taugi, double taugd, double taupi, double tauw, double seq){
        this.tau=new double[4];
        this.tau[0]=taugi;
        this.tau[1]=taugd;
        this.tau[2]=taupi;
        this.tau[3]=tauw;
        this.seq=seq;
        this.isgs=true;
    }
    
    //constructor for photosynthesis use
    public Dfun(double[] tau, double rmax, double tmax, double kr,
            double kt, double gamma, double psi, double vj, double vfeq, 
            double vceq, double ci, double wc){
        this.tau=tau;
        this.rmax=rmax;
        this.tmax=tmax;
        this.kr=kr;
        this.kt=kt;
        this.gamma=gamma;
        this.psi=psi;
        this.vj=vj;
        this.vfeq=vfeq;
        this.vceq=vceq;
        this.ci=ci;
        this.wc=wc;
        this.isgs=false;     
    }
    
    public double[] getdf(double[] y0){
        double[] df;
        if (this.isgs==true){
            df=new double[3];
            if (this.seq>y0[0]){
                df[0]=(this.seq-y0[0])/this.tau[0];
            }else{
                df[0]=(this.seq-y0[0])/this.tau[1];
            }
            df[1]=(y0[0]-y0[1])/this.tau[2];
            df[2]=(y0[1]-y0[2])/this.tau[3]; 
        }else{
            df=new double[5];
            double tauc;
            double tauf;
            if (vceq>y0[0]){
                tauc=tau[0];
            }else{
                tauc=tau[1];
            }
            if (vfeq>y0[1]){
                tauf=tau[2];
            }else{
                tauf=tau[3];
            }
            
            df[0]=(vceq-y0[0])/tauc; //eqn4 vc
            df[1]=(vfeq-y0[1])/tauf; //eqn5 vf
            df[2]=vj*(1-y0[2]/tmax)-(y0[1]*y0[2]/(kt+y0[2]))*(1-y0[3]/rmax); //eqn1 T, no5/3 in units of RuBP
            df[3]=(y0[1]*y0[2]/(kt+y0[2]))*(1-y0[3]/rmax)-(1+2*gamma/ci)*wc*(y0[3]/(kr+y0[3])); //eqn2a R
            df[4]=(2*gamma/ci)*wc*(y0[3]/(kr+y0[3]))-psi*y0[4];//eqn3 G     
        }

        return df; 
    }
    
    public double[] dfphoto(double[] y0, double vj, double vceq, double vfeq,
            double ci,double wc){
        double[] df;
        df=new double[5];
        double tauc;
        double tauf;
        if (vceq>y0[0]){
             tauc=tau[0];
        }else{
             tauc=tau[1];
        }
        if (vfeq>y0[1]){
             tauf=tau[2];
        }else{
             tauf=tau[3];
        }
        
        df[0]=(vceq-y0[0])/tauc; //eqn4 vc
        df[1]=(vfeq-y0[1])/tauf; //eqn5 vf
        df[2]=vj*(1-y0[2]/tmax)-(y0[1]*y0[2]/(kt+y0[2]))*(1-y0[3]/rmax); //eqn1 T, no5/3 in units of RuBP
        df[3]=(y0[1]*y0[2]/(kt+y0[2]))*(1-y0[3]/rmax)-(1+2*gamma/ci)*wc*(y0[3]/(kr+y0[3])); //eqn2a R
        df[4]=(2*gamma/ci)*wc*(y0[3]/(kr+y0[3]))-psi*y0[4];//eqn3 G
        return df;
    
    }
    
}
