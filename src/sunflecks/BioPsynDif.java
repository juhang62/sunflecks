/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package sunflecks;

/**
 *
 * @author lfhan
 */
public class BioPsynDif {
    double tauci;
    double taucd;
    double taufi;
    double taufd;
    double rmax;
    double tmax;
    double kr;
    double kt;
    double gamma;
    double psi;
    
    public BioPsynDif(double[] tau, double rmax, double tmax, double kr,
            double kt, double gamma, double psi){
        this.tauci=tau[0];
        this.taucd=tau[1];
        this.taufi=tau[2];
        this.taufd=tau[3];
        this.rmax=rmax;
        this.tmax=tmax;
        this.kr=kr;
        this.kt=kt;
        this.gamma=gamma;
        this.psi=psi;
    }
    
    public double[] diffun(double[] y0, double vj, double vceq, double vfeq,
            double ci,double wc){
        double[] df;
        df=new double[5];
        double tauc;
        double tauf;
        if (vceq>y0[0]){
             tauc=tauci;
        }else{
             tauc=taucd;
        }
        if (vfeq>y0[1]){
             tauf=taufi;
        }else{
             tauf=taufd;
        }
        
        df[0]=(vceq-y0[0])/tauc; //eqn4 vc
        df[1]=(vfeq-y0[1])/tauf; //eqn5 vf
        df[2]=vj*(1-y0[2]/tmax)-(y0[1]*y0[2]/(kt+y0[2]))*(1-y0[3]/rmax); //eqn1 T, no5/3 in units of RuBP
        df[3]=(y0[1]*y0[2]/(kt+y0[2]))*(1-y0[3]/rmax)-(1+2*gamma/ci)*wc*(y0[3]/(kr+y0[3])); //eqn2a R
        df[4]=(2*gamma/ci)*wc*(y0[3]/(kr+y0[3]))-psi*y0[4];//eqn3 G
        return df;
    }
    
    
}
