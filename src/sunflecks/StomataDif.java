/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package sunflecks;

/**
 *
 * @author lfhan
 */
public class StomataDif {
    double taugi;
    double taugd;
    double taupi;
    double tauw;
            
    public StomataDif(double taugi, double taugd, double taupi, double tauw){
        this.taugi=taugi;
        this.taugd=taugd;
        this.taupi=taupi;
        this.tauw=tauw;      
    }
    
    
    public double[] diffun(double[] y0, double seq){
        double[] df;
        df=new double[3];
        if (seq>y0[0]){
            df[1]=(seq-y0[0])/this.taugi;
        }else{
            df[1]=(seq-y0[0])/this.taugd;
        }
        df[2]=(y0[0]-y0[1])/this.taupi;
        df[3]=(y0[1]-y0[2])/this.tauw;
        return df;     
    }
}
